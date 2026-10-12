// High-utility sequential pattern mining.
//
// Semantics (unchanged from 0.1.1):
//   An embedding of pattern P in sequence s is an increasing tuple of positions
//   p1 < p2 < ... < pL with items[pi] == P[i].
//   u(P, s) = MAX over embeddings of P in s of the sum of utils at those positions.
//   u(P)    = SUM over sequences of u(P, s).
//   support(P) = number of sequences in which P has at least one embedding.
//   Utilities are strictly positive.
//
// Search design:
//   * A node keeps, per sequence, a Pareto frontier of (end_pos, utility) states:
//     the best embedding of the current pattern ending at each position. Because a
//     state (q1,u1) dominates (q2,u2) whenever q1 <= q2 and u1 >= u2, the retained
//     frontier has end_pos and utility both strictly increasing. That makes the best
//     state with end_pos < p simply the last one before p, so a single left-to-right
//     walk scores EVERY candidate at once -- O(sequence length) per node instead of
//     O(#candidates * #states * length).
//   * g[p][b] = max utility of a subsequence of positions >= p with length <= b.
//     A state ending at q with b slots left cannot gain more than g[q+1][b], which is
//     strictly tighter than the plain suffix sum because it respects max_length.
//   * Top-K uses a bounded heap instead of re-sorting the result vector per insert.
//     The heap comparator is the same total order used for the final output
//     (utility desc, length asc, pattern lex asc), so "evict the worst" keeps exactly
//     the first K patterns of the final order. Deterministic, independent of visit order.

#include <Rcpp.h>
#include <algorithm>
#include <cmath>
#include <string>
#include <vector>
#include <future>
#include <mutex>

using namespace Rcpp;

namespace {

struct State {
  int end_pos;
  double utility;
};

struct Result {
  int offset;  // into the pattern arena
  int length;
  double utility;
  int support;
};

struct HeapItem {
  std::vector<int> pattern;
  double utility;
  int support;
};

// Total order of the final output: utility desc, length asc, pattern lex asc.
// Used both to sort the output and as the heap comparator. Note that passing this
// to std::push_heap makes the heap top the WORST element, which is what we evict.
inline bool ranks_before(const std::vector<int>& pa, double ua,
                         const std::vector<int>& pb, double ub) {
  if (ua != ub) return ua > ub;
  if (pa.size() != pb.size()) return pa.size() < pb.size();
  return std::lexicographical_compare(pa.begin(), pa.end(), pb.begin(), pb.end());
}

inline bool heap_cmp(const HeapItem& a, const HeapItem& b) {
  return ranks_before(a.pattern, a.utility, b.pattern, b.utility);
}

// Sparse, grouped, per-sequence state table. Groups are the non-empty sequences
// only, in ascending sequence order; inside a group both end_pos and utility are
// strictly increasing.
struct StateTable {
  std::vector<State> data;
  std::vector<int> gstart;  // ngroups + 1 boundaries into data
  std::vector<int> gseq;    // sequence id of each group
};

struct Builder {
  int item;  // the candidate item this builder belongs to
  std::vector<State> data;
  std::vector<int> gstart;
  std::vector<int> gseq;
};

}  // namespace

class HUSPM {
 public:
  HUSPM(const IntegerVector& seq, const IntegerVector& item,
        const NumericVector& utility, int n_sequences, int n_items,
        double min_utility, int top_k_, int max_length, int min_support_,
        bool include_singletons_)
      : nseq(n_sequences),
        nitems(n_items),
        min_util(min_utility),
        max_len(max_length),
        min_support(min_support_),
        top_k(top_k_ > 0 ? top_k_ : 0),
        use_topk(top_k_ > 0),
        include_singletons(include_singletons_) {
    items.resize(nseq);
    utils.resize(nseq);
    for (int i = 0; i < seq.size(); ++i) {
      items[seq[i]].push_back(item[i]);
      utils[seq[i]].push_back(utility[i]);
    }
    build_bounds();

    seq_seen.assign(nitems, 0);
    seq_best_util.assign(nitems, 0.0);
    seq_best_ub.assign(nitems, 0.0);
    cand_seen.assign(nitems, 0);
    cand_util.assign(nitems, 0.0);
    cand_ub.assign(nitems, 0.0);
    cand_supp.assign(nitems, 0);
    surv_index.assign(nitems, -1);
    surv_best.assign(nitems, -1.0);
  }

  void run() {
    // Occurrence lists. Built by ascending (sequence, position), so each list is
    // already sorted by sequence and then by position.
    std::vector<std::vector<std::pair<int, int> > > occ(nitems);
    for (int s = 0; s < nseq; ++s) {
      const int n = (int)items[s].size();
      for (int p = 0; p < n; ++p) occ[items[s][p]].push_back(std::make_pair(s, p));
    }

    const int root_budget = max_len - 1;

    // Collect items that need processing (non-empty occurrence list and pass threshold)
    std::vector<std::pair<int, StateTable>> root_jobs;
    std::vector<double> root_util;
    std::vector<int> root_support;
    for (int x = 0; x < nitems; ++x) {
      if (occ[x].empty()) continue;

      StateTable st;
      double util_sum = 0.0;
      int support = 0;
      size_t i = 0;
      while (i < occ[x].size()) {
        const int s = occ[x][i].first;
        st.gstart.push_back((int)st.data.size());
        st.gseq.push_back(s);
        // Pareto frontier of a singleton = running max of its own utilities.
        double best = -1.0;
        size_t k = i;
        while (k < occ[x].size() && occ[x][k].first == s) {
          const int p = occ[x][k].second;
          const double v = utils[s][p];
          if (v > best) {
            best = v;
            st.data.push_back(State{p, v});
          }
          ++k;
        }
        util_sum += best;
        ++support;
        i = k;
      }
      st.gstart.push_back((int)st.data.size());

      // Check if we need to process this root at all
      bool need_process = false;
      if (include_singletons || max_len > 1) {
        double ub = util_sum;
        if (max_len > 1) {
          ub += table_ub(st, root_budget);
        }
        if (ub + eps() >= threshold()) {
          need_process = true;
        }
      }
      if (need_process || include_singletons) {
        root_jobs.emplace_back(x, std::move(st));
        root_util.push_back(util_sum);
        root_support.push_back(support);
      }
    }

    // If we only have a few root jobs, process serially
    if (root_jobs.empty()) {
      return;
    }

    const size_t n_jobs = root_jobs.size();
    if (n_jobs == 1) {
      // Single job - process serially
      const int x = root_jobs[0].first;
      StateTable& st = root_jobs[0].second;
      const double util_sum = root_util[0];
      const int support = root_support[0];
      path.assign(1, x);
      if (include_singletons) emit(util_sum, support);
      if (max_len > 1 && table_ub(st, root_budget) + eps() >= threshold()) node(st);
      path.clear();
      return;
    }

    // Parallel processing of root jobs: each root job gets its own thread because
    // processing is completely independent and there's no shared mutable state.
    std::vector<std::future<HUSPM>> futures;
    futures.reserve(n_jobs);

    for (size_t j = 0; j < n_jobs; ++j) {
      const int x = root_jobs[j].first;
      StateTable st = std::move(root_jobs[j].second);
      const double util_sum = root_util[j];
      const int support = root_support[j];

      // Copy the original instance immediately, before any merging has happened,
      // so we get a pristine copy with empty results/heap. If we don't do this now,
      // a concurrent merge could happen before the async lambda runs, and we'd get
      // duplicate items from earlier processing.
      HUSPM worker_original = *this;

      futures.push_back(std::async(std::launch::async, [worker_original, x, st = std::move(st), util_sum, support, root_budget]() mutable -> HUSPM {
        // We already have a pristine copy of the original HUSPM instance so this
        // thread has its own copy. All parameters are copied from the original, and
        // we just need to replace the scratch buffers with newly allocated ones that are
        // thread-local.

        // We need our own copies of scratch buffers because they get modified during
        // node processing. Each thread works on its own roots independently, so
        // they don't interfere.
        std::vector<int> local_path;
        local_path.assign(1, x);
        worker_original.path.swap(local_path);

        // Reset the scratch buffers to be empty/cleared for this thread's use
        worker_original.seq_seen.assign(worker_original.nitems, 0);
        worker_original.seq_best_util.assign(worker_original.nitems, 0.0);
        worker_original.seq_best_ub.assign(worker_original.nitems, 0.0);
        worker_original.cand_seen.assign(worker_original.nitems, 0);
        worker_original.cand_util.assign(worker_original.nitems, 0.0);
        worker_original.cand_ub.assign(worker_original.nitems, 0.0);
        worker_original.cand_supp.assign(worker_original.nitems, 0);
        worker_original.cand_list.clear();
        worker_original.survivors.clear();
        worker_original.surv_best.assign(worker_original.nitems, -1.0);
        worker_original.surv_index.assign(worker_original.nitems, -1);
        worker_original.surv_touched.clear();

        // Emit singleton if requested
        if (worker_original.include_singletons) {
          if (support >= worker_original.min_support && util_sum + 1e-9 * std::max(1.0, std::fabs(worker_original.threshold())) >= worker_original.threshold()) {
            if (worker_original.use_topk) {
              HeapItem hi;
              hi.pattern = worker_original.path;
              hi.utility = util_sum;
              hi.support = support;
              worker_original.heap.push_back(std::move(hi));
              std::push_heap(worker_original.heap.begin(), worker_original.heap.end(), heap_cmp);
              if ((int)worker_original.heap.size() > worker_original.top_k) {
                std::pop_heap(worker_original.heap.begin(), worker_original.heap.end(), heap_cmp);
                worker_original.heap.pop_back();
              }
            } else {
              const int off = (int)worker_original.arena.size();
              worker_original.arena.insert(worker_original.arena.end(), worker_original.path.begin(), worker_original.path.end());
              worker_original.results.push_back(Result{off, (int)worker_original.path.size(), util_sum, support});
            }
          }
        }

        // Recurse if there's room and we pass the bound
        if (worker_original.max_len > 1 && worker_original.table_ub(st, root_budget) + worker_original.eps() >= worker_original.threshold()) {
          worker_original.node(st);
        }

        return worker_original;
      }));
    }

    // Collect results from all threads and merge into final output
    std::mutex results_mutex;
    for (auto& future : futures) {
      HUSPM worker = future.get();
      std::lock_guard<std::mutex> lock(results_mutex);
      if (use_topk) {
        // Merge heap items into the global heap
        for (auto& item : worker.heap) {
          heap.push_back(std::move(item));
          std::push_heap(heap.begin(), heap.end(), heap_cmp);
          if ((int)heap.size() > top_k) {
            std::pop_heap(heap.begin(), heap.end(), heap_cmp);
            heap.pop_back();
          }
        }
      } else {
        // Append results and arena
        size_t offset = arena.size();
        for (auto& res : worker.results) {
          // The pattern offsets in worker results are relative to worker.arena. We need
          // to adjust them to point into the global arena.
          Result shifted = res;
          shifted.offset += offset;
          results.push_back(shifted);
        }
        arena.insert(arena.end(), worker.arena.begin(), worker.arena.end());
        // Check for excessive memory usage in exhaustive mode
        const size_t max_arena_bytes = 500000000; // 500MB
        if (arena.size() * sizeof(int) > max_arena_bytes) {
          stop("Excessive memory usage in exhaustive mining: "
               "pattern arena exceeds 500MB. "
               "Use the 'top_k' parameter to limit the number of output patterns.");
        }
      }
    }
  }

  DataFrame to_dataframe() {
    std::vector<std::string> pats;
    std::vector<int> lens, sups;
    std::vector<double> us;

    if (use_topk) {
      std::sort(heap.begin(), heap.end(), heap_cmp);
      for (size_t i = 0; i < heap.size(); ++i) {
        pats.push_back(join(heap[i].pattern));
        lens.push_back((int)heap[i].pattern.size());
        us.push_back(heap[i].utility);
        sups.push_back(heap[i].support);
      }
    } else {
      const int* arena_p = arena.empty() ? NULL : &arena[0];
      std::sort(results.begin(), results.end(),
                [arena_p](const Result& a, const Result& b) {
                  if (a.utility != b.utility) return a.utility > b.utility;
                  if (a.length != b.length) return a.length < b.length;
                  return std::lexicographical_compare(
                      arena_p + a.offset, arena_p + a.offset + a.length,
                      arena_p + b.offset, arena_p + b.offset + b.length);
                });
      for (size_t i = 0; i < results.size(); ++i) {
        pats.push_back(join(results[i].offset, results[i].length));
        lens.push_back(results[i].length);
        us.push_back(results[i].utility);
        sups.push_back(results[i].support);
      }
    }

    const int m = (int)pats.size();
    CharacterVector pattern_id(m);
    IntegerVector length(m), support(m);
    NumericVector util(m);
    for (int i = 0; i < m; ++i) {
      pattern_id[i] = pats[i];
      length[i] = lens[i];
      util[i] = us[i];
      support[i] = sups[i];
    }
    return DataFrame::create(_["pattern_id"] = pattern_id, _["length"] = length,
                             _["utility"] = util, _["support"] = support,
                             _["stringsAsFactors"] = false);
  }

 private:
  // ---- parameters -----------------------------------------------------------
  int nseq, nitems;
  double min_util;
  int max_len, min_support, top_k;
  bool use_topk, include_singletons;

  // ---- input ----------------------------------------------------------------
  std::vector<std::vector<int> > items;
  std::vector<std::vector<double> > utils;

  // ---- bounds ---------------------------------------------------------------
  std::vector<std::vector<double> > gtab;  // gtab[s][b*stride + p], stride = len+1
  int bmax = 0;

  // ---- reusable scratch (allocated once, never per node) --------------------
  std::vector<int> seq_seen, seq_touched;
  std::vector<double> seq_best_util, seq_best_ub;
  std::vector<int> cand_seen, cand_supp, cand_list, survivors;
  std::vector<double> cand_util, cand_ub;
  std::vector<int> surv_index;
  std::vector<double> surv_best;
  std::vector<int> surv_touched;
  std::vector<int> path;

  // ---- output ---------------------------------------------------------------
  std::vector<Result> results;
  std::vector<int> arena;
  std::vector<HeapItem> heap;

  void build_bounds() {
    int max_seq_len = 0;
    for (int s = 0; s < nseq; ++s)
      max_seq_len = std::max(max_seq_len, (int)items[s].size());
    // g[p][b] saturates once b >= len - p <= max_seq_len, so the budget axis only
    // needs min(max_len, max_seq_len) + 1 columns. Without this clamp a large
    // max_length on a long database allocates O(max_length * events).
    bmax = std::min(max_len, max_seq_len);
    gtab.resize(nseq);
    // Check total allocation size to avoid out-of-memory
    // 1GB (1e9 bytes) is a reasonable upper bound
    const size_t max_total_bytes = 1000000000;
    size_t total_bytes = 0;
    for (int s = 0; s < nseq; ++s) {
      const int n = (int)items[s].size();
      const size_t bytes = sizeof(double) * (size_t)(bmax + 1) * (n + 1);
      if (bytes > max_total_bytes) {
        stop("A single sequence is too long for the given max_length. "
             "Total allocation for bound table exceeds 1GB.");
      }
      total_bytes += bytes;
      if (total_bytes > max_total_bytes) {
        stop("Total allocation for bound table exceeds 1GB. "
             "Try reducing max_length or filtering shorter sequences.");
      }
    }
    // Allocate now that we've checked the size
    for (int s = 0; s < nseq; ++s) {
      const int n = (int)items[s].size();
      const int stride = n + 1;
      std::vector<double>& G = gtab[s];
      G.assign((size_t)(bmax + 1) * stride, 0.0);
      for (int p = n - 1; p >= 0; --p) {
        for (int b = 1; b <= bmax; ++b) {
          const double skip = G[(size_t)b * stride + (p + 1)];
          const double take = utils[s][p] + G[(size_t)(b - 1) * stride + (p + 1)];
          G[(size_t)b * stride + p] = skip > take ? skip : take;
        }
      }
    }
  }

  // Max utility obtainable from positions >= pos using at most `budget` items.
  inline double bound_at(int s, int pos, int budget) const {
    const int stride = (int)items[s].size() + 1;
    const int b = budget < bmax ? budget : bmax;
    return gtab[s][(size_t)b * stride + pos];
  }

  // ---- thresholds -----------------------------------------------------------
  double threshold() const {
    if (!use_topk) return min_util;
    if ((int)heap.size() < top_k) return min_util;
    return std::max(min_util, heap.front().utility);
  }

  // Relative slack: an absolute 1e-12 stops covering an ulp once utilities exceed
  // ~1e4, and a bound that rounds one ulp below a real descendant would drop it.
  double eps() const {
    const double t = threshold();
    return 1e-9 * std::max(1.0, std::fabs(t));
  }

  void emit(double utility, int support) {
    if (support < min_support) return;  // never let these raise the Top-K threshold
    const double t = threshold();
    if (utility + 1e-9 * std::max(1.0, std::fabs(t)) < t) return;

    if (use_topk) {
      HeapItem hi;
      hi.pattern = path;
      hi.utility = utility;
      hi.support = support;
      heap.push_back(std::move(hi));
      std::push_heap(heap.begin(), heap.end(), heap_cmp);
      if ((int)heap.size() > top_k) {
        std::pop_heap(heap.begin(), heap.end(), heap_cmp);
        heap.pop_back();
      }
    } else {
      const int off = (int)arena.size();
      arena.insert(arena.end(), path.begin(), path.end());
      results.push_back(Result{off, (int)path.size(), utility, support});
      // Check for excessive memory usage in exhaustive mode
      const size_t max_arena_bytes = 500000000; // 500MB
      if (arena.size() * sizeof(int) > max_arena_bytes) {
        stop("Excessive memory usage in exhaustive mining: "
             "pattern arena exceeds 500MB. "
             "Use the 'top_k' parameter to limit the number of output patterns.");
      }
    }
  }

  double table_ub(const StateTable& st, int budget) const {
    double ub = 0.0;
    for (size_t g = 0; g < st.gseq.size(); ++g) {
      const int s = st.gseq[g];
      double best = 0.0;
      for (int i = st.gstart[g]; i < st.gstart[g + 1]; ++i) {
        const double v =
            st.data[i].utility + bound_at(s, st.data[i].end_pos + 1, budget);
        if (v > best) best = v;
      }
      ub += best;
    }
    return ub;
  }

  void node(const StateTable& st) {
    checkUserInterrupt();
    const int L = (int)path.size();
    if (L >= max_len) return;
    const int budget = max_len - L - 1;  // slots left AFTER appending one item

    // ---- one walk per sequence scores every candidate at once ----------------
    cand_list.clear();
    for (size_t g = 0; g < st.gseq.size(); ++g) {
      const int s = st.gseq[g];
      const std::vector<int>& it = items[s];
      const std::vector<double>& ut = utils[s];
      const int n = (int)it.size();
      const int gb = st.gstart[g], ge = st.gstart[g + 1];

      seq_touched.clear();
      int ptr = gb;
      for (int p = 0; p < n; ++p) {
        // best state with end_pos < p; valid because the frontier is sorted and
        // its utilities increase with end_pos.
        while (ptr < ge && st.data[ptr].end_pos < p) ++ptr;
        if (ptr == gb) continue;

        const int x = it[p];
        const double val = st.data[ptr - 1].utility + ut[p];
        const double vub = val + bound_at(s, p + 1, budget);
        if (seq_seen[x] == 0) {
          seq_seen[x] = 1;
          seq_best_util[x] = val;
          seq_best_ub[x] = vub;
          seq_touched.push_back(x);
        } else {
          // val and vub are maximised independently: they need not come from the
          // same position p, and each is a valid per-sequence max on its own.
          if (val > seq_best_util[x]) seq_best_util[x] = val;
          if (vub > seq_best_ub[x]) seq_best_ub[x] = vub;
        }
      }
      // accumulate in ascending sequence order, matching the legacy summation order
      for (size_t t = 0; t < seq_touched.size(); ++t) {
        const int x = seq_touched[t];
        if (cand_seen[x] == 0) {
          cand_seen[x] = 1;
          cand_list.push_back(x);
        }
        cand_util[x] += seq_best_util[x];
        cand_ub[x] += seq_best_ub[x];
        cand_supp[x] += 1;
        seq_seen[x] = 0;
      }
    }
    if (cand_list.empty()) return;

    std::sort(cand_list.begin(), cand_list.end(), [this](int a, int b) {
      if (cand_ub[a] != cand_ub[b]) return cand_ub[a] > cand_ub[b];
      return a < b;
    });

    survivors.clear();
    for (size_t i = 0; i < cand_list.size(); ++i) {
      const int x = cand_list[i];
      // Emission is gated on the pattern's own utility, recursion on its upper
      // bound. These are different predicates and must not be merged.
      path.push_back(x);
      emit(cand_util[x], cand_supp[x]);
      path.pop_back();

      // Support is monotone non-increasing down the tree, so a child below
      // min_support cannot have a qualifying descendant.
      if (cand_supp[x] < min_support) continue;
      if (cand_ub[x] + eps() >= threshold()) survivors.push_back(x);
    }

    // Release the accumulators before recursing; children reuse the same arrays.
    for (size_t i = 0; i < cand_list.size(); ++i) {
      const int x = cand_list[i];
      cand_seen[x] = 0;
      cand_util[x] = 0.0;
      cand_ub[x] = 0.0;
      cand_supp[x] = 0;
    }
    cand_list.clear();

    // Clear the whole map before recursing, not one entry per iteration: a child
    // node reads surv_index while materializing its own states, and would
    // otherwise pick up our entries for items it does not track. Clearing all of
    // them up front also keeps the invariant that surv_index is entirely -1 on
    // entry to (and exit from) every node.
    for (int j = 0; j < nitems; ++j) surv_index[j] = -1;

    int ns = (int)survivors.size();
    std::vector<Builder> bd(ns);  // local: `survivors` itself is not safe to read below
    surv_best.assign(ns, -1.0);
    for (int j = 0; j < ns; ++j) {
      bd[j].item = survivors[j];
      surv_index[survivors[j]] = j;
    }

    for (size_t g = 0; g < st.gseq.size(); ++g) {
      const int s = st.gseq[g];
      const std::vector<int>& it = items[s];
      const std::vector<double>& ut = utils[s];
      const int n = (int)it.size();
      const int gb = st.gstart[g], ge = st.gstart[g + 1];

      surv_touched.clear();
      int ptr = gb;
      for (int p = 0; p < n; ++p) {
        while (ptr < ge && st.data[ptr].end_pos < p) ++ptr;
        if (ptr == gb) continue;

        const int j = surv_index[it[p]];
        if (j < 0) continue;
        const double val = st.data[ptr - 1].utility + ut[p];
        if (surv_best[j] < 0.0) {
          // first touch in this sequence: open a group
          surv_best[j] = val;
          surv_touched.push_back(j);
          bd[j].gstart.push_back((int)bd[j].data.size());
          bd[j].gseq.push_back(s);
          bd[j].data.push_back(State{p, val});
        } else if (val > surv_best[j]) {
          // running max in ascending position order == Pareto frontier
          surv_best[j] = val;
          bd[j].data.push_back(State{p, val});
        }
      }
      for (size_t t = 0; t < surv_touched.size(); ++t)
        surv_best[surv_touched[t]] = -1.0;
    }

    // Clear the whole map before recursing, not one entry per iteration: a child
    // node reads surv_index while materializing its own states, and would
    // otherwise pick up our entries for items it does not track. Clearing all of
    // them up front also keeps the invariant that surv_index is entirely -1 on
    // entry to (and exit from) every node.
    for (int j = 0; j < ns; ++j) surv_index[bd[j].item] = -1;

    for (int j = 0; j < ns; ++j) {
      bd[j].gstart.push_back((int)bd[j].data.size());

      StateTable child;
      child.data.swap(bd[j].data);
      child.gstart.swap(bd[j].gstart);
      child.gseq.swap(bd[j].gseq);

      // Read the item from the (local) builder. `survivors` is a member, and the
      // recursive call below overwrites it with the child's own survivor list --
      // reading survivors[j] here after the first recursion pairs the wrong item
      // id with the wrong state table.
      path.push_back(bd[j].item);
      node(child);
      path.pop_back();
    }
  }

  std::string join(const std::vector<int>& v) const {
    std::string s;
    for (size_t i = 0; i < v.size(); ++i) {
      if (i) s += ',';
      s += std::to_string(v[i]);
    }
    return s;
  }

  std::string join(int off, int len) const {
    std::string s;
    for (int i = 0; i < len; ++i) {
      if (i) s += ',';
      s += std::to_string(arena[off + i]);
    }
    return s;
  }
};

// [[Rcpp::export]]
DataFrame huspm_cpp(IntegerVector seq, IntegerVector item, NumericVector utility,
                    int n_sequences, int n_items, double min_utility,
                    int top_k = 0, int max_length = 5, int min_support = 1,
                    bool include_singletons = true) {
  if (seq.size() != item.size() || seq.size() != utility.size())
    stop("seq, item, utility must have equal length.");
  if (n_sequences < 1 || n_items < 1)
    stop("n_sequences and n_items must be positive.");
  if (max_length < 1) stop("max_length must be >= 1.");
  if (min_support < 1) stop("min_support must be >= 1.");
  if (top_k < 0) stop("top_k must be >= 0.");
  if (!R_finite(min_utility) || min_utility < 0)
    stop("min_utility must be finite and non-negative.");

  for (int i = 0; i < seq.size(); ++i) {
    if (seq[i] < 0 || seq[i] >= n_sequences) stop("Invalid sequence id.");
    if (item[i] < 0 || item[i] >= n_items) stop("Invalid item id.");
    if (!R_finite(utility[i]) || utility[i] <= 0)
      stop("Utilities must be finite and strictly positive.");
  }

  HUSPM miner(seq, item, utility, n_sequences, n_items, min_utility, top_k,
              max_length, min_support, include_singletons);
  miner.run();
  return miner.to_dataframe();
}
