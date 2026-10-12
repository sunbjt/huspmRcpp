# huspmRcpp: 高效用频繁序列模式挖掘 R 包

<!-- badges: start -->
<!-- badges: end -->

`huspmRcpp` 是一个基于 Rcpp 实现的高效用频繁序列模式挖掘 (High-Utility Sequential Pattern Mining, HUSPM) 算法包。

## 算法特性

- 实现了高效的深度优先搜索算法，结合边界剪枝快速剪去不可能产生有效模式的分支
- 自动并行处理：当有多个根项目（不同的初始项目）时，每个根项目会在单独的线程中并行处理，自动利用多个 CPU 核心加速计算，对用户完全透明
- 支持限制输出模式数量 (`top_k` 参数)，只输出满足条件中效用最高的前 K 个模式
- 支持限制模式的最大长度 (`max_length` 参数)
- 支持最小支持度阈值 (`min_support` 参数)，只保留出现在至少一定数量序列中的模式
- 可以选择是否输出单项目模式

## 安装

使用 `pak` 安装：

```r
# install.packages("pak")
pak::pkg_install("sunbjt/huspmRcpp")
```

或者使用 `devtools` 安装：

```r
# install.packages("devtools")
devtools::install_github("sunbjt/huspmRcpp")
```

## 使用示例

```r
library(huspmRcpp)

# 构造一个简单的例子
x <- data.frame(
  sequence_id = c(1, 1, 1, 2, 2, 2, 3, 3),
  item_id     = c("A", "B", "C", "A", "B", "D", "A", "C"),
  utility     = c(  5,   2,  10,   4,   3,   8,   6,   7)
)

# 挖掘所有满足最小效用阈值 >= 15 的模式
ans <- huspm(x, min_utility = 15, max_length = 3, min_support = 1)
print(ans)
```

输出结果:
```
<huspm>  5  patterns
     pattern length utility support utility_ratio
      A -> C      2      28       2     0.6222222
           C      1      17       2     0.3777778
 A -> B -> C      3      17       1     0.3777778
           A      1      15       3     0.3333333
 A -> B -> D      3      15       1     0.3333333
```

## 参数说明

### `huspm()` 主要参数

| 参数 | 说明 | 默认值 |
|------|------|--------|
| `data` | 输入数据框，必须包含 `sequence_id`, `item_id`, `utility` 三列 | - |
| `min_utility` | 最小总模式效用，所有满足条件的模式都会被输出 | `0` |
| `top_k` | 如果提供，则只保留效用最高的前 K 个模式。`NULL` 表示输出所有满足条件的模式 | `NULL` |
| `max_length` | 模式中允许的最大项目个数，必须是正整数 | `5` |
| `min_support` | 一个模式至少需要出现在多少个不同序列中才能被保留 | `1` |
| `include_singletons` | 是否输出长度为 1 的单项目模式 | `TRUE` |
| `sequence_id` | 输入数据框中序列 ID 列的列名 | `"sequence_id"` |
| `item_id` | 输入数据框中项目 ID 列的列名 | `"item_id"` |
| `utility` | 输入数据框中项目效用列的列名 | `"utility"` |
| `time` | 可选列名，用于排序同一个序列中的事件顺序，如果不提供则按输入顺序排序 | `NULL` |

## 输出说明

输出是一个 `data.frame`，包含以下列：

| 列名 | 说明 |
|------|------|
| `pattern` | 模式标签，项目之间用 ` -> ` 分隔 |
| `length` | 模式中包含的项目个数 |
| `utility` | 模式的总效用（每个序列中最大嵌入效用之和） |
| `support` | 包含这个模式的序列个数 |
| `utility_ratio` | 模式总效用除以所有输入事件效用之和。注意：由于模式效用是每个序列中的最大嵌入，分母是所有输入事件的总和，因此这个比例不会加起来等于 1.0，并且可能会对出现在多个模式中的事件重复计数。 |

## 并行计算说明

`huspmRcpp` 使用 C++11 标准库中的 `std::async` 和 `std::future` 实现自动并行。每个根项目（即长度为 1 的单项目模式）在单独的线程中处理，每个线程都有自己独立的 scratch 空间，处理过程中不需要锁，只有在最后合并结果的时候才使用一个互斥锁保护共享的结果结构。

当输入数据中有很多独特的根项目时，并行计算可以显著提升性能，算法会自动利用多个 CPU 核心，不需要用户做任何修改。

## 数据格式要求

- 输入数据是长格式（long format），每一行代表一个序列中的一个项目事件
- 每个项目在同一个序列中可以出现多次
- 所有项目的效用必须是正的有限数值
- `sequence_id` 可以是任何类型（整数、字符、因子都可以），算法内部会自动转换为整数索引

## 参考文献

本算法基于以下论文中的思想实现：

- **High-Utility Sequential Pattern Mining: A Brief Survey**
- 算法使用了动态规划预处理边界，以及 Pareto  frontier 优化状态空间，保证了算法的效率。

## 许可证

MIT License
