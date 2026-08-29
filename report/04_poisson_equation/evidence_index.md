# 第四题报告证据索引

本索引用于说明第四题报告中的关键结论分别来自哪个题面、配置、代码、日志、数据或图片。

## 1. 题面与问题定义

| 内容 | 证据 | 用途 |
|---|---|---|
| 原始训练题集 | `../../pdf/training_examples_incomp.pdf` | 第四题 Poisson 方程题面、制造解、源项和边界条件 |
| 统一问题定义 | `../../0-caseDict/caseDict` | 题意、假设、验收标准、报告输出路径 |

## 2. 数值配置与实现

| 内容 | 证据 | 用途 |
|---|---|---|
| 四边形配置 | `../../scripts/configs/04_poisson_equation/01_poisson_manufactured_quad.json` | 域、网格类型、线性求解器、容差、边界条件 |
| 三角形配置 | `../../scripts/configs/04_poisson_equation/02_poisson_manufactured_tri.json` | 域、网格类型、线性求解器、容差、边界条件 |
| 求解器源码 | `../../UDF/solver/04_poisson_equation/poissonFoamStudent/poissonFoamStudent.C` | 拉普拉斯组装、源项读取和求解流程 |
| 求解器可执行文件 | `../../build/04_poisson_equation/bin/poissonFoamStudent` | 编译结果 |

## 3. 四边形网格结果

| 内容 | 证据 | 用途 |
|---|---|---|
| 汇总分析文字 | `../../data/04_poisson_equation/analysis/01_poisson_manufactured_quad/analysis.md` | 误差表和收敛阶结论 |
| 原始汇总数据 | `../../data/04_poisson_equation/analysis/01_poisson_manufactured_quad/convergence_summary.csv` | 可复核的汇总表 |
| 原始逐项结果 | `../../data/04_poisson_equation/analysis/01_poisson_manufactured_quad/raw_results.csv` | 误差、范数和附加统计量 |
| 汇总图 | `../../figures/04_poisson_equation/analysis/01_poisson_manufactured_quad/all_N_comparison.png` | 不同分辨率的场对比 |
| 误差图 | `../../figures/04_poisson_equation/analysis/01_poisson_manufactured_quad/convergence_errors.png` | 误差下降趋势 |
| 收敛阶图 | `../../figures/04_poisson_equation/analysis/01_poisson_manufactured_quad/convergence_order.png` | 观察收敛阶 |
| N10 单案例图 | `../../figures/04_poisson_equation/cases/01_poisson_manufactured_quad/N10/field_comparison.png` | 最粗网格场对比 |
| N20 单案例图 | `../../figures/04_poisson_equation/cases/01_poisson_manufactured_quad/N20/field_comparison.png` | 次粗网格场对比 |
| N40 单案例图 | `../../figures/04_poisson_equation/cases/01_poisson_manufactured_quad/N40/field_comparison.png` | 中等网格场对比 |
| N80 单案例图 | `../../figures/04_poisson_equation/cases/01_poisson_manufactured_quad/N80/field_comparison.png` | 最细网格场对比 |

每个网格分辨率的详细文件为：

```text
data/04_poisson_equation/cases/01_poisson_manufactured_quad/Nxx/summary.json
data/04_poisson_equation/cases/01_poisson_manufactured_quad/Nxx/time_history.csv
data/04_poisson_equation/cases/01_poisson_manufactured_quad/Nxx/field_data.csv
data/04_poisson_equation/cases/01_poisson_manufactured_quad/Nxx/error_field.csv
```

## 4. 三角形网格结果

| 内容 | 证据 | 用途 |
|---|---|---|
| 汇总分析文字 | `../../data/04_poisson_equation/analysis/02_poisson_manufactured_tri/analysis.md` | 误差表和收敛阶结论 |
| 原始汇总数据 | `../../data/04_poisson_equation/analysis/02_poisson_manufactured_tri/convergence_summary.csv` | 可复核的汇总表 |
| 原始逐项结果 | `../../data/04_poisson_equation/analysis/02_poisson_manufactured_tri/raw_results.csv` | 误差、范数和附加统计量 |
| 汇总图 | `../../figures/04_poisson_equation/analysis/02_poisson_manufactured_tri/all_N_comparison.png` | 不同分辨率的场对比 |
| 误差图 | `../../figures/04_poisson_equation/analysis/02_poisson_manufactured_tri/convergence_errors.png` | 误差下降趋势 |
| 收敛阶图 | `../../figures/04_poisson_equation/analysis/02_poisson_manufactured_tri/convergence_order.png` | 观察收敛阶 |
| N10 单案例图 | `../../figures/04_poisson_equation/cases/02_poisson_manufactured_tri/N10/field_comparison.png` | 最粗网格场对比 |
| N20 单案例图 | `../../figures/04_poisson_equation/cases/02_poisson_manufactured_tri/N20/field_comparison.png` | 次粗网格场对比 |
| N40 单案例图 | `../../figures/04_poisson_equation/cases/02_poisson_manufactured_tri/N40/field_comparison.png` | 中等网格场对比 |
| N80 单案例图 | `../../figures/04_poisson_equation/cases/02_poisson_manufactured_tri/N80/field_comparison.png` | 最细网格场对比 |

每个网格分辨率的详细文件为：

```text
data/04_poisson_equation/cases/02_poisson_manufactured_tri/Nxx/summary.json
data/04_poisson_equation/cases/02_poisson_manufactured_tri/Nxx/time_history.csv
data/04_poisson_equation/cases/02_poisson_manufactured_tri/Nxx/field_data.csv
data/04_poisson_equation/cases/02_poisson_manufactured_tri/Nxx/error_field.csv
```

## 5. 日志与完整性

| 内容 | 证据 | 用途 |
|---|---|---|
| 四边形求解日志 | `../../cases/04_poisson_equation/01_poisson_manufactured_quad/N80/log.poissonFoamStudent` | 求解结束、残差与运行轨迹 |
| 三角形求解日志 | `../../cases/04_poisson_equation/02_poisson_manufactured_tri/N80/log.poissonFoamStudent` | 求解结束、残差与运行轨迹 |
| 四边形网格检查 | `../../cases/04_poisson_equation/01_poisson_manufactured_quad/N80/log.checkMesh` | `meshOK=true` 的来源 |
| 三角形网格检查 | `../../cases/04_poisson_equation/02_poisson_manufactured_tri/N80/log.checkMesh` | `meshOK=true` 的来源 |

## 6. 证据使用规则

报告中的数值结论优先追溯到 `summary.json`、`convergence_summary.csv` 或 `analysis.md`；
图像结论追溯到 `figures/`；
程序行为追溯到求解器源码、OpenFOAM case 字典和运行日志。
