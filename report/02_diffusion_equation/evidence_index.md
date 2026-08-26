# 第二题报告证据索引

本索引用于说明 `report.md` 中的关键结论分别来自哪个题目文件、推导文档、配置、代码、
OpenFOAM case、日志、数据或图片。报告中的数值结论应优先追溯到 `summary.json`、
`convergence_summary.csv` 或后处理 CSV；图像结论应追溯到 `figures/`；程序行为应追溯到
求解器源码、OpenFOAM 字典和运行日志。

## 1. 参数源与问题定义

| 内容 | 证据 | 用途 |
|---|---|---|
| 原始第二题 | `../../pdf/training_examples_incomp.pdf` | 控制方程、两个扩散算例、网格和误差要求 |
| 第二题自包含题面 | `../../pdf/tex/第二题_二维扩散方程_自包含题目.tex` | 中文化题面和问题整理 |
| 有限体积推导 | `../../pdf/tex/diffusion_fvm_explicit_solver_derivation.tex` | 控制体积分、扩散通量、显式 Euler 和扩散时间步推导 |
| 间断四边形配置 | `../../scripts/configs/02_diffusion_equation/01_discontinuous_quad.json` | 区域、初值、Neumann 边界、网格、终止时间和扩散格式 |
| 间断三角形配置 | `../../scripts/configs/02_diffusion_equation/02_discontinuous_tri.json` | 三角形网格下的同一间断扩散算例 |
| Gaussian 四边形配置 | `../../scripts/configs/02_diffusion_equation/03_gaussian_quad.json` | Gaussian 初值、解析解 Dirichlet 边界、四边形网格 |
| Gaussian 三角形配置 | `../../scripts/configs/02_diffusion_equation/04_gaussian_tri.json` | Gaussian 初值、解析解 Dirichlet 边界、三角形网格 |

## 2. 开发与数值实现

| 内容 | 证据 | 用途 |
|---|---|---|
| 学生版扩散求解器 | `../../UDF/solver/02_diffusion_equation/explicitDiffusionFoamStudent/` | 扩散残差、显式时间推进、边界条件更新和运行输出 |
| 求解器可执行文件 | `../../build/02_diffusion_equation/bin/explicitDiffusionFoamStudent` | 编译结果和运行入口 |
| 构建脚本 | `../../scripts/build_student_solver.sh` | 从源码复现编译 |
| 统一案例生成脚本 | `../../scripts/common/foam_case.py` | 根据 JSON 生成 OpenFOAM case、Allrun、字典和边界 |
| 扩散解析解与初值工具 | `../../scripts/common/diffusion_tools.py` | 间断初值、Gaussian 初值、解析解和 codedFixedValue 边界 |
| 后处理脚本 | `../../scripts/common/postprocess_case.py` | 误差计算、CSV 输出和诊断图生成 |
| 研究汇总脚本 | `../../scripts/common/study_analysis.py` | 间断算例的结果收集、收敛阶计算和收敛图 |
| 三角形网格生成器 | `../../scripts/common/gmsh_tri_mesh.py` | Gmsh 三角形棱柱网格生成 |

## 3. OpenFOAM Case 与运行日志

| 内容 | 证据 | 用途 |
|---|---|---|
| 间断四边形 case | `../../cases/02_diffusion_equation/01_discontinuous_quad/` | 四边形 `blockMesh`、`zeroGradient` 边界和运行日志 |
| 间断三角形 case | `../../cases/02_diffusion_equation/02_discontinuous_tri/` | Gmsh 三角形网格、`createPatch` 和运行日志 |
| Gaussian 四边形 case | `../../cases/02_diffusion_equation/03_gaussian_quad/` | 四边形 Gaussian 初值和 `codedFixedValue` 边界 |
| Gaussian 三角形 case | `../../cases/02_diffusion_equation/04_gaussian_tri/` | 三角形 Gaussian 初值、普通外边界和 `codedFixedValue` 边界 |
| 四边形网格日志 | `../../cases/02_diffusion_equation/*_quad/N*/log.blockMesh` | `blockMesh` 网格生成过程 |
| 三角形网格转换日志 | `../../cases/02_diffusion_equation/*_tri/N*/log.gmshToFoam` | Gmsh 网格转 OpenFOAM 网格过程 |
| 三角形边界重建日志 | `../../cases/02_diffusion_equation/*_tri/N*/log.createPatch` | `zMin/zMax` empty 边界和外边界 patch 构造 |
| 网格检查日志 | `../../cases/02_diffusion_equation/*/N*/log.checkMesh` | 网格质量和 `Mesh OK` |
| 求解器日志 | `../../cases/02_diffusion_equation/*/N*/log.explicitDiffusionFoamStudent` | 时间推进、扩散 Co、残差、质量和场范围 |
| 单案例元数据 | `../../cases/02_diffusion_equation/*/N*/metadata.json` | 配置、分辨率、边界条件、格式和输出字段 |

## 4. 间断初值扩散结果

| 内容 | 证据 |
|---|---|
| 四边形汇总数据 | `../../data/02_diffusion_equation/analysis/01_discontinuous_quad/convergence_summary.csv` |
| 四边形分析文字 | `../../data/02_diffusion_equation/analysis/01_discontinuous_quad/analysis.md` |
| 四边形原始汇总 | `../../data/02_diffusion_equation/analysis/01_discontinuous_quad/raw_results.csv` |
| 四边形误差图 | `../../figures/02_diffusion_equation/analysis/01_discontinuous_quad/convergence_errors.png` |
| 四边形收敛阶图 | `../../figures/02_diffusion_equation/analysis/01_discontinuous_quad/convergence_order.png` |
| 四边形全分辨率场对比 | `../../figures/02_diffusion_equation/analysis/01_discontinuous_quad/all_N_comparison.png` |
| 四边形 N80 场图 | `../../figures/02_diffusion_equation/cases/01_discontinuous_quad/N80/field_comparison.png` |
| 四边形 N80 中线剖面 | `../../figures/02_diffusion_equation/cases/01_discontinuous_quad/N80/midline_profile.png` |
| 三角形汇总数据 | `../../data/02_diffusion_equation/analysis/02_discontinuous_tri/convergence_summary.csv` |
| 三角形分析文字 | `../../data/02_diffusion_equation/analysis/02_discontinuous_tri/analysis.md` |
| 三角形原始汇总 | `../../data/02_diffusion_equation/analysis/02_discontinuous_tri/raw_results.csv` |
| 三角形误差图 | `../../figures/02_diffusion_equation/analysis/02_discontinuous_tri/convergence_errors.png` |
| 三角形收敛阶图 | `../../figures/02_diffusion_equation/analysis/02_discontinuous_tri/convergence_order.png` |
| 三角形全分辨率场对比 | `../../figures/02_diffusion_equation/analysis/02_discontinuous_tri/all_N_comparison.png` |
| 三角形 N80 场图 | `../../figures/02_diffusion_equation/cases/02_discontinuous_tri/N80/field_comparison.png` |
| 三角形 N80 对角线剖面 | `../../figures/02_diffusion_equation/cases/02_discontinuous_tri/N80/diagonal_profile.png` |

每个分辨率的详细数据还包括：

```text
data/02_diffusion_equation/cases/01_discontinuous_quad/Nxx/summary.json
data/02_diffusion_equation/cases/01_discontinuous_quad/Nxx/time_history.csv
data/02_diffusion_equation/cases/01_discontinuous_quad/Nxx/field_data.csv
data/02_diffusion_equation/cases/01_discontinuous_quad/Nxx/error_field.csv

data/02_diffusion_equation/cases/02_discontinuous_tri/Nxx/summary.json
data/02_diffusion_equation/cases/02_discontinuous_tri/Nxx/time_history.csv
data/02_diffusion_equation/cases/02_discontinuous_tri/Nxx/field_data.csv
data/02_diffusion_equation/cases/02_discontinuous_tri/Nxx/error_field.csv
```

## 5. Gaussian 扩散结果

| 内容 | 证据 |
|---|---|
| 四边形 N10 结果 | `../../data/02_diffusion_equation/cases/03_gaussian_quad/N10/summary.json` |
| 四边形 N20 结果 | `../../data/02_diffusion_equation/cases/03_gaussian_quad/N20/summary.json` |
| 四边形 N40 结果 | `../../data/02_diffusion_equation/cases/03_gaussian_quad/N40/summary.json` |
| 四边形 N80 结果 | `../../data/02_diffusion_equation/cases/03_gaussian_quad/N80/summary.json` |
| 四边形 N80 场图 | `../../figures/02_diffusion_equation/cases/03_gaussian_quad/N80/field_comparison.png` |
| 四边形 N80 中线剖面 | `../../figures/02_diffusion_equation/cases/03_gaussian_quad/N80/midline_profile.png` |
| 四边形 N80 时间历史 | `../../figures/02_diffusion_equation/cases/03_gaussian_quad/N80/diffusion_step_history.png` |
| 三角形 N10 结果 | `../../data/02_diffusion_equation/cases/04_gaussian_tri/N10/summary.json` |
| 三角形 N20 结果 | `../../data/02_diffusion_equation/cases/04_gaussian_tri/N20/summary.json` |
| 三角形 N40 结果 | `../../data/02_diffusion_equation/cases/04_gaussian_tri/N40/summary.json` |
| 三角形 N80 结果 | `../../data/02_diffusion_equation/cases/04_gaussian_tri/N80/summary.json` |
| 三角形 N80 场图 | `../../figures/02_diffusion_equation/cases/04_gaussian_tri/N80/field_comparison.png` |
| 三角形 N80 对角线剖面 | `../../figures/02_diffusion_equation/cases/04_gaussian_tri/N80/diagonal_profile.png` |
| 三角形 N80 振幅历史 | `../../figures/02_diffusion_equation/cases/04_gaussian_tri/N80/amplitude_history.png` |

Gaussian 算例目前没有单独生成 `data/02_diffusion_equation/analysis/03_gaussian_quad/` 和
`data/02_diffusion_equation/analysis/04_gaussian_tri/` 下的跨分辨率汇总图；报告中的
Gaussian 收敛阶表由各分辨率 `summary.json` 的 $L_1$ 误差计算得到。

每个分辨率的详细数据还包括：

```text
data/02_diffusion_equation/cases/03_gaussian_quad/Nxx/summary.json
data/02_diffusion_equation/cases/03_gaussian_quad/Nxx/time_history.csv
data/02_diffusion_equation/cases/03_gaussian_quad/Nxx/field_data.csv
data/02_diffusion_equation/cases/03_gaussian_quad/Nxx/error_field.csv

data/02_diffusion_equation/cases/04_gaussian_tri/Nxx/summary.json
data/02_diffusion_equation/cases/04_gaussian_tri/Nxx/time_history.csv
data/02_diffusion_equation/cases/04_gaussian_tri/Nxx/field_data.csv
data/02_diffusion_equation/cases/04_gaussian_tri/Nxx/error_field.csv
```

## 6. 错误、修复与限制

| 内容 | 证据 | 说明 |
|---|---|---|
| Gaussian 三角形 `createPatch` 修复 | `../../scripts/common/foam_case.py` | 只有 `periodicXY` 使用 cyclic；Gaussian Dirichlet 外边界改为普通 patch |
| 修复后的三角形日志 | `../../cases/02_diffusion_equation/04_gaussian_tri/N10/log.createPatch` | `createPatch` 能完成普通外边界重建 |
| 运行 Bug 记录 | `../../docs/bug_log.md` | 项目运行中遇到的问题和修复记录 |
| Gaussian 汇总图限制 | `../../data/02_diffusion_equation/cases/03_gaussian_quad/`、`../../data/02_diffusion_equation/cases/04_gaussian_tri/` | 已有逐分辨率数据和图，尚未生成独立 analysis 汇总图 |

## 7. 证据使用规则

1. 数值误差、收敛阶、质量误差和时间步数应优先追溯到 `summary.json` 或
   `convergence_summary.csv`。
2. 具体时程、最大值变化和稳定性应追溯到 `time_history.csv` 和求解器日志。
3. 场形状、剖面和图像现象应追溯到 `figures/02_diffusion_equation/`。
4. 数学公式应追溯到原始题目和有限体积推导文档。
5. 程序行为应追溯到求解器源码、统一脚本、OpenFOAM case 字典和运行日志。
