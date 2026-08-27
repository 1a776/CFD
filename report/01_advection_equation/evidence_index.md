# 第一题报告证据索引

本索引用于说明报告中的关键结论分别来自哪个题目文件、配置、代码、日志、数据或图片。

## 1. 参数源与问题定义

| 内容 | 证据 | 用途 |
|---|---|---|
| 原始第一题 | `../../pdf/01/first_advection_problem.pdf` | 控制方程、正弦波算例、固体旋转算例和题目要求 |
| 有限体积推导 | `../../pdf/01/advection_fvm_explicit_solver_derivation.pdf` | 控制体积分、面通量、迎风格式、显式 Euler 和 CFL 推导 |
| 四边形一阶迎风配置 | `../../scripts/configs/01_advection_equation/01_sine_wave_quad_upwind.json` | 网格、速度、终止时间、CFL 和空间格式 |
| 四边形线性迎风配置 | `../../scripts/configs/01_advection_equation/02_sine_wave_quad_linearUpwind.json` | 高阶格式扩展 |
| 三角形一阶迎风配置 | `../../scripts/configs/01_advection_equation/03_sine_wave_tri_upwind.json` | 三角形网格正弦波案例 |
| 三角形线性迎风配置 | `../../scripts/configs/01_advection_equation/03_sine_wave_tri_linearUpwind.json` | 三角形高阶格式扩展 |
| 四边形旋转配置 | `../../scripts/configs/01_advection_equation/04_solid_rotation_quad_upwind.json` | 四边形固体旋转多网格研究 |
| 三角形旋转配置 | `../../scripts/configs/01_advection_equation/04_solid_rotation_tri_upwind.json` | 三角形固体旋转多网格研究 |

## 2. 开发与数值实现

| 内容 | 证据 | 用途 |
|---|---|---|
| 学生版求解器 | `../../UDF/solver/01_advection_equation/explicitAdvectionFoamStudent/explicitAdvectionFoamStudent.C` | 面通量、CFL、显式残差和时间推进 |
| UDF 开发说明 | `../../UDF/README.md` | 公式到 OpenFOAM 代码的对应关系 |
| 构建脚本 | `../../scripts/build_student_solver.sh` | 编译复现 |
| 求解器可执行文件 | `../../build/01_advection_equation/bin/explicitAdvectionFoamStudent` | 编译结果 |
| 四边形网格定义 | `../../cases/01_advection_equation/01_sine_wave_quad/N20/system/blockMeshDict` | 由 JSON 直接生成的 `blockMesh` 网格定义 |
| 三角形网格生成器 | `../../scripts/common/gmsh_tri_mesh.py` | Gmsh 三角形棱柱网格生成 |
| 统一案例脚本 | `../../scripts/common/foam_case.py` | 根据 JSON 准备 OpenFOAM case |

## 3. 正弦波结果

| 内容 | 证据 |
|---|---|
| 四边形汇总数据 | `../../data/01_advection_equation/analysis/01_sine_wave_quad/convergence_summary.csv` |
| 四边形分析文字 | `../../data/01_advection_equation/analysis/01_sine_wave_quad/analysis.md` |
| 四边形误差图 | `../../figures/01_advection_equation/analysis/01_sine_wave_quad/convergence_errors.png` |
| 四边形收敛阶图 | `../../figures/01_advection_equation/analysis/01_sine_wave_quad/convergence_order.png` |
| 四边形全分辨率场对比 | `../../figures/01_advection_equation/analysis/01_sine_wave_quad/all_N_comparison.png` |
| 四边形一阶迎风 N80 诊断图 | `../../figures/01_advection_equation/cases/01_sine_wave_quad/N80/` |
| 四边形线性迎风汇总图 | `../../figures/01_advection_equation/analysis/02_sine_wave_quad_linearUpwind/` |
| 四边形线性迎风 N80 诊断图 | `../../figures/01_advection_equation/cases/02_sine_wave_quad_linearUpwind/N80/` |
| 三角形汇总数据 | `../../data/01_advection_equation/analysis/03_sine_wave_tri_upwind/convergence_summary.csv` |
| 三角形分析文字 | `../../data/01_advection_equation/analysis/03_sine_wave_tri_upwind/analysis.md` |
| 三角形误差图 | `../../figures/01_advection_equation/analysis/03_sine_wave_tri_upwind/convergence_errors.png` |
| 三角形收敛阶图 | `../../figures/01_advection_equation/analysis/03_sine_wave_tri_upwind/convergence_order.png` |
| 三角形一阶迎风全分辨率场对比 | `../../figures/01_advection_equation/analysis/03_sine_wave_tri_upwind/all_N_comparison.png` |
| 三角形一阶迎风 N80 诊断图 | `../../figures/01_advection_equation/cases/03_sine_wave_tri_upwind/N80/` |
| 三角形线性迎风汇总图 | `../../figures/01_advection_equation/analysis/03_sine_wave_tri_linearUpwind/` |
| 三角形线性迎风 N80 诊断图 | `../../figures/01_advection_equation/cases/03_sine_wave_tri_linearUpwind/N80/` |

每个分辨率的详细数据还包括：

```text
data/01_advection_equation/cases/<caseName>/Nxx/summary.json
data/01_advection_equation/cases/<caseName>/Nxx/time_history.csv
data/01_advection_equation/cases/<caseName>/Nxx/field_data.csv
data/01_advection_equation/cases/<caseName>/Nxx/error_field.csv
```

## 4. 固体旋转结果

| 内容 | 证据 |
|---|---|
| 四边形旋转汇总 | `../../data/01_advection_equation/cases/04_solid_rotation_quad_upwind/N50/summary.json` |
| 四边形旋转汇总 | `../../data/01_advection_equation/cases/04_solid_rotation_quad_upwind/N100/summary.json` |
| 四边形旋转汇总 | `../../data/01_advection_equation/cases/04_solid_rotation_quad_upwind/N200/summary.json` |
| 四边形旋转图 | `../../figures/01_advection_equation/cases/04_solid_rotation_quad_upwind/N100/field_comparison.png` |
| 四边形旋转 N50/N100/N200 全部图 | `../../figures/01_advection_equation/cases/04_solid_rotation_quad_upwind/` |
| 三角形旋转汇总 | `../../data/01_advection_equation/cases/04_solid_rotation_tri_upwind/N50/summary.json` |
| 三角形旋转汇总 | `../../data/01_advection_equation/cases/04_solid_rotation_tri_upwind/N100/summary.json` |
| 三角形旋转汇总 | `../../data/01_advection_equation/cases/04_solid_rotation_tri_upwind/N200/summary.json` |
| 三角形旋转图 | `../../figures/01_advection_equation/cases/04_solid_rotation_tri_upwind/N100/field_comparison.png` |
| 三角形旋转 N50/N100/N200 全部图 | `../../figures/01_advection_equation/cases/04_solid_rotation_tri_upwind/` |

## 5. 错误、修复与限制

| 内容 | 证据 |
|---|---|
| 运行 Bug 记录 | `../../docs/bug_log.md` |
| 时间精度问题修复 | `../../UDF/solver/01_advection_equation/explicitAdvectionFoamStudent/explicitAdvectionFoamStudent.C` |

## 6. 证据使用规则

报告中的数值结论应优先追溯到 `summary.json` 或 `convergence_summary.csv`；
图像结论应追溯到 `figures/`；数学公式应追溯到题目和有限体积推导 PDF；
程序行为应追溯到求解器源码、OpenFOAM case 字典和运行日志。
