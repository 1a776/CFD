# 第二题报告证据索引

本索引用于说明第二题报告中的关键结论分别来自哪个题目文件、配置、代码、日志、数据或图片。

## 1. 参数源与问题定义

| 内容 | 证据 | 用途 |
|---|---|---|
| 原始题目 | `../../pdf/training_examples_incomp.pdf` | 第二题原题 |
| 自包含题面 | `../../pdf/tex/第二题_二维扩散方程_自包含题目.tex` | 便于统一符号和题意 |
| 题目推导 | `../../pdf/tex/diffusion_fvm_explicit_solver_derivation.tex` | 有限体积离散推导 |
| 四边形配置 | `../../scripts/configs/02_diffusion_equation/01_discontinuous_quad.json` | 四边形算例入口 |
| 三角形配置 | `../../scripts/configs/02_diffusion_equation/02_discontinuous_tri.json` | 三角形算例入口 |

## 2. 开发与数值实现

| 内容 | 证据 | 用途 |
|---|---|---|
| 学生版求解器 | `../../UDF/solver/02_diffusion_equation/explicitDiffusionFoamStudent/explicitDiffusionFoamStudent.C` | 显式扩散项、时间推进、日志输出 |
| UDF 开发说明 | `../../UDF/README.md` | 公式到 OpenFOAM 代码的对应关系 |
| 构建脚本 | `../../scripts/build_student_solver.sh` | 编译复现 |
| 求解器可执行文件 | `../../build/02_diffusion_equation/bin/explicitDiffusionFoamStudent` | 编译结果 |
| 四边形模板 | `../../cases/02_diffusion_equation/01_discontinuous_quad/N20/system/blockMeshDict` | `blockMesh` 网格定义 |
| 三角形网格生成器 | `../../scripts/common/gmsh_tri_mesh.py` | Gmsh 三角形棱柱网格生成 |
| 统一案例脚本 | `../../scripts/common/foam_case.py` | 根据 JSON 准备、运行、后处理 case |

## 3. 四边形结果

| 内容 | 证据 |
|---|---|
| 汇总数据 | `../../data/02_diffusion_equation/analysis/01_discontinuous_quad/convergence_summary.csv` |
| 分析文字 | `../../data/02_diffusion_equation/analysis/01_discontinuous_quad/analysis.md` |
| 误差图 | `../../figures/02_diffusion_equation/analysis/01_discontinuous_quad/convergence_errors.png` |
| 收敛阶图 | `../../figures/02_diffusion_equation/analysis/01_discontinuous_quad/convergence_order.png` |
| 全分辨率场对比 | `../../figures/02_diffusion_equation/analysis/01_discontinuous_quad/all_N_comparison.png` |
| N10 诊断图 | `../../figures/02_diffusion_equation/cases/01_discontinuous_quad/N10/` |
| N20 诊断图 | `../../figures/02_diffusion_equation/cases/01_discontinuous_quad/N20/` |
| N40 诊断图 | `../../figures/02_diffusion_equation/cases/01_discontinuous_quad/N40/` |
| N80 诊断图 | `../../figures/02_diffusion_equation/cases/01_discontinuous_quad/N80/` |

## 4. 三角形结果

| 内容 | 证据 |
|---|---|
| 汇总数据 | `../../data/02_diffusion_equation/analysis/02_discontinuous_tri/convergence_summary.csv` |
| 分析文字 | `../../data/02_diffusion_equation/analysis/02_discontinuous_tri/analysis.md` |
| 误差图 | `../../figures/02_diffusion_equation/analysis/02_discontinuous_tri/convergence_errors.png` |
| 收敛阶图 | `../../figures/02_diffusion_equation/analysis/02_discontinuous_tri/convergence_order.png` |
| 全分辨率场对比 | `../../figures/02_diffusion_equation/analysis/02_discontinuous_tri/all_N_comparison.png` |
| N10 诊断图 | `../../figures/02_diffusion_equation/cases/02_discontinuous_tri/N10/` |
| N20 诊断图 | `../../figures/02_diffusion_equation/cases/02_discontinuous_tri/N20/` |
| N40 诊断图 | `../../figures/02_diffusion_equation/cases/02_discontinuous_tri/N40/` |
| N80 诊断图 | `../../figures/02_diffusion_equation/cases/02_discontinuous_tri/N80/` |

## 5. 错误、修复与限制

| 内容 | 证据 |
|---|---|
| 扩散算例中的边界和初值修复 | `../../scripts/common/foam_case.py` |
| 三角形网格域支持修复 | `../../scripts/common/gmsh_tri_mesh.py` |
| 三角形扩散后处理修复 | `../../scripts/common/postprocess_case.py` |
| 运行日志 | `../../cases/02_diffusion_equation/<case>/Nxx/log.*` |

## 6. 证据使用规则

报告中的数值结论优先追溯到 `summary.json` 或 `convergence_summary.csv`；
图像结论优先追溯到 `figures/`；数学公式优先追溯到题目与推导 PDF；
程序行为优先追溯到求解器源码、OpenFOAM case 字典和运行日志。

