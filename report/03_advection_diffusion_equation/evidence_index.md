# 第三题报告证据索引

本索引用于说明第三题报告中的关键结论分别来自哪个题目文件、配置、代码、日志、数据或图片。

## 1. 参数源与问题定义

| 内容 | 证据 | 用途 |
|---|---|---|
| 第三题题面与解答 | `../../pdf/题目解答.pdf` | 控制方程、两个算例族、网格与时间要求 |
| 正弦波四边形配置 | `../../scripts/configs/03_advection_diffusion_equation/01_sine_wave_quad_upwind.json` | 周期边界、四边形网格、误差定义、终止时间 |
| 正弦波三角形配置 | `../../scripts/configs/03_advection_diffusion_equation/02_sine_wave_tri_upwind.json` | 周期边界、三角形网格、误差定义、终止时间 |
| 旋转尖峰四边形近似边界配置 | `../../scripts/configs/03_advection_diffusion_equation/03_rotating_peak_quad_upwind.json` | 零 Dirichlet 近似边界、四边形网格 |
| 旋转尖峰三角形近似边界配置 | `../../scripts/configs/03_advection_diffusion_equation/04_rotating_peak_tri_upwind.json` | 零 Dirichlet 近似边界、三角形网格 |
| 旋转尖峰四边形解析边界配置 | `../../scripts/configs/03_advection_diffusion_equation/05_rotating_peak_quad_analyticDirichlet_upwind.json` | 解析 Dirichlet 边界对照 |
| 旋转尖峰三角形解析边界配置 | `../../scripts/configs/03_advection_diffusion_equation/06_rotating_peak_tri_analyticDirichlet_upwind.json` | 解析 Dirichlet 边界对照 |
| caseDict 主记录 | `../../0-caseDict/caseDict` | 统一参数源、报表口径、工作流状态 |

## 2. 数值实现与复现入口

| 内容 | 证据 | 用途 |
|---|---|---|
| 学生版求解器源码 | `../../UDF/solver/03_advection_diffusion_equation/explicitAdvectionDiffusionFoamStudent/explicitAdvectionDiffusionFoamStudent.C` | 对流项、扩散项、显式时间推进 |
| 编译脚本 | `../../scripts/build_student_solver.sh` | 复现求解器编译 |
| 编译产物 | `../../build/03_advection_diffusion_equation/bin/explicitAdvectionDiffusionFoamStudent` | 编译结果 |
| 批量运行脚本 | `../../scripts/run_study.py` | 根据配置运行多分辨率研究 |
| 汇总分析脚本 | `../../scripts/analyze_study.py` | 生成收敛表 |
| 绘图脚本 | `../../scripts/plot_study.py` | 生成收敛图和场图 |
| 通用工具 | `../../scripts/common/advection_diffusion_tools.py` `../../scripts/common/advection_sine.py` `../../scripts/common/advection_rotation.py` `../../scripts/common/study_analysis.py` | 正弦波、旋转尖峰、误差定义和结果汇总 |

## 3. 正弦波结果

| 内容 | 证据 | 用途 |
|---|---|---|
| 四边形汇总数据 | `../../data/03_advection_diffusion_equation/analysis/01_sine_wave_quad_upwind/convergence_summary.csv` | 主误差、旧误差、收敛阶 |
| 四边形分析文字 | `../../data/03_advection_diffusion_equation/analysis/01_sine_wave_quad_upwind/analysis.md` | 对病态相对误差的解释 |
| 四边形误差图 | `../../figures/03_advection_diffusion_equation/analysis/01_sine_wave_quad_upwind/convergence_errors.png` | 误差随网格变化 |
| 四边形收敛阶图 | `../../figures/03_advection_diffusion_equation/analysis/01_sine_wave_quad_upwind/convergence_order.png` | 观察收敛阶 |
| 四边形全分辨率场图 | `../../figures/03_advection_diffusion_equation/analysis/01_sine_wave_quad_upwind/all_N_comparison.png` | 各分辨率场对比 |
| 四边形 N80 诊断图 | `../../figures/03_advection_diffusion_equation/cases/01_sine_wave_quad_upwind/N80/` | 最细网格场与时间历史 |
| 三角形汇总数据 | `../../data/03_advection_diffusion_equation/analysis/02_sine_wave_tri_upwind/convergence_summary.csv` | 主误差、旧误差、收敛阶 |
| 三角形分析文字 | `../../data/03_advection_diffusion_equation/analysis/02_sine_wave_tri_upwind/analysis.md` | 对病态相对误差的解释 |
| 三角形误差图 | `../../figures/03_advection_diffusion_equation/analysis/02_sine_wave_tri_upwind/convergence_errors.png` | 误差随网格变化 |
| 三角形收敛阶图 | `../../figures/03_advection_diffusion_equation/analysis/02_sine_wave_tri_upwind/convergence_order.png` | 观察收敛阶 |
| 三角形全分辨率场图 | `../../figures/03_advection_diffusion_equation/analysis/02_sine_wave_tri_upwind/all_N_comparison.png` | 各分辨率场对比 |
| 三角形 N80 诊断图 | `../../figures/03_advection_diffusion_equation/cases/02_sine_wave_tri_upwind/N80/` | 最细网格场与时间历史 |
| 病态误差归档 | `../../data/03_advection_diffusion_equation/pathological_relative_error/README.md` | 旧 exact-relative 误差的历史记录 |

## 4. 旋转尖峰结果

| 内容 | 证据 | 用途 |
|---|---|---|
| 四边形近似边界汇总 | `../../data/03_advection_diffusion_equation/analysis/03_rotating_peak_quad_upwind/convergence_summary.csv` | 零边界近似下的误差与收敛阶 |
| 四边形近似边界分析 | `../../data/03_advection_diffusion_equation/analysis/03_rotating_peak_quad_upwind/analysis.md` | 边界解释与结论 |
| 四边形近似边界图 | `../../figures/03_advection_diffusion_equation/analysis/03_rotating_peak_quad_upwind/all_N_comparison.png` | 各分辨率场对比 |
| 四边形 N80 诊断图 | `../../figures/03_advection_diffusion_equation/cases/03_rotating_peak_quad_upwind/N80/` | 最细网格场、轮廓和剖面 |
| 三角形近似边界汇总 | `../../data/03_advection_diffusion_equation/analysis/04_rotating_peak_tri_upwind/convergence_summary.csv` | 零边界近似下的误差与收敛阶 |
| 三角形近似边界分析 | `../../data/03_advection_diffusion_equation/analysis/04_rotating_peak_tri_upwind/analysis.md` | 边界解释与结论 |
| 三角形近似边界图 | `../../figures/03_advection_diffusion_equation/analysis/04_rotating_peak_tri_upwind/all_N_comparison.png` | 各分辨率场对比 |
| 三角形 N80 诊断图 | `../../figures/03_advection_diffusion_equation/cases/04_rotating_peak_tri_upwind/N80/` | 最细网格场、轮廓和剖面 |
| 四边形解析边界汇总 | `../../data/03_advection_diffusion_equation/analysis/05_rotating_peak_quad_analyticDirichlet_upwind/convergence_summary.csv` | 解析边界下的误差与收敛阶 |
| 四边形解析边界分析 | `../../data/03_advection_diffusion_equation/analysis/05_rotating_peak_quad_analyticDirichlet_upwind/analysis.md` | 与近似边界的对照 |
| 四边形解析边界图 | `../../figures/03_advection_diffusion_equation/analysis/05_rotating_peak_quad_analyticDirichlet_upwind/all_N_comparison.png` | 各分辨率场对比 |
| 四边形 N80 诊断图 | `../../figures/03_advection_diffusion_equation/cases/05_rotating_peak_quad_analyticDirichlet_upwind/N80/` | 最细网格场、轮廓和剖面 |
| 三角形解析边界汇总 | `../../data/03_advection_diffusion_equation/analysis/06_rotating_peak_tri_analyticDirichlet_upwind/convergence_summary.csv` | 解析边界下的误差与收敛阶 |
| 三角形解析边界分析 | `../../data/03_advection_diffusion_equation/analysis/06_rotating_peak_tri_analyticDirichlet_upwind/analysis.md` | 与近似边界的对照 |
| 三角形解析边界图 | `../../figures/03_advection_diffusion_equation/analysis/06_rotating_peak_tri_analyticDirichlet_upwind/all_N_comparison.png` | 各分辨率场对比 |
| 三角形 N80 诊断图 | `../../figures/03_advection_diffusion_equation/cases/06_rotating_peak_tri_analyticDirichlet_upwind/N80/` | 最细网格场、轮廓和剖面 |

## 5. 完整性与限制

| 内容 | 证据 | 用途 |
|---|---|---|
| 代表性状态位 | `../../data/03_advection_diffusion_equation/cases/*/N80/summary.json` | `meshOK`、`solverEnded`、`solverFatal`、`finalTimeError` |
| 工作流记录 | `../../0-caseDict/caseDict` | `workflowState`、`reportData`、`provenance.updatedAt` |
| 报告正文 | `./report.md` | 结果叙述和结论 |
| 病态误差限制 | `../../data/03_advection_diffusion_equation/pathological_relative_error/` | 旧误差为何不能作为主指标 |

## 6. 证据使用规则

报告中的数值结论优先追溯到 `summary.json` 和 `convergence_summary.csv`；  
图像结论优先追溯到 `figures/`；  
误差定义和边界处理优先追溯到 `caseDict`、配置文件和对应分析文字。  
