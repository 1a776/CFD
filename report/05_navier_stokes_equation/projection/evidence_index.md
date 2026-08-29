# 第五题压力投影法报告证据索引

本索引用于说明压力投影法报告中的关键结论来自哪个题面、配置、代码、数据、图件或日志。

## 1. 题面与统一定义

| 内容 | 证据 | 用途 |
|---|---|---|
| 原始训练题集 | `../../../pdf/training_examples_incomp.pdf` | Navier-Stokes 方程、投影法要求和腔流算例 |
| 统一问题定义 | `../../../0-caseDict/caseDict` | 第五题的物理问题、边界、参考数据和局限性 |

## 2. 算法实现

| 内容 | 证据 | 用途 |
|---|---|---|
| 投影法源码 | `../../../UDF/solver/05_navier_stokes_equation/projectionFoamStudent/projectionFoamStudent.C` | 预测速度、压力泊松方程、通量和速度修正 |
| 构建产物 | `../../../build/05_navier_stokes_equation/bin/projectionFoamStudent` | 求解器构建结果 |
| 方腔后处理 | `../../../scripts/postprocess_lid_cavity.py` | Ghia 中心线采样、误差和场图 |
| 三角腔后处理 | `../../../scripts/postprocess_triangular_cavity.py` | 主涡、流函数、中心线和涡结构 |
| 方腔混合网格生成器 | `../../../scripts/common/gmsh_hybrid_cavity_mesh.py` | 方腔混合网格和边界物理组 |
| 三角腔混合网格生成器 | `../../../scripts/common/gmsh_triangular_cavity_hybrid_mesh.py` | 等边三角腔几何、壁面层和物理组 |

## 3. 方腔配置与结果

| 内容 | 证据 |
|---|---|
| Re1000 hybrid40 配置 | `../../../scripts/configs/05_navier_stokes_equation/07_lid_driven_cavity_projection_Re1000_hybrid40.json` |
| Re1000 hybrid80 配置 | `../../../scripts/configs/05_navier_stokes_equation/08_lid_driven_cavity_projection_Re1000_hybrid80.json` |
| Re3200 hybrid80 配置 | `../../../scripts/configs/05_navier_stokes_equation/11_lid_driven_cavity_projection_Re3200_hybrid80.json` |
| Re1000 hybrid40 摘要 | `../../../data/05_navier_stokes_equation/cases/07_lid_driven_cavity_projection_Re1000_hybrid40/summary.json` |
| Re1000 hybrid80 摘要 | `../../../data/05_navier_stokes_equation/cases/08_lid_driven_cavity_projection_Re1000_hybrid80/summary.json` |
| Re3200 hybrid80 摘要 | `../../../data/05_navier_stokes_equation/cases/11_lid_driven_cavity_projection_Re3200_hybrid80/summary.json` |
| Re1000 hybrid40 场图 | `../../../figures/05_navier_stokes_equation/cases/07_lid_driven_cavity_projection_Re1000_hybrid40/field_and_streamlines.png` |
| Re1000 hybrid80 场图 | `../../../figures/05_navier_stokes_equation/cases/08_lid_driven_cavity_projection_Re1000_hybrid80/field_and_streamlines.png` |
| Re1000 hybrid80 中心线图 | `../../../figures/05_navier_stokes_equation/cases/08_lid_driven_cavity_projection_Re1000_hybrid80/centerline_comparison.png` |
| Re1000 hybrid40 日志 | `../../../cases/05_navier_stokes_equation/07_lid_driven_cavity_projection_Re1000_hybrid40/log.projectionFoamStudent` |
| Re1000 hybrid80 日志 | `../../../cases/05_navier_stokes_equation/08_lid_driven_cavity_projection_Re1000_hybrid80/log.projectionFoamStudent` |
| Re3200 hybrid80 网格检查 | `../../../cases/05_navier_stokes_equation/11_lid_driven_cavity_projection_Re3200_hybrid80/log.checkMesh` |

每个方腔摘要旁边的中心线数据包括：

```text
u_centerline.csv
v_centerline.csv
summary.json
```

## 4. 三角腔配置与结果

| 内容 | 证据 |
|---|---|
| Re100 配置 | `../../../scripts/configs/05_navier_stokes_equation/26_triangular_cavity_projection_Re100_hybrid80.json` |
| Re200 配置 | `../../../scripts/configs/05_navier_stokes_equation/27_triangular_cavity_projection_Re200_hybrid80.json` |
| Re500 配置 | `../../../scripts/configs/05_navier_stokes_equation/28_triangular_cavity_projection_Re500_hybrid80.json` |
| Re100 摘要 | `../../../data/05_navier_stokes_equation/cases/26_triangular_cavity_projection_Re100_hybrid80/summary.json` |
| Re200 摘要 | `../../../data/05_navier_stokes_equation/cases/27_triangular_cavity_projection_Re200_hybrid80/summary.json` |
| Re500 摘要 | `../../../data/05_navier_stokes_equation/cases/28_triangular_cavity_projection_Re500_hybrid80/summary.json` |
| Re200 场图 | `../../../figures/05_navier_stokes_equation/cases/27_triangular_cavity_projection_Re200_hybrid80/field_streamlines.png` |
| Re200 流函数图 | `../../../figures/05_navier_stokes_equation/cases/27_triangular_cavity_projection_Re200_hybrid80/streamfunction_vortices.png` |
| Re200 网格检查 | `../../../cases/05_navier_stokes_equation/27_triangular_cavity_projection_Re200_hybrid80/N80/log.checkMesh` |
| Re500 网格检查 | `../../../cases/05_navier_stokes_equation/28_triangular_cavity_projection_Re500_hybrid80/N80/log.checkMesh` |

三角腔结果数据包括：

```text
u_centerline.csv
v_horizontal.csv
summary.json
```

## 5. 已知证据缺口

| 问题 | 现状 | 报告处理 |
|---|---|---|
| Re100 三角腔 | 摘要终止于 `t≈0.6`，没有当前 case 目录 | 只作为瞬态历史诊断，不作稳态结论 |
| Re200/Re500 三角腔 | 保留网格和 `0` 时刻字段，缺少最终正时间目录与 solver log | 只作为摘要级结果证据 |
| Re3200 方腔 | 有摘要和中心线数据，当前未找到投影 solver log | 不把日志级稳态结论扩大到该案例 |
| 完整矩阵 | 配置文件多于当前数据摘要 | 明确报告范围为已归档子集 |

数值结论优先来自 `summary.json`；
中心线曲线来自 `u_centerline.csv` 和 `v_centerline.csv`；
图像结论来自 `figures/`；
稳态性结论必须同时参考 case-local solver log。
