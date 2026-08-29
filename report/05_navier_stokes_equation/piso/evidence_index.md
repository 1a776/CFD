# 第五题 PISO 法报告证据索引

本索引用于说明 PISO 报告中的关键结论来自哪个题面、配置、代码、数据、图件或运行日志。

## 1. 题面与统一定义

| 内容 | 证据 | 用途 |
|---|---|---|
| 原始训练题集 | `../../../pdf/training_examples_incomp.pdf` | Navier-Stokes 方程、PISO 法要求和腔流算例 |
| 统一问题定义 | `../../../0-caseDict/caseDict` | 第五题物理问题、边界、参考数据和验收标准 |

## 2. 算法实现

| 内容 | 证据 | 用途 |
|---|---|---|
| PISO 求解器源码 | `../../../UDF/solver/06_piso_navier_stokes_equation/pisoFoamStudent/pisoFoamStudent.C` | 动量方程、rAU、HbyA、压力方程和通量修正 |
| PISO 构建产物 | `../../../build/06_piso_navier_stokes_equation/bin/pisoFoamStudent` | 求解器构建结果 |
| 方腔后处理 | `../../../scripts/postprocess_lid_cavity.py` | Ghia 中心线采样、误差和图件 |
| 三角腔后处理 | `../../../scripts/postprocess_triangular_cavity.py` | 三角腔主涡、流函数和中心线 |
| 方腔混合网格生成器 | `../../../scripts/common/gmsh_hybrid_cavity_mesh.py` | 壁面层和方腔边界物理组 |
| 三角腔混合网格生成器 | `../../../scripts/common/gmsh_triangular_cavity_hybrid_mesh.py` | 三角腔几何、壁面层和混合拓扑 |

## 3. 方腔配置与结果

| 内容 | 证据 |
|---|---|
| Re1000 hybrid40 配置 | `../../../scripts/configs/05_navier_stokes_equation/16_lid_driven_cavity_piso_Re1000_hybrid40.json` |
| Re1000 hybrid80 配置 | `../../../scripts/configs/05_navier_stokes_equation/17_lid_driven_cavity_piso_Re1000_hybrid80.json` |
| Re3200 hybrid40 配置 | `../../../scripts/configs/05_navier_stokes_equation/18_lid_driven_cavity_piso_Re3200_hybrid40.json` |
| Re3200 hybrid80 配置 | `../../../scripts/configs/05_navier_stokes_equation/19_lid_driven_cavity_piso_Re3200_hybrid80.json` |
| Re1000 hybrid40 摘要 | `../../../data/05_navier_stokes_equation/cases/16_lid_driven_cavity_piso_Re1000_hybrid40/summary.json` |
| Re1000 hybrid80 摘要 | `../../../data/05_navier_stokes_equation/cases/17_lid_driven_cavity_piso_Re1000_hybrid80/summary.json` |
| Re3200 hybrid40 摘要 | `../../../data/05_navier_stokes_equation/cases/18_lid_driven_cavity_piso_Re3200_hybrid40/summary.json` |
| Re3200 hybrid80 摘要 | `../../../data/05_navier_stokes_equation/cases/19_lid_driven_cavity_piso_Re3200_hybrid80/summary.json` |
| Re1000 hybrid80 场图 | `../../../figures/05_navier_stokes_equation/cases/17_lid_driven_cavity_piso_Re1000_hybrid80/field_and_streamlines.png` |
| Re1000 hybrid80 中心线图 | `../../../figures/05_navier_stokes_equation/cases/17_lid_driven_cavity_piso_Re1000_hybrid80/centerline_comparison.png` |
| Re3200 hybrid80 场图 | `../../../figures/05_navier_stokes_equation/cases/19_lid_driven_cavity_piso_Re3200_hybrid80/field_and_streamlines.png` |
| Re3200 hybrid80 中心线图 | `../../../figures/05_navier_stokes_equation/cases/19_lid_driven_cavity_piso_Re3200_hybrid80/centerline_comparison.png` |
| Re1000 hybrid40 solver log | `../../../cases/05_navier_stokes_equation/16_lid_driven_cavity_piso_Re1000_hybrid40/log.pisoFoamStudent` |
| Re1000 hybrid80 solver log | `../../../cases/05_navier_stokes_equation/17_lid_driven_cavity_piso_Re1000_hybrid80/log.pisoFoamStudent` |
| Re3200 hybrid40 solver log | `../../../cases/05_navier_stokes_equation/18_lid_driven_cavity_piso_Re3200_hybrid40/N40/log.pisoFoamStudent` |
| Re3200 hybrid80 solver log | `../../../cases/05_navier_stokes_equation/19_lid_driven_cavity_piso_Re3200_hybrid80/N80/log.pisoFoamStudent` |

每个方腔结果目录还包括：

```text
u_centerline.csv
v_centerline.csv
summary.json
```

## 4. 三角腔配置与结果

| 内容 | 证据 |
|---|---|
| Re100 配置 | `../../../scripts/configs/05_navier_stokes_equation/29_triangular_cavity_piso_Re100_hybrid80.json` |
| Re200 配置 | `../../../scripts/configs/05_navier_stokes_equation/30_triangular_cavity_piso_Re200_hybrid80.json` |
| Re500 配置 | `../../../scripts/configs/05_navier_stokes_equation/31_triangular_cavity_piso_Re500_hybrid80.json` |
| Re100 摘要 | `../../../data/05_navier_stokes_equation/cases/29_triangular_cavity_piso_Re100_hybrid80/summary.json` |
| Re200 摘要 | `../../../data/05_navier_stokes_equation/cases/30_triangular_cavity_piso_Re200_hybrid80/summary.json` |
| Re500 摘要 | `../../../data/05_navier_stokes_equation/cases/31_triangular_cavity_piso_Re500_hybrid80/summary.json` |
| Re100 场图 | `../../../figures/05_navier_stokes_equation/cases/29_triangular_cavity_piso_Re100_hybrid80/field_streamlines.png` |
| Re200 场图 | `../../../figures/05_navier_stokes_equation/cases/30_triangular_cavity_piso_Re200_hybrid80/field_streamlines.png` |
| Re500 场图 | `../../../figures/05_navier_stokes_equation/cases/31_triangular_cavity_piso_Re500_hybrid80/field_streamlines.png` |
| Re100 流函数图 | `../../../figures/05_navier_stokes_equation/cases/29_triangular_cavity_piso_Re100_hybrid80/streamfunction_vortices.png` |
| Re200 流函数图 | `../../../figures/05_navier_stokes_equation/cases/30_triangular_cavity_piso_Re200_hybrid80/streamfunction_vortices.png` |
| Re500 流函数图 | `../../../figures/05_navier_stokes_equation/cases/31_triangular_cavity_piso_Re500_hybrid80/streamfunction_vortices.png` |
| Re100 solver log | `../../../cases/05_navier_stokes_equation/29_triangular_cavity_piso_Re100_hybrid80/N80/log.pisoFoamStudent` |
| Re200 solver log | `../../../cases/05_navier_stokes_equation/30_triangular_cavity_piso_Re200_hybrid80/N80/log.pisoFoamStudent` |
| Re500 solver log | `../../../cases/05_navier_stokes_equation/31_triangular_cavity_piso_Re500_hybrid80/N80/log.pisoFoamStudent` |

三角腔结果目录还包括：

```text
u_centerline.csv
v_horizontal.csv
summary.json
```

## 5. 网格检查与稳态证据

| 内容 | 证据 |
|---|---|
| 方腔 hybrid40 网格检查 | `../../../cases/05_navier_stokes_equation/16_lid_driven_cavity_piso_Re1000_hybrid40/log.checkMesh` |
| 方腔 hybrid80 网格检查 | `../../../cases/05_navier_stokes_equation/17_lid_driven_cavity_piso_Re1000_hybrid80/N80/log.checkMesh` |
| 三角腔 Re100 网格检查 | `../../../cases/05_navier_stokes_equation/29_triangular_cavity_piso_Re100_hybrid80/N80/log.checkMesh` |
| 三角腔 Re200 网格检查 | `../../../cases/05_navier_stokes_equation/30_triangular_cavity_piso_Re200_hybrid80/N80/log.checkMesh` |
| 三角腔 Re500 网格检查 | `../../../cases/05_navier_stokes_equation/31_triangular_cavity_piso_Re500_hybrid80/N80/log.checkMesh` |

所有上述 PISO solver log 的结尾均包含 `Steady state reached` 和 `End`，
并记录了最终的 `max |div(phi)|` 与 `max |U-Uprevious|`。

## 6. 证据使用规则

数值结论优先追溯到 `summary.json`；
中心线结果追溯到 `u_centerline.csv`、`v_centerline.csv`；
主涡和流函数结果追溯到三角腔摘要及 `streamfunction_vortices.png`；
稳态结论必须追溯到对应的 PISO solver log；
网格质量结论追溯到 `log.checkMesh`。
