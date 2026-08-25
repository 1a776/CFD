# `cases/` 案例目录说明

`cases/` 存放本项目每一个题目的 OpenFOAM 运行具体案例。
## 1. 目录树

```text
cases/
├── README.md
│
├── 01_sine_wave_quad/
│   ├── N10/
│   ├── N20/
│   ├── N40/
│   └── N80/
│
├── 02_sine_wave_quad_linearUpwind/
│   ├── N10/
│   ├── N20/
│   ├── N40/
│   └── N80/
│
├── 03_sine_wave_tri_upwind/
│   ├── N10/
│   ├── N20/
│   ├── N40/
│   └── N80/
│
├── 03_sine_wave_tri_linearUpwind/
│   ├── N10/
│   ├── N20/
│   ├── N40/
│   └── N80/
│
├── 04_solid_rotation_quad_upwind/
│   ├── N50/
│   ├── N100/
│   └── N200/
│
└── 04_solid_rotation_tri_upwind/
    ├── N50/
    ├── N100/
    └── N200/
```



## 2. 案例族、配置和结果索引

`cases/` 是 OpenFOAM 原始案例目录；真正用于报告分析的数据和图片放在 `data/`、`figures/` 和 `report/` 中。

- 第一题总报告和证据索引
  - 第一题报告系统分析二维线性对流方程中正弦波平移和复杂轮廓固体旋转两个算例在不同网格与格式下的数值表现。
  - 汇总报告：[`report/01/report.md`](../report/01/report.md)，分析各案例的题目定义、网格设置、格式选择、运行结果和主要结论。
  - 证据索引：[`report/01/evidence_index.md`](../report/01/evidence_index.md)，把报告结论对应到 JSON 配置、OpenFOAM 运行目录、日志文件、数据表和图片。
  - 覆盖的案例目录：
    - `01_sine_wave_quad/`
    - `02_sine_wave_quad_linearUpwind/`
    - `03_sine_wave_tri_upwind/`
    - `03_sine_wave_tri_linearUpwind/`
    - `04_solid_rotation_quad_upwind/`
    - `04_solid_rotation_tri_upwind/`

- `01_sine_wave_quad/`
  - 题目/算例：第一题，正弦波平移。
  - 网格类型：结构化四边形网格。
  - 空间插值格式：`Gauss upwind`，一阶迎风。
  - `N`：`10, 20, 40, 80`。
  - 速度或速度模型：常速度 `U=(1,1,0)`。
  - 终止时间：`1.0`。
  - 主要检查内容：`L1` 误差、网格收敛阶、数值耗散。
  - JSON 配置：[`01_sine_wave_quad_upwind.json`](../scripts/configs/01_sine_wave_quad_upwind.json)。
  - 结果数据：[`data/cases/01_sine_wave_quad/`](../data/cases/01_sine_wave_quad/)，[`data/analysis/01_sine_wave_quad/`](../data/analysis/01_sine_wave_quad/)。
  - 结果图片：[`figures/cases/01_sine_wave_quad/`](../figures/cases/01_sine_wave_quad/)，[`figures/analysis/01_sine_wave_quad/`](../figures/analysis/01_sine_wave_quad/)。
  - 分析说明：[`analysis.md`](../data/analysis/01_sine_wave_quad/analysis.md) 主要记录正弦波一阶迎风四边形网格的 `L1/L2/Linf` 误差、质量守恒、CFL 历史、振幅衰减和观察收敛阶。

- `02_sine_wave_quad_linearUpwind/`
  - 题目/算例：第一题，正弦波平移。
  - 网格类型：结构化四边形网格。
  - 空间插值格式：`Gauss linearUpwind grad(T)`，线性迎风。
  - `N`：`10, 20, 40, 80`。
  - 速度或速度模型：常速度 `U=(1,1,0)`。
  - 终止时间：`1.0`。
  - 主要检查内容：与一阶迎风比较误差、耗散和过冲。
  - JSON 配置：[`02_sine_wave_quad_linearUpwind.json`](../scripts/configs/02_sine_wave_quad_linearUpwind.json)。
  - 结果数据：[`data/cases/02_sine_wave_quad_linearUpwind/`](../data/cases/02_sine_wave_quad_linearUpwind/)，[`data/analysis/02_sine_wave_quad_linearUpwind/`](../data/analysis/02_sine_wave_quad_linearUpwind/)。
  - 结果图片：[`figures/cases/02_sine_wave_quad_linearUpwind/`](../figures/cases/02_sine_wave_quad_linearUpwind/)，[`figures/analysis/02_sine_wave_quad_linearUpwind/`](../figures/analysis/02_sine_wave_quad_linearUpwind/)。
  - 分析说明：[`analysis.md`](../data/analysis/02_sine_wave_quad_linearUpwind/analysis.md) 主要记录正弦波线性迎风四边形网格的误差、收敛阶、耗散变化和与一阶迎风相比的精度表现。

- `03_sine_wave_tri_upwind/`
  - 题目/算例：第一题，正弦波平移。
  - 网格类型：Gmsh 三角形棱柱网格。
  - 空间插值格式：`Gauss upwind`，一阶迎风。
  - `N`：`10, 20, 40, 80`。
  - 速度或速度模型：常速度 `U=(1,1,0)`。
  - 终止时间：`1.0`。
  - 主要检查内容：三角形网格上的误差和收敛阶。
  - JSON 配置：[`03_sine_wave_tri_upwind.json`](../scripts/configs/03_sine_wave_tri_upwind.json)。
  - 结果数据：[`data/cases/03_sine_wave_tri_upwind/`](../data/cases/03_sine_wave_tri_upwind/)，[`data/analysis/03_sine_wave_tri_upwind/`](../data/analysis/03_sine_wave_tri_upwind/)。
  - 结果图片：[`figures/cases/03_sine_wave_tri_upwind/`](../figures/cases/03_sine_wave_tri_upwind/)，[`figures/analysis/03_sine_wave_tri_upwind/`](../figures/analysis/03_sine_wave_tri_upwind/)。
  - 分析说明：[`analysis.md`](../data/analysis/03_sine_wave_tri_upwind/analysis.md) 主要记录正弦波一阶迎风三角形棱柱网格的误差、收敛阶、真实 cell centre/cell volume 积分和网格类型影响。

- `03_sine_wave_tri_linearUpwind/`
  - 题目/算例：第一题，正弦波平移。
  - 网格类型：Gmsh 三角形棱柱网格。
  - 空间插值格式：`Gauss linearUpwind grad(T)`。
  - `N`：`10, 20, 40, 80`。
  - 速度或速度模型：常速度 `U=(1,1,0)`。
  - 终止时间：`1.0`。
  - 主要检查内容：三角形网格上线性迎风的精度和稳定性。
  - JSON 配置：[`03_sine_wave_tri_linearUpwind.json`](../scripts/configs/03_sine_wave_tri_linearUpwind.json)。
  - 结果数据：[`data/cases/03_sine_wave_tri_linearUpwind/`](../data/cases/03_sine_wave_tri_linearUpwind/)，[`data/analysis/03_sine_wave_tri_linearUpwind/`](../data/analysis/03_sine_wave_tri_linearUpwind/)。
  - 结果图片：[`figures/cases/03_sine_wave_tri_linearUpwind/`](../figures/cases/03_sine_wave_tri_linearUpwind/)，[`figures/analysis/03_sine_wave_tri_linearUpwind/`](../figures/analysis/03_sine_wave_tri_linearUpwind/)。
  - 分析说明：[`analysis.md`](../data/analysis/03_sine_wave_tri_linearUpwind/analysis.md) 主要记录正弦波线性迎风三角形棱柱网格的误差、收敛阶、真实 cell centre/cell volume 积分和稳定性表现。

- `04_solid_rotation_quad_upwind/`
  - 题目/算例：第一题，复杂轮廓固体旋转。
  - 网格类型：结构化四边形网格。
  - 空间插值格式：`Gauss upwind`，一阶迎风。
  - `N`：`50, 100, 200`。
  - 速度或速度模型：绕 `(0.5,0.5)` 旋转，角速度 `1`。
  - 终止时间：`2π`。
  - 主要检查内容：最终轮廓、振幅衰减、数值耗散。
  - JSON 配置：[`04_solid_rotation_quad_upwind.json`](../scripts/configs/04_solid_rotation_quad_upwind.json)。
  - 结果数据：[`data/cases/04_solid_rotation_quad_upwind/`](../data/cases/04_solid_rotation_quad_upwind/)。
  - 结果图片：[`figures/cases/04_solid_rotation_quad_upwind/`](../figures/cases/04_solid_rotation_quad_upwind/)。

- `04_solid_rotation_tri_upwind/`
  - 题目/算例：第一题，复杂轮廓固体旋转。
  - 网格类型：Gmsh 三角形棱柱网格。
  - 空间插值格式：`Gauss upwind`，一阶迎风。
  - `N`：`50, 100, 200`。
  - 速度或速度模型：绕 `(0.5,0.5)` 旋转，角速度 `1`。
  - 终止时间：`2π`。
  - 主要检查内容：三角形和四边形网格的旋转结果比较。
  - JSON 配置：[`04_solid_rotation_tri_upwind.json`](../scripts/configs/04_solid_rotation_tri_upwind.json)。
  - 结果数据：[`data/cases/04_solid_rotation_tri_upwind/`](../data/cases/04_solid_rotation_tri_upwind/)。
  - 结果图片：[`figures/cases/04_solid_rotation_tri_upwind/`](../figures/cases/04_solid_rotation_tri_upwind/)。
