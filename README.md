# 学生版线性对流 CFD 工程

本项目用于学习并实现 PDF 第 1 题：

$$\frac{\partial \phi}{\partial t}+\nabla\cdot(\boldsymbol{u}\phi)=0.$$

当前已经实现的是二维周期正弦波算例：

- 结构化四边形网格；
- 显式时间推进；
- `Gauss upwind` 和 `Gauss linearUpwind grad(T)` 两种空间格式；
- `N10/N20/N40/N80` 多网格误差与收敛阶分析。

三角形网格和固体旋转算例已经建立配置文件，但目前仍标记为规划项，尚未接入可运行的网格生成和后处理逻辑。

## 目录结构

```text
student_project/
├── README.md
├── UDF/
│   └── solver/
│       └── explicitAdvectionFoamStudent/
│
├── build/
│   └── bin/
│       └── explicitAdvectionFoamStudent
│
├── cases/
│   ├── 01_sine_wave_quad/
│   │   ├── N10/
│   │   ├── N20/
│   │   ├── N40/
│   │   └── N80/
│   └── 02_sine_wave_quad_linearUpwind/
│       ├── N10/
│       ├── N20/
│       ├── N40/
│       └── N80/
│
├── scripts/
│   ├── build_student_solver.sh
│   ├── prepare_case.py
│   ├── run_case.py
│   ├── postprocess_case.py
│   ├── run_study.py
│   ├── collect_results.py
│   ├── analyze_study.py
│   ├── plot_study.py
│   ├── configs/
│   │   ├── 01_sine_wave_quad_upwind.json
│   │   ├── 02_sine_wave_quad_linearUpwind.json
│   │   ├── 03_sine_wave_tri_upwind.json
│   │   └── 04_solid_rotation_quad_upwind.json
│   └── common/
│       ├── paths.py
│       ├── case_config.py
│       ├── foam_case.py
│       ├── mesh_tools.py
│       ├── advection_sine.py
│       ├── advection_rotation.py
│       ├── metrics.py
│       └── plotting.py
│
├── data/
│   ├── cases/
│   └── analysis/
│
├── figures/
│   ├── cases/
│   └── analysis/
│
└── docs/
    ├── 01/
    └── compare/
```

旧的 `scripts/01_sine_wave_quad/` 和 `scripts/02_sine_wave_quad_linearUpwind/` 已删除。现在所有案例统一通过顶层脚本和 JSON 配置运行。

## 配置文件

同一个题目内部，不同案例通过 JSON 配置文件区分。配置文件记录：

| 参数 | 含义 |
|---|---|
| `caseName` | 对应 `cases/` 下的案例目录 |
| `problem` | 物理问题类型 |
| `meshType` | 网格类型 |
| `schemeName` | 人类可读的格式名称 |
| `divScheme` | 写入 `system/fvSchemes` 的对流离散格式 |
| `solver` | 要调用的 OpenFOAM 求解器 |
| `resolutions` | 默认的 N 序列 |
| `endTime` | 终止时间 |
| `maxCo` | 目标 Courant 数 |

当前可运行配置：

```text
scripts/configs/01_sine_wave_quad_upwind.json
scripts/configs/02_sine_wave_quad_linearUpwind.json
```

规划配置：

```text
scripts/configs/03_sine_wave_tri_upwind.json
scripts/configs/04_solid_rotation_quad_upwind.json
```

规划配置中的 `"implemented"` 为 `false`，脚本会拒绝运行，避免误以为已经完成。

## 编译求解器

```bash
cd /home/a776/workdocuments/上交船舶/slover/student_project
sh scripts/build_student_solver.sh
```

编译产物：

```text
build/bin/explicitAdvectionFoamStudent
```

## 只准备案例

准备 01 案例的 `N40`，不运行 OpenFOAM：

```bash
python3 scripts/prepare_case.py \
    --config scripts/configs/01_sine_wave_quad_upwind.json \
    --N 40 \
    --overwrite
```

准备全部默认 N：

```bash
python3 scripts/run_study.py \
    --config scripts/configs/01_sine_wave_quad_upwind.json \
    --prepare-only \
    --overwrite
```

02 案例只替换配置文件：

```bash
python3 scripts/run_study.py \
    --config scripts/configs/02_sine_wave_quad_linearUpwind.json \
    --prepare-only \
    --overwrite
```

## 运行一个 N

运行 OpenFOAM 前：

```bash
source /opt/openfoam14/etc/bashrc
```

运行并后处理 01 案例的 `N40`：

```bash
python3 scripts/run_case.py \
    --config scripts/configs/01_sine_wave_quad_upwind.json \
    --N 40 \
    --overwrite
```

如果案例已经准备好，只想运行：

```bash
python3 scripts/run_case.py \
    --config scripts/configs/01_sine_wave_quad_upwind.json \
    --N 40 \
    --no-prepare
```

如果已有 OpenFOAM 结果，只做后处理：

```bash
python3 scripts/postprocess_case.py \
    --config scripts/configs/01_sine_wave_quad_upwind.json \
    --N 40
```

## 运行全部 N

```bash
source /opt/openfoam14/etc/bashrc

python3 scripts/run_study.py \
    --config scripts/configs/01_sine_wave_quad_upwind.json \
    --resolutions 10,20,40,80 \
    --overwrite
```

这个命令会对每个 N 执行：

```text
prepare_case
    -> Allclean
    -> 生成 0.orig/T
    -> blockMesh
    -> checkMesh
    -> 求解器
    -> 单网格后处理
```

全部 N 完成后自动执行：

```text
collect_results
    -> analyze_study
    -> plot_study
```

## 只重新分析已有结果

如果每个 N 已经有：

```text
data/cases/<caseName>/Nxx/summary.json
```

可以只收集：

```bash
python3 scripts/collect_results.py \
    --config scripts/configs/01_sine_wave_quad_upwind.json \
    --resolutions 10,20,40,80
```

计算收敛表和分析报告：

```bash
python3 scripts/analyze_study.py \
    --config scripts/configs/01_sine_wave_quad_upwind.json
```

重新绘制总体图：

```bash
python3 scripts/plot_study.py \
    --config scripts/configs/01_sine_wave_quad_upwind.json
```

这些命令不会运行 OpenFOAM，只读取已有结果。

## 输出位置

单个 N 的数据：

```text
data/cases/<caseName>/N40/
├── summary.json
├── time_history.csv
├── field_data.csv
└── error_field.csv
```

单个 N 的图片：

```text
figures/cases/<caseName>/N40/
├── field_comparison.png
├── diagonal_profile.png
├── amplitude_history.png
└── cfl_history.png
```

所有 N 的分析：

```text
data/analysis/<caseName>/
├── raw_results.csv
├── convergence_summary.csv
├── run_manifest.json
└── analysis.md
```

所有 N 的总体图片：

```text
figures/analysis/<caseName>/
├── convergence_errors.png
├── convergence_order.png
└── all_N_comparison.png
```

## 收敛阶

最终时刻的归一化误差为：

$$L_1=\frac{\sum_c V_c|T_c-T_c^{exact}|}{\sum_c V_c|T_c^{exact}|}.$$

当网格由 $N$ 加密到 $2N$ 时：

$$p=\frac{\log(E_N/E_{2N})}{\log(2)}.$$

实际代码位置：

```text
scripts/common/advection_tools.py
    normalized_errors()

scripts/common/study_analysis.py
    observed_order()
    analyse()
```

## 当前实现边界

已实现：

- 线性对流方程；
- 二维周期正弦波；
- 结构化四边形网格；
- `Gauss upwind`；
- `Gauss linearUpwind grad(T)`；
- 单个 N 后处理；
- 多个 N 的误差和收敛阶分析。

暂未实现：

- 三角形网格生成和读取；
- 固体旋转复杂剖面；
- slotted disk、cone、smooth hump 初始场；
- 非结构网格上的专用可视化；
- 第 1 题之外的扩散、Poisson 和 Navier-Stokes 求解器。

清理单个案例：

```bash
sh cases/01_sine_wave_quad/N40/Allclean
```

`Allclean` 只删除可重新生成的运行产物，不删除 `0.orig/`、`system/` 和常量输入文件。
