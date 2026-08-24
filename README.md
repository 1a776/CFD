# 学生版线性对流 CFD 工程

本项目对应教学材料中的二维周期正弦波线性对流算例。当前求解器使用有限体积法、显式前向 Euler 时间推进，并通过 `CFL=0.2` 自动调整时间步。

当前案例：

- `01_sine_wave_quad`：四边形网格、一阶迎风插值；
- `02_sine_wave_quad_linearUpwind`：四边形网格、`linearUpwind` 插值，用于格式对比。

## 目录结构

```text
student_project/
├── README.md
├── UDF/
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
│   ├── 02_sine_wave_quad_linearUpwind/
│   │   ├── N10/
│   │   ├── N20/
│   │   ├── N40/
│   │   └── N80/
│   └── 后续案例/
│
├── scripts/
│   ├── build_student_solver.sh
│   ├── 01_sine_wave_quad/
│   │   ├── create_initial_fields.py
│   │   ├── measure_error.py
│   │   ├── plot_results.py
│   │   ├── run_all_N.py
│   │   ├── collect_results.py
│   │   ├── plot_convergence.py
│   │   ├── analyze_results.py
│   │   └── README.md
│   ├── 02_sine_wave_quad_linearUpwind/
│   │   ├── create_initial_fields.py
│   │   ├── measure_error.py
│   │   ├── plot_results.py
│   │   ├── run_all_N.py
│   │   ├── collect_results.py
│   │   ├── analyze_results.py
│   │   ├── plot_convergence.py
│   │   └── README.md
│   └── common/
│       └── 通用辅助函数
│
├── data/
│   ├── cases/
│   │   ├── 01_sine_wave_quad/N10/...
│   │   ├── 01_sine_wave_quad/N20/...
│   │   └── 02_sine_wave_quad_linearUpwind/N20/...
│   └── analysis/
│       ├── 01_sine_wave_quad/
│       └── 02_sine_wave_quad_linearUpwind/
│
├── figures/
│   ├── cases/
│   │   ├── 01_sine_wave_quad/N10/...
│   │   └── 02_sine_wave_quad_linearUpwind/N20/...
│   └── analysis/
│       ├── 01_sine_wave_quad/
│       └── 02_sine_wave_quad_linearUpwind/
│
└── docs/
    ├── 01/
    └── compare/
```

`cases/` 只保存 OpenFOAM 算例输入和运行产物，不再放 Python 分析脚本。每个 `Nxx` 都是一个可以独立运行的完整算例，包含 `0.orig/`、`system/`、`constant/`、`Allrun` 和 `Allclean`。

## 编译求解器

```bash
cd /home/a776/workdocuments/上交船舶/slover/student_project
sh scripts/build_student_solver.sh
```

编译产物：

```text
build/bin/explicitAdvectionFoamStudent
```

## 运行 01 案例的全部网格

```bash
cd /home/a776/workdocuments/上交船舶/slover/student_project
source /opt/openfoam14/etc/bashrc
python3 scripts/01_sine_wave_quad/run_all_N.py --overwrite
```

默认运行：

```text
N10, N20, N40, N80
```

脚本会自动生成对应网格、生成初始场、运行 `blockMesh`、`checkMesh`、学生求解器，并进行单网格后处理和所有网格的收敛分析。

## 运行 02 案例的全部网格

```bash
cd /home/a776/workdocuments/上交船舶/slover/student_project
source /opt/openfoam14/etc/bashrc
python3 scripts/02_sine_wave_quad_linearUpwind/run_all_N.py --overwrite
```

## 单个 N 的运行和后处理

例如 `01` 案例的 `N20`：

```bash
source /opt/openfoam14/etc/bashrc
sh cases/01_sine_wave_quad/N20/Allrun

python3 scripts/01_sine_wave_quad/plot_results.py \
    --case-dir cases/01_sine_wave_quad/N20
```

快速查看误差：

```bash
python3 scripts/01_sine_wave_quad/measure_error.py \
    --case-dir cases/01_sine_wave_quad/N20
```

## 输出位置

每个 N 的数据：

```text
data/cases/01_sine_wave_quad/N20/
├── summary.json
├── time_history.csv
├── field_data.csv
└── error_field.csv
```

每个 N 的图：

```text
figures/cases/01_sine_wave_quad/N20/
├── field_comparison.png
├── diagonal_profile.png
├── amplitude_history.png
└── cfl_history.png
```

全部 N 的分析数据：

```text
data/analysis/01_sine_wave_quad/
├── raw_results.csv
├── convergence_summary.csv
├── run_manifest.json
└── analysis.md
```

全部 N 的分析图：

```text
figures/analysis/01_sine_wave_quad/
├── convergence_errors.png
├── convergence_order.png
└── all_N_comparison.png
```

## 只重新分析已有结果

如果四个 N 已经运行完成，可以跳过求解器：

```bash
python3 scripts/01_sine_wave_quad/collect_results.py
python3 scripts/01_sine_wave_quad/analyze_results.py
```

重新绘制收敛图：

```bash
python3 scripts/01_sine_wave_quad/plot_convergence.py
```

## 清理规则

每个分辨率算例自己的运行产物由对应的 `Allclean` 管理：

```bash
sh cases/01_sine_wave_quad/N20/Allclean
```

它会清理时间目录、`0/`、`constant/polyMesh/`、`postProcessing/` 和日志，但保留 `0.orig/`、`system/` 和手工维护的常量文件。

## 学习文档

```text
docs/01/00_learning_path.md
docs/01/09_stage6_visualization_and_convergence.md
docs/01/10_plotting_commands.md
docs/01/11_periodic_boundary_explanation.md
docs/01/12_project_structure.md
docs/compare/0102.md
docs/compare/命令.md
```
