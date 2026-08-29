# 不可压 CFD 有限体积案例项目

本项目围绕 [pdf/training_examples_incomp.pdf](pdf/training_examples_incomp.pdf) 中的五个不可压 CFD 数值测试题目展开，目标是完成对应方程的有限体积离散、OpenFOAM 学生版求解器实现、JSON 驱动案例生成、计算运行、后处理、数据分析、画图和报告归档。


## 题目索引

- ★ **每个案例的 OpenFOAM 目录、数据目录、图片目录、报告和证据索引的完整链接见 [cases/README.md](cases/README.md)。**
- ★ **题目解答整理见 [pdf/题目解答.pdf](pdf/题目解答.pdf)。**

| 题目 | 求解器 | 案例配置 | 报告 | 证据索引 |
|---|---|---|---|---|
| [1. 对流方程](cases/README.md#advection-equation) | [explicitAdvectionFoamStudent](UDF/README.md#udf-advection-equation) | [正弦波平移](cases/README.md#advection-sine-wave-translation-config)、[复杂轮廓固体旋转](cases/README.md#advection-solid-rotation-config) | [report/01_advection_equation/report.md](report/01_advection_equation/report.md) | [report/01_advection_equation/evidence_index.md](report/01_advection_equation/evidence_index.md) |
| [2. 扩散方程](cases/README.md#diffusion-equation) | [explicitDiffusionFoamStudent](UDF/README.md#udf-diffusion-equation) | [间断初值扩散](cases/README.md#diffusion-discontinuous-initial-config)、[Gaussian 扩散](cases/README.md#diffusion-gaussian-config) | [report/02_diffusion_equation/report.md](report/02_diffusion_equation/report.md) | [report/02_diffusion_equation/evidence_index.md](report/02_diffusion_equation/evidence_index.md) |
| [3. 对流-扩散方程](cases/README.md#advection-diffusion-equation) | [explicitAdvectionDiffusionFoamStudent](UDF/README.md#udf-advection-diffusion-equation) | [正弦波平移](cases/README.md#advection-diffusion-sine-wave-translation-config)、[旋转尖峰平流扩散](cases/README.md#advection-diffusion-rotating-peak-config) | [report/03_advection_diffusion_equation/report.md](report/03_advection_diffusion_equation/report.md) | [report/03_advection_diffusion_equation/evidence_index.md](report/03_advection_diffusion_equation/evidence_index.md) |
| [4. Poisson 方程](cases/README.md#poisson-equation) | [poissonFoamStudent](UDF/README.md#udf-poisson-equation) | [制造解 Poisson 方程](cases/README.md#poisson-manufactured-solution-config) | [report/04_poisson_equation/report.md](report/04_poisson_equation/report.md) | [report/04_poisson_equation/evidence_index.md](report/04_poisson_equation/evidence_index.md) |
| [5. Navier-Stokes 方程](cases/README.md#navier-stokes-equation) | [projectionFoamStudent / pisoFoamStudent](UDF/README.md#udf-navier-stokes-equation) | [方腔顶盖驱动流](cases/README.md#navier-stokes-lid-driven-cavity-config)、[等边三角腔顶盖驱动流](cases/README.md#navier-stokes-triangular-cavity-config) | [投影法报告](report/05_navier_stokes_equation/projection/report.md)、[PISO 报告](report/05_navier_stokes_equation/piso/report.md) | [投影法证据](report/05_navier_stokes_equation/projection/evidence_index.md)、[PISO 证据](report/05_navier_stokes_equation/piso/evidence_index.md) |


## 项目架构

| 目录 | 内容 | 入口说明 |
|---|---|---|
| [UDF/](UDF/) | 学生版 OpenFOAM 求解器源码 | 求解器和题目/案例的对应关系见 [UDF/README.md](UDF/README.md) |
| [scripts/](scripts/) | JSON 解析、case 构建、运行、后处理、数据汇总和画图脚本 | 配置文件在 [scripts/configs/](scripts/configs/) |
| [cases/](cases/) | 由 JSON 生成或保留的 OpenFOAM 案例目录 | 案例配置、数据和图片索引见 [cases/README.md](cases/README.md) |
| [data/](data/) | 每个案例的后处理数据、误差表、剖面数据和汇总结果 | 由 `scripts/postprocess_*.py`、`scripts/collect_results.py`、`scripts/analyze_study.py` 生成 |
| [figures/](figures/) | 每个案例的场图、剖面图、误差曲线和收敛阶图 | 由后处理和绘图脚本生成 |
| [report/](report/) | 实验报告和证据索引 | 每题报告与证据索引见下方题目索引 |
| [pdf/](pdf/) | 原始题目和题目解答 | [原题](pdf/training_examples_incomp.pdf)、[解答](pdf/题目解答.pdf) |
| [docs/](docs/) | 推导、调试和阶段性记录 | 第一题推导记录见 [docs/01_advection_equation/](docs/01_advection_equation/) |
| [build/](build/) | 本项目内的求解器编译输出和 Matplotlib 缓存 | 由构建脚本和绘图脚本自动生成 |

脚本中的项目根目录由脚本文件位置推导，不需要在脚本里写本机绝对路径。OpenFOAM 和 Gmsh Python 环境通过 `OPENFOAM_BASHRC`、`VIBEFLOW_PYTHON` 等环境变量指定。


## 复现环境检查

从项目根目录开始执行。默认假设 OpenFOAM 14 安装在 `/opt/openfoam14`；如果安装位置不同，先设置 `OPENFOAM_BASHRC`。

```bash
export OPENFOAM_BASHRC="${OPENFOAM_BASHRC:-/opt/openfoam14/etc/bashrc}"
test -f "$OPENFOAM_BASHRC"
. "$OPENFOAM_BASHRC"

command -v python3
command -v wmake
command -v blockMesh
command -v checkMesh
command -v foamPostProcess
command -v gmshToFoam
```

Python 后处理和 Gmsh 网格脚本需要 `numpy`、`scipy`、`matplotlib` 和 `gmsh`：

```bash
python3 -c "import numpy, scipy, matplotlib, gmsh; print('python dependencies ok')"
```

如果 `gmsh` 安装在单独的 Python 环境中，设置：

```bash
export VIBEFLOW_PYTHON=/path/to/python-with-gmsh
```

不设置 `VIBEFLOW_PYTHON` 时，生成的 `Allrun.pre` 会使用当前环境中的 `python3`。

## 编译求解器

```bash
. "$OPENFOAM_BASHRC"
sh scripts/build_student_solver.sh
```

编译产物写入本项目的 `build/` 目录，不写入 OpenFOAM 安装目录：

```text
build/01_advection_equation/bin/explicitAdvectionFoamStudent
build/02_diffusion_equation/bin/explicitDiffusionFoamStudent
build/03_advection_diffusion_equation/bin/explicitAdvectionDiffusionFoamStudent
build/04_poisson_equation/bin/poissonFoamStudent
build/05_navier_stokes_equation/bin/projectionFoamStudent
build/06_piso_navier_stokes_equation/bin/pisoFoamStudent
```

## 从 JSON 复现案例

下面用第四题制造解 Poisson 方程的四边形网格案例作为例子。这个例子从一个 JSON 配置开始，自动生成 `cases/` 中的 OpenFOAM 案例，运行求解器，然后生成 `data/` 和 `figures/` 中的分析结果。

```bash
python3 scripts/run_study.py \
  --config scripts/configs/04_poisson_equation/01_poisson_manufactured_quad.json \
  --overwrite
```

这条命令的含义如下：

- `--config scripts/configs/04_poisson_equation/01_poisson_manufactured_quad.json`：指定实验配置。JSON 中记录案例名、求解器、网格类型、离散格式、边界条件、分辨率列表和后处理设置。
- `--overwrite`：如果目标 case 已经存在，先删除旧 case，再按当前 JSON 重新生成。复现实验时建议加上这个参数，避免旧文件影响结果。
- 脚本会按 JSON 中的 `resolutions` 依次生成多个网格分辨率的 case；本例对应 `N=10,20,40,80`。
- 每个分辨率都会生成 `cases/04_poisson_equation/01_poisson_manufactured_quad/N*/`，并在其中写入 `0/`、`constant/`、`system/`、`Allrun` 和 `metadata.json`。
- 运行结束后，单个分辨率的误差、场数据和图片写入 `data/04_poisson_equation/cases/01_poisson_manufactured_quad/` 与 `figures/04_poisson_equation/cases/01_poisson_manufactured_quad/`。
- 多分辨率汇总和收敛图写入 `data/04_poisson_equation/analysis/01_poisson_manufactured_quad/` 与 `figures/04_poisson_equation/analysis/01_poisson_manufactured_quad/`。

快速检查时可以只跑一个分辨率，例如只生成并运行 `N=10`：

```bash
python3 scripts/run_study.py \
  --config scripts/configs/04_poisson_equation/01_poisson_manufactured_quad.json \
  --resolutions 10 \
  --overwrite
```

如果只想检查 JSON 能否正确生成 OpenFOAM case，而不运行求解器，可以加 `--prepare-only`：

```bash
python3 scripts/run_study.py \
  --config scripts/configs/04_poisson_equation/01_poisson_manufactured_quad.json \
  --prepare-only \
  --overwrite
```

如果 case 已经运行完，只想重新收集数据、分析和画图，可以使用：

```bash
python3 scripts/collect_results.py --config scripts/configs/04_poisson_equation/01_poisson_manufactured_quad.json
python3 scripts/analyze_study.py --config scripts/configs/04_poisson_equation/01_poisson_manufactured_quad.json
python3 scripts/plot_study.py --config scripts/configs/04_poisson_equation/01_poisson_manufactured_quad.json
```

其它题目按同样方式复现：把 `--config` 后面的 JSON 换成 [scripts/configs/](scripts/configs/) 中对应案例的配置文件即可。第一题到第四题通常直接使用 JSON 中的 `resolutions`；第五题每个 JSON 通常只对应一个 Reynolds 数、一个算法和一个网格等级，必要时用 `--resolutions` 明确指定 `40` 或 `80`。

例如第五题方腔顶盖驱动流的一个投影法案例：

```bash
python3 scripts/run_study.py \
  --config scripts/configs/05_navier_stokes_equation/07_lid_driven_cavity_projection_Re1000_hybrid40.json \
  --resolutions 40 \
  --overwrite
```

第五题的后处理产物主要是各 case 的 `summary.json`、中心线/水平剖面 CSV、速度场图、流线图、流函数图、涡量图和主涡数据；报告汇总见 [report/05_navier_stokes_equation/](report/05_navier_stokes_equation/)。

## 常用脚本入口

| 命令 | 作用 |
|---|---|
| `python3 scripts/prepare_case.py --config <json> --N <N> --overwrite` | 从 JSON 生成单个分辨率的 OpenFOAM case |
| `python3 scripts/run_case.py --config <json> --N <N> --overwrite` | 生成、运行并后处理单个分辨率 |
| `python3 scripts/run_study.py --config <json> --overwrite` | 按 JSON 中的 `resolutions` 批量运行 |
| `python3 scripts/postprocess_case.py --config <json> --N <N>` | 对已运行的前四题/对流-扩散 case 重新后处理 |
| `python3 scripts/collect_results.py --config <json>` | 收集多个分辨率的 `summary.json` |
| `python3 scripts/analyze_study.py --config <json>` | 生成误差和观察收敛阶表 |
| `python3 scripts/plot_study.py --config <json>` | 生成误差曲线和收敛阶图 |

## 结果阅读顺序

1. 先看 [pdf/training_examples_incomp.pdf](pdf/training_examples_incomp.pdf) 和 [pdf/题目解答.pdf](pdf/题目解答.pdf)，确认题目和离散思路。
2. 再看 [UDF/README.md](UDF/README.md)，确认每个求解器对应哪个题目和案例。
3. 然后看 [cases/README.md](cases/README.md)，进入具体实验配置、OpenFOAM case、数据、图片和报告链接。
4. 最后看 [report/](report/) 中各题报告和证据索引，核对图表来源、运行记录和结论。
