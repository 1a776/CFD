# CFD 学生版显式有限体积求解器项目

## 1. 项目总览

本项目围绕题目资料：

[`pdf/training_examples_incomp.pdf`](pdf/training_examples_incomp.pdf)

开展有限体积法、OpenFOAM 求解器开发、网格研究、误差分析和数值验证。
该 PDF 是整组训练题的题目与参考材料集合，内容覆盖对流、扩散、对流扩散以及
后续不可压流体相关的数值问题。

当前仓库重点完成的是第一题二维线性对流方程：

$$
\frac{\partial T}{\partial t}
+
\nabla\cdot(\boldsymbol{U}T)=0.
$$

当前已经围绕第一题完成：

- OpenFOAM 14 学生版显式对流求解器；
- 正弦波平移算例；
- 复杂轮廓固体旋转算例；
- 结构化四边形网格；
- 非结构三角形网格；
- 一阶迎风插值；
- `linearUpwind` 梯度重构插值；
- 多个 $N$ 的网格研究；
- $L_1$、$L_2$、$L_\infty$ 误差计算；
- 网格观察收敛阶计算；
- CFL、质量守恒和最终场极值检查；
- 自动生成结果表格、分析文件和图片；
- 第一题完整验证报告。

第一题正式报告：

[`report/01/report.md`](report/01/report.md)

报告证据索引：

[`report/01/evidence_index.md`](report/01/evidence_index.md)

本项目目前不是把所有训练题都混在一个求解器里，而是按照题目分阶段开发。
当前第一题已经形成了“数学推导 -> OpenFOAM 求解器 -> 可运行案例 -> 多分辨率实验
-> 后处理 -> 验证报告”的完整链条。后续题目可以复用项目中的脚本组织方式，
但通常需要根据新的控制方程重新设计求解器核心。

## 2. 项目整体架构

项目根目录结构如下：

```text
student_project/
├── README.md
├── UDF/
│   ├── README.md
│   └── solver/
│       └── explicitAdvectionFoamStudent/
│           ├── explicitAdvectionFoamStudent.C
│           └── Make/
│               ├── files
│               └── options
│
├── build/
│   └── bin/
│       └── explicitAdvectionFoamStudent
│
├── cases/
│   ├── 01_sine_wave_quad/
│   ├── 02_sine_wave_quad_linearUpwind/
│   ├── 03_sine_wave_tri_upwind/
│   ├── 03_sine_wave_tri_linearUpwind/
│   ├── 04_solid_rotation_quad_upwind/
│   └── 04_solid_rotation_tri_upwind/
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
│   └── common/
│
├── data/
│   ├── cases/
│   └── analysis/
│
├── figures/
│   ├── cases/
│   └── analysis/
│
├── docs/
│   ├── 01/
│   ├── compare/
│   ├── bug_log.md
│   ├── changelog.md
│   └── changelog_intro.md
│
├── report/
│   └── 01/
│       ├── report.md
│       └── evidence_index.md
│
└── pdf/
    ├── training_examples_incomp.pdf
    ├── 01/
    ├── 02/
    ├── 03/
    ├── tex/
    └── 参考资料/
```

## 3. 各目录的职责

### 3.1 `UDF/`：学生版求解器源代码

`UDF/` 存放本项目自己开发的 OpenFOAM 代码，不修改
`/opt/openfoam14/src` 中的 OpenFOAM 源码。

```text
UDF/
├── README.md
└── solver/
    └── explicitAdvectionFoamStudent/
        ├── explicitAdvectionFoamStudent.C
        └── Make/
            ├── files
            └── options
```

主要作用：

- `explicitAdvectionFoamStudent.C`：读取网格、速度场和标量场，计算面通量，
  根据 CFL 设置时间步，计算显式对流残差并更新标量场；
- `Make/files`：告诉 `wmake` 编译哪个源文件以及生成哪个可执行文件；
- `Make/options`：指定 OpenFOAM 头文件路径和链接库；
- `UDF/README.md`：记录求解器的数学公式、代码结构、编译方法和验证方法。

### 3.2 `build/`：编译产物

```text
build/
└── bin/
    └── explicitAdvectionFoamStudent
```

这里存放编译后的学生版求解器。编译脚本会把输出放到项目自己的
`build/bin/`，不会把可执行文件写入 OpenFOAM 系统目录。

如果重新克隆仓库后 `build/bin/` 中没有可执行文件，需要重新编译：

```bash
sh scripts/build_student_solver.sh
```

### 3.3 `cases/`：OpenFOAM 可运行案例和题目数据

`cases/` 是项目的运行案例目录。每个案例对应一种题目物理问题、网格类型和空间格式，
每个案例下面再按分辨率划分为 `N10`、`N20` 等目录。

```text
cases/
├── 01_sine_wave_quad/
│   ├── N10/
│   ├── N20/
│   ├── N40/
│   └── N80/
├── 02_sine_wave_quad_linearUpwind/
├── 03_sine_wave_tri_upwind/
├── 03_sine_wave_tri_linearUpwind/
├── 04_solid_rotation_quad_upwind/
└── 04_solid_rotation_tri_upwind/
```

以一个分辨率案例为例：

```text
cases/01_sine_wave_quad/N20/
├── 0.orig/
│   ├── U
│   └── T
├── constant/
│   └── polyMesh/
├── system/
│   ├── controlDict
│   ├── fvSchemes
│   ├── fvSolution
│   └── blockMeshDict
├── Allrun
└── Allclean
```

各部分含义：

- `0.orig/U`：速度场初值；
- `0.orig/T`：被输运标量的初值；
- `constant/polyMesh/`：`blockMesh` 或 `gmshToFoam` 生成的网格；
- `system/controlDict`：起止时间、输出控制、`maxCo` 和求解器名称；
- `system/fvSchemes`：时间离散、梯度、散度和面插值格式；
- `system/fvSolution`：OpenFOAM 运行所需的求解设置；
- `system/blockMeshDict`：四边形结构化网格定义；
- `Allrun`：单个案例的网格、检查和求解运行入口；
- `Allclean`：清理该案例的运行产物。

`cases/` 是由 JSON 配置和脚本生成的运行区域。长期维护时，JSON 和脚本是参数源，
`cases/` 中的运行文件是根据参数生成的具体 OpenFOAM 实例。

### 3.4 `scripts/configs/`：案例参数入口

每个 JSON 文件描述一组具有相同物理问题、网格类型和空间格式的案例族。

例如：

```text
scripts/configs/
├── 01_sine_wave_quad_upwind.json
├── 02_sine_wave_quad_linearUpwind.json
├── 03_sine_wave_tri_upwind.json
├── 03_sine_wave_tri_linearUpwind.json
├── 04_solid_rotation_quad_upwind.json
└── 04_solid_rotation_tri_upwind.json
```

JSON 不直接存储每一个 cell 的 `T` 数值，而是描述如何生成这些数据。
脚本读取 JSON 后，根据 `problem`、`meshType`、`velocityModel` 和 `initialProfile`
生成网格、初始场和 OpenFOAM 字典。

常用字段如下：

| 字段 | 作用 |
|---|---|
| `caseName` | 生成到 `cases/<caseName>/` 的案例名称 |
| `description` | 案例说明 |
| `problem` | 物理问题类型，如正弦波或固体旋转 |
| `meshType` | `quad` 或 `tri` |
| `meshBackend` | 网格生成方式，如 `blockMesh` |
| `schemeName` | 空间格式名称 |
| `divScheme` | OpenFOAM 散度格式，如 `Gauss upwind` 或 `Gauss linearUpwind grad(T)` |
| `solver` | 使用的求解器名称 |
| `templateCaseName` | 从哪个案例模板复制字典 |
| `templateResolution` | 模板案例的分辨率 |
| `resolutions` | 默认运行的网格分辨率列表 |
| `velocity` | 常速度场 `[Ux, Uy, Uz]` |
| `velocityModel` | 非常速度场模型，如固体旋转 |
| `initialProfile` | 初始标量场形状 |
| `boundaryCondition` | 边界条件类型 |
| `endTime` | 终止时间 |
| `maxCo` | 目标最大 CFL 数 |
| `thickness` | 二维棱柱网格的厚度 |
| `postprocess` | 后处理类型和比较方式 |
| `implemented` | 是否允许脚本实际运行 |

### 3.5 `scripts/`：自动化工作流

顶层脚本的职责如下：

| 脚本 | 作用 |
|---|---|
| `build_student_solver.sh` | 编译学生版求解器 |
| `prepare_case.py` | 根据 JSON 和一个 `N` 生成单个 OpenFOAM 案例 |
| `run_case.py` | 准备、运行并后处理单个案例 |
| `postprocess_case.py` | 对已经运行的单个案例进行后处理 |
| `run_study.py` | 对多个 `N` 自动准备、运行、收集、分析和绘图 |
| `collect_results.py` | 收集多个分辨率的 `summary.json` |
| `analyze_study.py` | 计算误差和观察收敛阶 |
| `plot_study.py` | 绘制误差、收敛阶、场对比和其他结果图 |

`scripts/common/` 是可复用的 Python 模块：

- `case_config.py`：读取和校验 JSON；
- `foam_case.py`：生成、运行和后处理 OpenFOAM case；
- `foam_fields.py`：读写 OpenFOAM 字段；
- `gmsh_tri_mesh.py`：生成三角形网格；
- `advection_sine.py`：正弦波初值和精确解；
- `advection_rotation.py`：固体旋转速度场和初始轮廓；
- `metrics.py`：误差、质量和场范围计算；
- `study_analysis.py`：多网格数据收集、收敛分析和绘图；
- `plotting.py`：通用绘图辅助函数。

### 3.6 `data/`：数值数据和分析数据

```text
data/
├── cases/
│   └── <caseName>/
│       └── N<N>/
│           ├── summary.json
│           ├── time_history.csv
│           ├── field_data.csv
│           └── error_field.csv
└── analysis/
    └── <caseName>/
        ├── raw_results.csv
        ├── convergence_summary.csv
        ├── run_manifest.json
        └── analysis.md
```

`data/cases/` 保存每个分辨率的详细结果：

- `summary.json`：该案例的总结果、误差、CFL、守恒和状态；
- `time_history.csv`：时间推进过程中 CFL、振幅、质量等监测量；
- `field_data.csv`：最终场的 cell centre、体积和数值场；
- `error_field.csv`：正弦波数值解与精确解的局部误差。

`data/analysis/` 保存一组网格研究的汇总：

- `raw_results.csv`：所有 `N` 的原始汇总；
- `convergence_summary.csv`：$L_1$、$L_2$、$L_\infty$ 和观察收敛阶；
- `run_manifest.json`：本次研究使用的配置和分辨率；
- `analysis.md`：自动生成的分析摘要。

### 3.7 `figures/`：图片结果

```text
figures/
├── cases/
│   └── <caseName>/N<N>/
│       ├── field_comparison.png
│       ├── diagonal_profile.png
│       ├── amplitude_history.png
│       ├── cfl_history.png
│       └── contour_final.png
└── analysis/
    └── <caseName>/
        ├── convergence_errors.png
        ├── convergence_order.png
        └── all_N_comparison.png
```

`figures/cases/` 用于查看一个具体分辨率的场和监测量，`figures/analysis/` 用于查看
一组分辨率的误差、收敛阶和全分辨率对比。

### 3.8 `docs/`：学习、调试和开发记录

- `docs/01/`：第一题从公式到代码的教学文档；
- `docs/compare/`：不同格式、网格和运行方式的对比说明；
- `docs/bug_log.md`：Bug、根因、修复和验证记录；
- `docs/changelog.md`：每轮重大改动记录；
- `docs/changelog_intro.md`：改动记录格式模板。

### 3.9 `report/`：正式实验报告

```text
report/
└── 01/
    ├── report.md
    └── evidence_index.md
```

`report/01/report.md` 是第一题的正式报告，包含问题定义、数学离散、案例设置、
实验表格、图片、图号引用、误差分析、网格和格式比较、结论及复现命令。

### 3.10 `pdf/`：题目、推导和参考资料

- `pdf/training_examples_incomp.pdf`：整组训练题的总题目资料；
- `pdf/01/`：第一题题目、有限体积推导和自包含解答；
- `pdf/02/`：第二题扩散方程题目与推导；
- `pdf/03/`：第三题对流扩散方程题目与推导；
- `pdf/tex/`：可编辑的 LaTeX 源文件；
- `pdf/参考资料/`：OpenFOAM、有限体积法和 C++ 参考书。

## 4. 从 GitHub 下载并复现项目

下面以 Ubuntu 22.04/24.04、OpenFOAM 14 和 Python 3 为例。

### 4.1 下载仓库

```bash
git clone https://github.com/1a776/CFD.git
cd CFD/slover/student_project
```

如果仓库默认目录结构与当前版本不同，先确认当前目录中存在：

```bash
ls
```

应能看到：

```text
README.md  UDF  cases  data  docs  figures  pdf  report  scripts
```

### 4.2 加载 OpenFOAM 环境

```bash
source /opt/openfoam14/etc/bashrc
```

可以检查关键命令：

```bash
command -v wmake
command -v blockMesh
command -v checkMesh
python3 --version
```

如果 OpenFOAM 安装在其他路径，可以临时指定：

```bash
export OPENFOAM_BASHRC=/your/path/to/openfoam/etc/bashrc
source "$OPENFOAM_BASHRC"
```

运行脚本时也可以显式传入：

```bash
python3 scripts/run_study.py \
    --config scripts/configs/01_sine_wave_quad_upwind.json \
    --bashrc "$OPENFOAM_BASHRC"
```

### 4.3 编译学生版求解器

```bash
sh scripts/build_student_solver.sh
```

检查编译产物：

```bash
ls -l build/bin/explicitAdvectionFoamStudent
build/bin/explicitAdvectionFoamStudent -help
```

编译脚本只修改项目自己的 `build/bin/`，不会修改 `/opt/openfoam14/src`。

### 4.4 先准备一个案例

建议第一次先只准备，不立即运行：

```bash
python3 scripts/run_study.py \
    --config scripts/configs/01_sine_wave_quad_upwind.json \
    --resolutions 10 \
    --prepare-only \
    --overwrite
```

检查生成的案例：

```bash
find cases/01_sine_wave_quad/N10 -maxdepth 2 -type f | sort
```

检查网格：

```bash
source /opt/openfoam14/etc/bashrc
cd cases/01_sine_wave_quad/N10
blockMesh
checkMesh
cd ../../..
```

通常正式运行时不需要手动执行这些步骤，因为 `run_study.py` 会自动准备案例。

### 4.5 一键运行一组网格研究

以第一题正弦波、四边形网格、一阶迎风为例：

```bash
python3 scripts/run_study.py \
    --config scripts/configs/01_sine_wave_quad_upwind.json \
    --resolutions 10,20,40,80 \
    --overwrite
```

这个命令会依次完成：

1. 根据 JSON 生成 `cases/01_sine_wave_quad/N10` 到 `N80`；
2. 生成或导入网格；
3. 写入初始速度场和标量场；
4. 运行 `checkMesh`；
5. 运行 `explicitAdvectionFoamStudent`；
6. 生成 `data/cases/` 中的单案例结果；
7. 收集为 `data/analysis/` 汇总表；
8. 计算误差和观察收敛阶；
9. 生成 `figures/analysis/` 和 `figures/cases/` 图片。

### 4.6 分步运行一个案例

如果希望观察每一步，可以使用单案例脚本：

```bash
python3 scripts/prepare_case.py \
    --config scripts/configs/01_sine_wave_quad_upwind.json \
    --N 40 \
    --overwrite
```

然后运行：

```bash
python3 scripts/run_case.py \
    --config scripts/configs/01_sine_wave_quad_upwind.json \
    --N 40 \
    --no-prepare
```

只重新生成后处理：

```bash
python3 scripts/postprocess_case.py \
    --config scripts/configs/01_sine_wave_quad_upwind.json \
    --N 40
```

### 4.7 单独重新收集、分析和绘图

当案例已经运行完成，只想重新生成汇总和图片时：

```bash
python3 scripts/collect_results.py \
    --config scripts/configs/01_sine_wave_quad_upwind.json \
    --resolutions 10,20,40,80

python3 scripts/analyze_study.py \
    --config scripts/configs/01_sine_wave_quad_upwind.json

python3 scripts/plot_study.py \
    --config scripts/configs/01_sine_wave_quad_upwind.json
```

## 5. 通过 JSON 修改案例

本项目的基本操作方式是：

```text
修改 scripts/configs/*.json
    -> run_study.py 读取 JSON
    -> 生成 cases/<caseName>/N<N>
    -> 运行 OpenFOAM
    -> 写入 data/
    -> 生成 figures/
```

大多数同一物理问题下的网格、分辨率、速度、终止时间和插值格式变化，
可以直接通过 JSON 完成，不需要手动逐个修改每个 `cases/Nxx/` 文件。

### 5.1 修改网格分辨率

修改：

```json
"resolutions": [10, 20, 40, 80]
```

例如只研究 `N=20,40,80`：

```json
"resolutions": [20, 40, 80]
```

也可以不修改 JSON，直接用命令行覆盖：

```bash
python3 scripts/run_study.py \
    --config scripts/configs/01_sine_wave_quad_upwind.json \
    --resolutions 20,40,80 \
    --overwrite
```

命令行的 `--resolutions` 优先于 JSON 中的 `resolutions`。

### 5.2 修改四边形和三角形网格

四边形案例：

```json
"meshType": "quad",
"meshBackend": "blockMesh"
```

三角形案例：

```json
"meshType": "tri"
```

不要只把已有 JSON 的 `meshType` 随意改成 `tri` 后就认为所有情况都能运行。
三角形案例还需要对应的三角形模板、网格生成器和字段初始化逻辑。
最稳妥的做法是复制已有配置作为新配置，再修改 `caseName`：

```bash
cp scripts/configs/01_sine_wave_quad_upwind.json \
   scripts/configs/my_sine_wave_tri_upwind.json
```

然后在新文件中修改：

```json
"caseName": "my_sine_wave_tri_upwind",
"meshType": "tri",
"divScheme": "Gauss upwind"
```

如果该网格类型和问题组合尚未在 `scripts/common/foam_case.py`、
`scripts/common/gmsh_tri_mesh.py` 或相关初值函数中实现，仅修改 JSON 不足以完成新案例。

### 5.3 修改一阶迎风和线性迎风

一阶迎风：

```json
"schemeName": "upwind",
"divScheme": "Gauss upwind"
```

线性迎风：

```json
"schemeName": "linearUpwind",
"divScheme": "Gauss linearUpwind grad(T)"
```

运行四边形线性迎风：

```bash
python3 scripts/run_study.py \
    --config scripts/configs/02_sine_wave_quad_linearUpwind.json \
    --resolutions 10,20,40,80 \
    --overwrite
```

运行三角形线性迎风：

```bash
python3 scripts/run_study.py \
    --config scripts/configs/03_sine_wave_tri_linearUpwind.json \
    --resolutions 10,20,40,80 \
    --overwrite
```

### 5.4 修改速度场和终止时间

正弦波常速度场：

```json
"velocity": [1.0, 1.0, 0.0],
"endTime": 1.0
```

固体旋转速度场：

```json
"velocityModel": {
  "type": "solidRotation",
  "center": [0.5, 0.5],
  "angularVelocity": 1.0
},
"endTime": 6.283185307179586
```

对于固体旋转，`endTime` 应写成一圈的时间 $2\pi/\omega$。如果修改
`angularVelocity`，应同步修改 `endTime`，否则就不再是完整一圈。

### 5.5 修改 CFL

```json
"maxCo": 0.2
```

例如改成 `0.1`：

```json
"maxCo": 0.1
```

然后重新运行：

```bash
python3 scripts/run_study.py \
    --config scripts/configs/01_sine_wave_quad_upwind.json \
    --resolutions 10,20,40,80 \
    --overwrite
```

`maxCo` 会影响自动时间步大小和时间步数量，但不会改变空间网格的单元数量。
进行格式或网格精度比较时，应尽量保持 CFL、终止时间和误差定义一致。

### 5.6 修改初始轮廓

正弦波初值和固体旋转初值由 `problem` 和 `initialProfile` 控制。
固体旋转案例中的典型字段包括：

```json
"initialProfile": {
  "type": "slottedDiskConeCosineHump",
  "radius": 0.15,
  "diskCenter": [0.5, 0.75],
  "coneCenter": [0.5, 0.25],
  "humpCenter": [0.25, 0.5],
  "slotHalfWidth": 0.025,
  "slotTopY": 0.85
}
```

这些字段会传递到初始场生成逻辑，而不是直接写出每一个 cell 的场值。

## 6. 常用案例命令

### 正弦波：四边形 + 一阶迎风

```bash
python3 scripts/run_study.py \
    --config scripts/configs/01_sine_wave_quad_upwind.json \
    --resolutions 10,20,40,80 \
    --overwrite
```

### 正弦波：四边形 + 线性迎风

```bash
python3 scripts/run_study.py \
    --config scripts/configs/02_sine_wave_quad_linearUpwind.json \
    --resolutions 10,20,40,80 \
    --overwrite
```

### 正弦波：三角形 + 一阶迎风

```bash
python3 scripts/run_study.py \
    --config scripts/configs/03_sine_wave_tri_upwind.json \
    --resolutions 10,20,40,80 \
    --overwrite
```

### 正弦波：三角形 + 线性迎风

```bash
python3 scripts/run_study.py \
    --config scripts/configs/03_sine_wave_tri_linearUpwind.json \
    --resolutions 10,20,40,80 \
    --overwrite
```

### 固体旋转：四边形

```bash
python3 scripts/run_study.py \
    --config scripts/configs/04_solid_rotation_quad_upwind.json \
    --resolutions 50,100,200 \
    --overwrite
```

### 固体旋转：三角形

```bash
python3 scripts/run_study.py \
    --config scripts/configs/04_solid_rotation_tri_upwind.json \
    --resolutions 50,100,200 \
    --overwrite
```

## 7. 如何检查结果

运行完成后，先检查数据：

```bash
cat data/analysis/01_sine_wave_quad/convergence_summary.csv
cat data/cases/01_sine_wave_quad/N80/summary.json
```

重点检查：

- `meshOK` 是否为 `true`；
- `solverEnded` 是否为 `true`；
- `solverFatal` 是否为 `false`；
- `maxCo` 是否不超过目标值；
- $L_1$ 是否随网格加密下降；
- 观察收敛阶是否逐渐稳定；
- 质量误差是否处于浮点舍入误差量级；
- 线性迎风是否出现过冲；
- 固体旋转在 $t=2\pi$ 时轮廓是否回到大致原位置。

查看图片：

```text
figures/analysis/<caseName>/
figures/cases/<caseName>/N<N>/
```

查看第一题完整报告：

[`report/01/report.md`](report/01/report.md)

## 8. 开发与学习入口

如果你要从公式开始学习第一题，建议按下面顺序阅读：

1. [`docs/01/00_learning_path.md`](docs/01/00_learning_path.md)
2. [`docs/01/01_formula_to_code_map.md`](docs/01/01_formula_to_code_map.md)
3. [`docs/01/04_stage1_face_flux.md`](docs/01/04_stage1_face_flux.md)
4. [`docs/01/05_stage2_cfl.md`](docs/01/05_stage2_cfl.md)
5. [`docs/01/06_stage3_convection_residual.md`](docs/01/06_stage3_convection_residual.md)
6. [`docs/01/07_stage4_forward_euler_update.md`](docs/01/07_stage4_forward_euler_update.md)
7. [`docs/01/08_stage5_time_loop.md`](docs/01/08_stage5_time_loop.md)
8. [`docs/01/09_stage6_visualization_and_convergence.md`](docs/01/09_stage6_visualization_and_convergence.md)

求解器开发说明：

[`UDF/README.md`](UDF/README.md)

Bug 记录：

[`docs/bug_log.md`](docs/bug_log.md)

改动记录模板：

[`docs/changelog_intro.md`](docs/changelog_intro.md)

## 9. 当前状态与后续扩展

当前第一题已经完成教学型求解器、四边形/三角形网格案例、多分辨率实验、
误差收敛分析、固体旋转验证和报告整理。

当前仍需谨慎解释的内容：

- 三角形和四边形相同名义 $N$ 对应的实际单元数不同，不能直接作为完全公平的网格性能排名；
- 一阶迎风有明显数值耗散；
- 线性迎风误差较小，但粗网格可能产生过冲；
- `cycleL1AgainstInitial` 是固体旋转回归指标，不等同于正弦波的解析解误差；
- 当前项目重点是第一题，后续扩散和对流扩散题目需要新增相应求解器或模型。

## 10. 许可证与说明

本项目是用于学习有限体积法、OpenFOAM 源码结构和数值验证流程的学生项目。
项目代码和案例应结合题目 PDF、推导文档、运行日志和结果数据一起阅读，
不能只根据单个图片或单个误差数字判断数值方法优劣。
