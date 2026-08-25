# 学生版线性对流 CFD 工程

本项目用于学习并实现 PDF 第 1 题：

$$\frac{\partial \phi}{\partial t}+\nabla\cdot(\boldsymbol{u}\phi)=0.$$

当前已经实现的是二维周期正弦波算例：

- 结构化四边形网格；
- Gmsh 生成的三角形棱柱网格；
- 显式时间推进；
- `Gauss upwind` 和 `Gauss linearUpwind grad(T)` 两种空间格式；
- `N10/N20/N40/N80` 多网格误差与收敛阶分析。

四边形和三角形网格、两种对流格式都已经接入 JSON 驱动流程。
固体旋转算例已经接入四边形和三角形两种网格。需要注意：
配置标记为 `"implemented": true` 表示脚本链路已经可以生成、运行和后处理
该配置，不等于所有网格分辨率和格式组合都已经完成物理验证。

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
│   ├── 02_sine_wave_quad_linearUpwind/
│       ├── N10/
│       ├── N20/
│       ├── N40/
│       └── N80/
│   ├── 03_sine_wave_tri_upwind/
│       ├── N10/
│       ├── N20/
│       ├── N40/
│       └── N80/
│   └── 03_sine_wave_tri_linearUpwind/
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
│   │   ├── 03_sine_wave_tri_linearUpwind.json
│   │   ├── 04_solid_rotation_quad_upwind.json
│   │   └── 04_solid_rotation_tri_upwind.json
│   └── common/
│       ├── paths.py
│       ├── case_config.py
│       ├── foam_case.py
│       ├── mesh_tools.py
│       ├── advection_sine.py
│       ├── advection_rotation.py
│       ├── foam_fields.py
│       ├── gmsh_tri_mesh.py
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
    ├── compare/
    └── bug_log.md
```

旧的 `scripts/01_sine_wave_quad/` 和 `scripts/02_sine_wave_quad_linearUpwind/` 已删除。现在所有案例统一通过顶层脚本和 JSON 配置运行。

开发、运行和后处理过程中遇到的 Bug、日志证据、解决方式和验证标准统一记录在：

```text
docs/bug_log.md
```

## 配置文件

同一个题目内部，不同网格、插值格式、速度、终止时间和网格分辨率都通过
JSON 配置文件区分。脚本入口不直接写死案例名称，而是先读取 JSON，
再把 JSON 转换成 OpenFOAM 案例目录和字典。

### 逐字段说明

| 字段 | 必需性 | 当前作用 | 写入或影响的位置 |
|---|---|---|---|
| `caseName` | 必需 | 案例族名称 | `cases/<caseName>/`、`data/<...>`、`figures/<...>` |
| `description` | 可选 | 给人看的案例说明 | 主要保留在 JSON 中 |
| `problem` | 必需 | 选择物理问题分支 | 初值和后处理分派；当前支持 `sine_wave_advection` |
| `meshType` | 必需 | 选择网格流程 | `quad` 走 `blockMesh`，`tri` 走 Gmsh |
| `meshBackend` | 可选 | 记录网格后端 | quad 默认 `blockMesh`；tri 当前使用 `gmsh` |
| `templateCaseName` | tri 推荐填写 | 指定基础案例族 | 与 `templateResolution` 共同定位模板 |
| `templateResolution` | 可选 | 指定基础模板的 N | 默认 `20`，例如 `cases/01_sine_wave_quad/N20/` |
| `gmshPython` | tri 可选 | 指定 Gmsh Python 解释器 | 未填写时使用 `/home/a776/vibeflow/python-env/bin/python` |
| `schemeName` | 可选 | 给人看的格式名称 | 写入 `metadata.json`，便于识别 |
| `divScheme` | 必需 | 对流散度离散格式 | 写入 `system/fvSchemes` 的 `div(phi,T)` |
| `gradTScheme` | linearUpwind 必需 | `grad(T)` 的梯度格式 | 写入 `system/fvSchemes` 的 `gradSchemes` |
| `solver` | 可选 | 指定求解器名称 | `controlDict.application`、`Allrun` |
| `resolutions` | 可选 | 默认的 N 列表 | `run_study.py` 未传 `--resolutions` 时使用 |
| `velocity` | 可选 | 常速度 $(u_x,u_y,u_z)$ | 写入 `0.orig/U`，也用于精确解 |
| `endTime` | 可选 | 终止时间 | 写入 `system/controlDict` |
| `maxCo` | 可选 | 目标 Courant 数 | 写入 `system/controlDict`，由求解器读取 |
| `thickness` | tri 可选 | 三角形棱柱厚度 | 传给 `gmsh_tri_mesh.py`，默认 `0.1` |
| `implemented` | 可选 | 是否允许脚本执行 | `false` 时在准备阶段直接拒绝 |

其中，`schemeName` 只是标签；真正决定 OpenFOAM 离散格式的是
`divScheme`，而 `linearUpwind` 还需要 `gradTScheme`。

### 当前配置总表

| 配置文件 | `caseName` | 网格 | 网格后端 | 对流格式 | 当前状态 |
|---|---|---|---|---|---|
| `01_sine_wave_quad_upwind.json` | `01_sine_wave_quad` | 四边形 | `blockMesh` | `Gauss upwind` | 已实现 |
| `02_sine_wave_quad_linearUpwind.json` | `02_sine_wave_quad_linearUpwind` | 四边形 | `blockMesh` | `Gauss linearUpwind grad(T)` | 已实现 |
| `03_sine_wave_tri_upwind.json` | `03_sine_wave_tri_upwind` | 三角形棱柱 | `gmsh` | `Gauss upwind` | 已实现 |
| `03_sine_wave_tri_linearUpwind.json` | `03_sine_wave_tri_linearUpwind` | 三角形棱柱 | `gmsh` | `Gauss linearUpwind grad(T)` | 已接入，需完成完整收敛验证 |
| `04_solid_rotation_quad_upwind.json` | `04_solid_rotation_quad_upwind` | 四边形 | `blockMesh` | `Gauss upwind` | 已接入 |
| `04_solid_rotation_tri_upwind.json` | `04_solid_rotation_tri_upwind` | 三角形棱柱 | `gmsh` | `Gauss upwind` | 已接入 |

当前可运行配置：

```text
scripts/configs/01_sine_wave_quad_upwind.json
scripts/configs/02_sine_wave_quad_linearUpwind.json
scripts/configs/03_sine_wave_tri_upwind.json
scripts/configs/03_sine_wave_tri_linearUpwind.json
```

旋转配置：

```text
scripts/configs/04_solid_rotation_quad_upwind.json
scripts/configs/04_solid_rotation_tri_upwind.json
```

两个旋转配置都已经接入统一脚本流程。

### 配置之间的关系

对于任意一个配置，目录关系可以按下面的链条理解：

```text
scripts/configs/<config>.json
    |
    v
CaseConfig
    |
    +-- caseName --------------------> cases/<caseName>/
    |
    +-- templateCaseName + N --------> 基础模板
    |
    +-- resolutions -----------------> N10/N20/N40/N80
    |
    +-- meshType --------------------> quad 或 tri 的网格分支
    |
    +-- divScheme/gradTScheme -------> system/fvSchemes
    |
    +-- endTime/maxCo ---------------> system/controlDict
    |
    +-- velocity --------------------> 0.orig/U 和精确解
    |
    +-- solver ----------------------> controlDict.application/Allrun
    |
    +-- problem ---------------------> 初值和后处理分支
```

当前最重要的关系是：

```text
meshType = quad
    -> 从模板复制
    -> 修改 blockMeshDict 的 (N N 1)
    -> blockMesh

meshType = tri
    -> 从四边形基础模板复制 0.orig/system/constant 输入
    -> 删除 blockMeshDict
    -> 写入 createPatchDict
    -> Gmsh 生成三角形棱柱网格
    -> gmshToFoam 导入
```

所以三角形案例没有一个必须长期维护的静态 `tri_template/` 目录。
当前三角形配置的基础模板是：

```text
templateCaseName = 01_sine_wave_quad
templateResolution = 20
实际模板 = cases/01_sine_wave_quad/N20/
```

这个模板只提供通用的 `0.orig/`、`system/` 和 `constant/` 输入。
三角形网格本身由 `scripts/common/gmsh_tri_mesh.py` 在运行时生成。

旋转三角形案例使用同一个 Gmsh 网格生成器，但采用不同的边界和场生成
分支：

```text
boundaryCondition=zeroScalarAtOuterBoundary
    -> xMin/xMax/yMin/yMax 为普通外边界 patch
    -> T 的外边界为 fixedValue 0

solid_rotation_advection
    -> 先生成真实三角形 cell centre
    -> 再在这些 centre 上生成旋转速度 U 和初始场 T
```

## 命令总览和前置条件

所有命令建议从项目根目录执行：

```bash
cd /home/a776/workdocuments/上交船舶/slover/student_project
```

### 不同命令需要什么环境

| 操作 | Python 3 | OpenFOAM 环境 | 已编译求解器 | Gmsh Python |
|---|---:|---:|---:|---:|
| `prepare_case.py` | 是 | 否 | 否 | 否 |
| `run_study.py --prepare-only` | 是 | 否 | 否 | 否 |
| 四边形实际运行 | 是 | 是 | 是 | 否 |
| 三角形实际运行 | 是 | 是 | 是 | 是 |
| 单案例后处理 | 是 | 否 | 否 | 否 |
| 收敛表和总体绘图 | 是 | 否 | 否 | 否 |

只准备目录时，下面的命令足够：

```bash
python3 scripts/run_study.py \
    --config scripts/configs/03_sine_wave_tri_linearUpwind.json \
    --resolutions 10,20,40,80 \
    --prepare-only \
    --overwrite
```

它不会调用 Gmsh、`gmshToFoam`、`checkMesh` 或学生求解器。
对于 tri 案例，真正的 `T` 初值要等网格生成并写出 `constant/C` 后，
由 `Allrun` 按真实 cell centre 重新生成。准备阶段目录中的 `0.orig/T`
不能作为三角形最终初值的验收依据。

实际运行前需要：

```bash
source /opt/openfoam14/etc/bashrc
sh scripts/build_student_solver.sh
```

三角形 `Allrun` 默认使用：

```text
/home/a776/vibeflow/python-env/bin/python
```

如果 Gmsh Python 环境不在这个位置，可以临时指定：

```bash
export VIBEFLOW_PYTHON=/path/to/python-with-gmsh
```

### 顶层脚本参数

| 脚本 | 主要参数 | 功能 |
|---|---|---|
| `build_student_solver.sh` | 无 | 编译 UDF 求解器到 `build/bin/` |
| `prepare_case.py` | `--config`、`--N`、`--overwrite` | 准备一个 N |
| `prepare_case.py` | `--refresh-initial-only` | 网格已生成后，只重写 tri 的 `0.orig/T` |
| `run_case.py` | `--config`、`--N`、`--overwrite` | 准备、运行、后处理一个 N |
| `run_case.py` | `--no-prepare` | 已有案例目录时跳过准备 |
| `run_case.py` | `--no-postprocess` | 运行后不做 Python 后处理 |
| `run_study.py` | `--config`、`--resolutions` | 批量运行多个 N |
| `run_study.py` | `--prepare-only` | 批量准备但不运行 |
| `run_study.py` | `--bashrc` | 指定 OpenFOAM bashrc |
| `collect_results.py` | `--config`、`--resolutions` | 收集各 N 的 `summary.json` |
| `analyze_study.py` | `--config` | 计算收敛表和报告 |
| `plot_study.py` | `--config` | 绘制总体误差、阶数和振幅图 |
| `postprocess_case.py` | `--config`、`--N` | 后处理一个已运行案例 |

`--resolutions` 会覆盖 JSON 中的 `resolutions`。例如：

```bash
python3 scripts/run_study.py \
    --config scripts/configs/03_sine_wave_tri_linearUpwind.json \
    --resolutions 20,40 \
    --prepare-only \
    --overwrite
```

如果不传这个参数，则使用 JSON 中的：

```json
"resolutions": [10, 20, 40, 80]
```

脚本会自动排序、去重，并拒绝非正整数。

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

03 三角形案例也可以只准备目录。它会写出 `Allrun`、`Allclean`、
`createPatchDict` 和 `metadata.json`，但不会在准备阶段生成 Gmsh 网格：

```bash
python3 scripts/run_study.py \
    --config scripts/configs/03_sine_wave_tri_upwind.json \
    --prepare-only \
    --overwrite
```

三角形网格使用 `linearUpwind` 时，只需要换用另一个 JSON：

```bash
python3 scripts/run_study.py \
    --config scripts/configs/03_sine_wave_tri_linearUpwind.json \
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

运行并后处理 03 三角形案例的 `N20`：

```bash
python3 scripts/run_case.py \
    --config scripts/configs/03_sine_wave_tri_upwind.json \
    --N 20 \
    --overwrite
```

运行三角形网格 + `linearUpwind` 的 `N20`：

```bash
python3 scripts/run_case.py \
    --config scripts/configs/03_sine_wave_tri_linearUpwind.json \
    --N 20 \
    --overwrite
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

三角形正弦波案例使用同一个入口，只替换 JSON：

```bash
source /opt/openfoam14/etc/bashrc

python3 scripts/run_study.py \
    --config scripts/configs/03_sine_wave_tri_upwind.json \
    --resolutions 10,20,40,80 \
    --overwrite
```

三角形网格 + `linearUpwind`：

```bash
source /opt/openfoam14/etc/bashrc

python3 scripts/run_study.py \
    --config scripts/configs/03_sine_wave_tri_linearUpwind.json \
    --resolutions 10,20,40,80 \
    --overwrite
```

三角形案例每个 N 的内部流程是：

```text
Gmsh Python
    -> mesh/mesh.msh
    -> gmshToFoam
    -> createPatch
    -> checkMesh
    -> foamPostProcess(writeCellCentres/writeCellVolumes)
    -> 按真实 cell centre 写 0.orig/T
    -> explicitAdvectionFoamStudent
    -> 按真实 cell centre 和 Vc 后处理
```

### 一个生成后的 tri 案例目录

以 `03_sine_wave_tri_linearUpwind/N20` 为例：

```text
cases/03_sine_wave_tri_linearUpwind/N20/
├── 0.orig/
│   ├── U                 # JSON velocity 写入的常速度
│   └── T                 # 运行时按真实 cell centre 生成
├── 0/                    # Allrun 复制的运行初始场
├── constant/
│   ├── polyMesh/         # gmshToFoam 导入后的 OpenFOAM 网格
│   ├── C                 # foamPostProcess 写出的 cell centre
│   ├── Vc                # foamPostProcess 写出的 cell volume
│   └── ...
├── mesh/
│   ├── mesh.msh          # Gmsh 原始网格
│   └── mesh_geometry.json# Gmsh 节点和三角形连接关系
├── system/
│   ├── controlDict
│   ├── createPatchDict
│   ├── fvSchemes
│   └── fvSolution
├── Allrun
├── Allclean
├── case.foam
└── metadata.json
```

`mesh/`、`constant/polyMesh/`、`constant/C`、`constant/Vc`、`0/`、
时间目录和 `log.*` 都属于可以重新生成的运行产物。`0.orig/`、
`system/` 和 JSON 配置才是长期需要理解和维护的输入。

### 两种空间格式最终写入什么

四边形和三角形只改变网格生成与几何读取方式；对流格式由
`fvSchemes` 中的同一项控制：

```foam
divSchemes
{
    default         none;
    div(phi,T)      Gauss upwind;
}
```

或：

```foam
gradSchemes
{
    default         Gauss linear;
    grad(T)         Gauss linear;
}

divSchemes
{
    default         none;
    div(phi,T)      Gauss linearUpwind grad(T);
}
```

因此：

```text
quad + upwind
    -> blockMesh + Gauss upwind

quad + linearUpwind
    -> blockMesh + grad(T) + Gauss linearUpwind grad(T)

tri + upwind
    -> Gmsh + gmshToFoam + Gauss upwind

tri + linearUpwind
    -> Gmsh + gmshToFoam + grad(T) + Gauss linearUpwind grad(T)
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

四边形案例的 `field_data.csv` 使用 `i,j` 网格索引；三角形案例使用真实
单元几何：

```text
cell,x,y,z,volume,initial,numerical,exact,error,absError
```

三角形案例的 `volume` 直接来自 OpenFOAM 写出的 `constant/Vc`，不会
使用四边形案例的统一体积假设。三角形图使用 Gmsh 的节点和三角形连接
关系绘制。

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
- Gmsh 三角形棱柱网格；
- `Gauss upwind`；
- `Gauss linearUpwind grad(T)`；
- 三角形网格的真实 cell centre、cell volume、误差和可视化；
- 单个 N 后处理；
- 多个 N 的误差和收敛阶分析。

### 当前验证状态

| 案例 | 配置状态 | 已有结果状态 | 说明 |
|---|---|---|---|
| `01_sine_wave_quad` | `implemented=true` | `N10/N20/N40/N80` 有结果 | 四边形 + upwind |
| `02_sine_wave_quad_linearUpwind` | `implemented=true` | `N10/N20/N40/N80` 有结果 | 四边形 + linearUpwind |
| `03_sine_wave_tri_upwind` | `implemented=true` | `N10/N20/N40/N80` 有结果 | 三角形 + upwind |
| `03_sine_wave_tri_linearUpwind` | `implemented=true` | 配置和案例目录已建立；完整四级结果待补 | `N10` 已做最小运行验证 |
| `04_solid_rotation_quad_upwind` | `implemented=true` | 已接入 | 四边形固体旋转 |
| `04_solid_rotation_tri_upwind` | `implemented=true` | 已接入 | 三角形固体旋转 |

三角形 `linearUpwind` 的最小 `N10` 验证结果为：

```text
网格单元数       200
最终时间         0.999999999996
maxCo            0.2
归一化 L1 误差   2.9825743812703637e-01
```

这个结果只证明配置字典、三角形网格、求解器时间循环和后处理链路能够
闭合，不代表 `N10/N20/N40/N80` 的最终收敛阶已经得到。

### 本次故障记录：tri + linearUpwind

第一次运行三角形 `linearUpwind` 时，外层 Python 报：

```text
subprocess.CalledProcessError
```

这只是 `Allrun` 返回非零状态后的包装异常。真正的 OpenFOAM 根因在：

```text
FOAM FATAL IO ERROR:
keyword divSchemes is undefined in dictionary .../system/fvSchemes
```

生成的 `fvSchemes` 曾经错误地变成：

```foam
gradSchemes
{
    default         Gauss linear;
}
    grad(T)         Gauss linear;
}
```

多出来的关闭大括号使 `divSchemes` 不再处于顶层字典。
问题已经在 `scripts/common/foam_case.py` 的 `_patch_fv_schemes()`
中修复。新的逻辑会识别 `gradSchemes` 的完整大括号范围，并把
`grad(T)` 插入唯一的关闭大括号之前。

如果某个案例是在修复前生成的，需要重新准备：

```bash
python3 scripts/run_study.py \
    --config scripts/configs/03_sine_wave_tri_linearUpwind.json \
    --resolutions 10,20,40,80 \
    --prepare-only \
    --overwrite
```

还要注意，`-overwrite` 是 `runApplication` 等 OpenFOAM 工具脚本的参数，
不是当前学生求解器本身的参数。不要手动执行：

```bash
./build/bin/explicitAdvectionFoamStudent -overwrite
```

暂未实现：

- 固体旋转复杂剖面；
- slotted disk、cone、smooth hump 初始场；
- 非结构网格上的专用可视化；
- 第 1 题之外的扩散、Poisson 和 Navier-Stokes 求解器。

清理单个案例：

```bash
sh cases/01_sine_wave_quad/N40/Allclean
```

`Allclean` 只删除可重新生成的运行产物，不删除 `0.orig/`、`system/` 和常量输入文件。

三角形案例的清理命令：

```bash
sh cases/03_sine_wave_tri_upwind/N20/Allclean
```

tri 版除四边形运行产物外，还会清理 `mesh/`、`constant/C`、
`constant/Cc*` 和 `constant/Vc`，但会保留 `0.orig/`、`system/` 和
`system/createPatchDict`。

## 代码入口速查

```text
scripts/run_study.py
    -> 读取 JSON、解析 N 列表、批量调度

scripts/prepare_case.py
    -> 准备一个 N 的案例目录

scripts/common/case_config.py
    -> JSON -> CaseConfig

scripts/common/foam_case.py
    -> 模板复制、字典修改、Allrun/Allclean 生成、运行调度

scripts/common/gmsh_tri_mesh.py
    -> 生成二维三角形底面并挤出为三角形棱柱网格

scripts/common/advection_sine.py
    -> 正弦波初值和任意 cell centre 的精确解

scripts/common/postprocess_case.py
    -> 读取场、几何、日志，计算误差并生成单案例图片

scripts/common/study_analysis.py
    -> 汇总 N、计算观察收敛阶、生成总体图

UDF/solver/explicitAdvectionFoamStudent/
    -> 学生版 OpenFOAM 求解器源码

build/bin/explicitAdvectionFoamStudent
    -> 编译后的可执行文件
```
