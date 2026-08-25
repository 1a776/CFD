# 第一题参数试验、JSON 配置和运行命令

本文说明如何在当前统一脚本框架下改变：

```text
网格分辨率 N
四边形/三角形网格
一阶迎风/linearUpwind
第一案例/第二案例
```

所有命令从项目根目录执行：

```bash
cd /home/a776/workdocuments/上交船舶/slover/student_project
```

## 1. 运行前的固定准备

### 1.1 加载 OpenFOAM

四边形和三角形实际运行都需要：

```bash
source /opt/openfoam14/etc/bashrc
```

### 1.2 编译学生版求解器

如果源码刚修改过，先编译：

```bash
sh scripts/build_student_solver.sh
```

检查可执行文件：

```bash
test -x build/01_advection_equation/bin/explicitAdvectionFoamStudent
echo $?
```

输出 `0` 表示文件存在且可执行。

### 1.3 三角形案例的 Gmsh Python

三角形案例还需要能导入 `gmsh` 的 Python：

```bash
export VIBEFLOW_PYTHON=/home/a776/vibeflow/python-env/bin/python
```

也可以检查：

```bash
$VIBEFLOW_PYTHON -c "import gmsh; print(gmsh.__file__)"
```

## 2. JSON 的核心字段

### 2.1 案例名称

```json
"caseName": "01_sine_wave_quad"
```

决定输出目录：

```text
cases/01_advection_equation/01_sine_wave_quad/
data/01_advection_equation/cases/01_sine_wave_quad/
figures/01_advection_equation/cases/01_sine_wave_quad/
data/01_advection_equation/analysis/01_sine_wave_quad/
figures/01_advection_equation/analysis/01_sine_wave_quad/
```

如果改变 `caseName`，就是一个新的案例族。

### 2.2 数学问题

第一案例：

```json
"problem": "sine_wave_advection"
```

第二案例：

```json
"problem": "solid_rotation_advection"
```

这个字段不是普通标签，它决定 `foam_case.py` 调用哪一组初值和后处理函数。

### 2.3 网格类型

四边形：

```json
"meshType": "quad",
"meshBackend": "blockMesh"
```

三角形：

```json
"meshType": "tri",
"meshBackend": "gmsh",
"templateCaseName": "01_sine_wave_quad",
"templateResolution": 20
```

四边形走：

```text
blockMeshDict -> blockMesh
```

三角形走：

```text
gmsh_tri_mesh.py
    -> mesh.msh
    -> gmshToFoam
    -> createPatch
```

### 2.4 对流格式

一阶迎风：

```json
"schemeName": "upwind",
"divScheme": "Gauss upwind"
```

`linearUpwind`：

```json
"schemeName": "linearUpwind",
"divScheme": "Gauss linearUpwind grad(T)",
"gradTScheme": "Gauss linear"
```

真正写入的是：

```foam
divSchemes
{
    div(phi,T) Gauss upwind;
}
```

或：

```foam
gradSchemes
{
    grad(T) Gauss linear;
}

divSchemes
{
    div(phi,T) Gauss linearUpwind grad(T);
}
```

`schemeName` 主要用于识别；`divScheme` 和 `gradTScheme` 才是 OpenFOAM
实际读取的内容。

### 2.5 分辨率

JSON 默认值：

```json
"resolutions": [10, 20, 40, 80]
```

命令行可以覆盖 JSON：

```bash
--resolutions 10,20,40,80
```

命令行优先级高于 JSON。

## 3. 第一案例的四种组合

当前配置文件：

| 组合 | 配置 |
|---|---|
| quad + upwind | `01_sine_wave_quad_upwind.json` |
| quad + linearUpwind | `02_sine_wave_quad_linearUpwind.json` |
| tri + upwind | `03_sine_wave_tri_upwind.json` |
| tri + linearUpwind | `03_sine_wave_tri_linearUpwind.json` |

### 3.1 quad + upwind

```bash
source /opt/openfoam14/etc/bashrc

python3 scripts/run_study.py \
    --config scripts/configs/01_advection_equation/01_sine_wave_quad_upwind.json \
    --resolutions 10,20,40,80 \
    --overwrite
```

流程：

```text
准备 N10/N20/N40/N80
    -> blockMesh
    -> checkMesh
    -> explicitAdvectionFoamStudent
    -> 单 N 后处理
    -> 收集结果
    -> 收敛阶分析
    -> 总体绘图
```

### 3.2 quad + linearUpwind

```bash
source /opt/openfoam14/etc/bashrc

python3 scripts/run_study.py \
    --config scripts/configs/01_advection_equation/02_sine_wave_quad_linearUpwind.json \
    --resolutions 10,20,40,80 \
    --overwrite
```

### 3.3 tri + upwind

```bash
source /opt/openfoam14/etc/bashrc
export VIBEFLOW_PYTHON=/home/a776/vibeflow/python-env/bin/python

python3 scripts/run_study.py \
    --config scripts/configs/01_advection_equation/03_sine_wave_tri_upwind.json \
    --resolutions 10,20,40,80 \
    --overwrite
```

### 3.4 tri + linearUpwind

```bash
source /opt/openfoam14/etc/bashrc
export VIBEFLOW_PYTHON=/home/a776/vibeflow/python-env/bin/python

python3 scripts/run_study.py \
    --config scripts/configs/01_advection_equation/03_sine_wave_tri_linearUpwind.json \
    --resolutions 10,20,40,80 \
    --overwrite
```

## 4. 只修改一个案例的 N 分布

例如，想让四边形正弦波只跑：

```text
N=20,40,80,160
```

有两种方式。

### 4.1 临时修改，不改 JSON

推荐先这样做：

```bash
python3 scripts/run_study.py \
    --config scripts/configs/01_advection_equation/01_sine_wave_quad_upwind.json \
    --resolutions 20,40,80,160 \
    --prepare-only \
    --overwrite
```

确认目录和字典正确后再实际运行：

```bash
source /opt/openfoam14/etc/bashrc

python3 scripts/run_study.py \
    --config scripts/configs/01_advection_equation/01_sine_wave_quad_upwind.json \
    --resolutions 20,40,80,160 \
    --overwrite
```

### 4.2 永久修改 JSON

把：

```json
"resolutions": [10, 20, 40, 80]
```

改成：

```json
"resolutions": [20, 40, 80, 160]
```

以后不传 `--resolutions` 时就会运行这组 N：

```bash
python3 scripts/run_study.py \
    --config scripts/configs/01_advection_equation/01_sine_wave_quad_upwind.json \
    --overwrite
```

## 5. 第二案例：当前四边形多 N

当前第二案例四边形配置已经包含：

```json
"problem": "solid_rotation_advection",
"meshType": "quad",
"resolutions": [50, 100, 200]
```

如果想跑：

```text
N=50,100,200
```

不改 JSON，直接：

```bash
source /opt/openfoam14/etc/bashrc

python3 scripts/run_study.py \
    --config scripts/configs/01_advection_equation/04_solid_rotation_quad_upwind.json \
    --resolutions 50,100,200 \
    --overwrite
```

如果希望保存为默认研究方案，修改 JSON：

```json
"resolutions": [50, 100, 200]
```

再运行：

```bash
python3 scripts/run_study.py \
    --config scripts/configs/01_advection_equation/04_solid_rotation_quad_upwind.json \
    --overwrite
```

输出目录：

```text
cases/01_advection_equation/04_solid_rotation_quad_upwind/N50/
cases/01_advection_equation/04_solid_rotation_quad_upwind/N100/
cases/01_advection_equation/04_solid_rotation_quad_upwind/N200/
```

第二案例的 `N` 含义是每条边的 cell 数：

```text
N=100 -> 100*100 个二维四边形 cell
```

## 6. 第二案例切换到 `linearUpwind`

当前没有单独的第二案例 `linearUpwind` JSON。可以复制：

```text
scripts/configs/01_advection_equation/04_solid_rotation_quad_upwind.json
```

改名为：

```text
scripts/configs/01_advection_equation/04_solid_rotation_quad_linearUpwind.json
```

并修改：

```json
{
  "caseName": "04_solid_rotation_quad_linearUpwind",
  "problem": "solid_rotation_advection",
  "meshType": "quad",
  "meshBackend": "blockMesh",
  "templateCaseName": "01_sine_wave_quad",
  "schemeName": "linearUpwind",
  "divScheme": "Gauss linearUpwind grad(T)",
  "gradTScheme": "Gauss linear",
  "resolutions": [50, 100, 150],
  "maxCo": 0.2,
  "endTime": 6.283185307179586,
  "implemented": true
}
```

其它字段：

```text
velocityModel
initialProfile
boundaryCondition
domain
```

保持原样。

试运行前先只准备一个 N：

```bash
python3 scripts/run_study.py \
    --config scripts/configs/01_advection_equation/04_solid_rotation_quad_linearUpwind.json \
    --resolutions 50 \
    --prepare-only \
    --overwrite
```

然后检查：

```bash
cat cases/01_advection_equation/04_solid_rotation_quad_linearUpwind/N50/system/fvSchemes
```

必须看到：

```foam
grad(T)         Gauss linear;
div(phi,T)      Gauss linearUpwind grad(T);
```

### 6.1 当前限制

这个配置修改在“脚本生成和字典写入层面”是合理的，但第二案例
`linearUpwind` 还需要实际运行验证。不能因为 JSON 能被解析，就直接认为
复杂轮廓经过一圈旋转后的结果已经可信。

## 7. 第二案例切换到三角形网格

题目最终要求第二案例也做三角形和四边形。目标 JSON 形式可以写成：

```json
{
  "caseName": "04_solid_rotation_tri_upwind",
  "problem": "solid_rotation_advection",
  "meshType": "tri",
  "meshBackend": "gmsh",
  "templateCaseName": "01_sine_wave_quad",
  "templateResolution": 20,
  "schemeName": "upwind",
  "divScheme": "Gauss upwind",
  "resolutions": [50, 100, 150],
  "endTime": 6.283185307179586,
  "maxCo": 0.2,
  "thickness": 0.1,
  "implemented": true
}
```

但是，**当前版本不能直接运行这个 JSON**。原因是：

```python
foam_case.py
    -> solid_rotation_advection + meshType=tri
    -> NotImplementedError
```

要真正支持第二案例三角形，至少需要扩展：

1. `advection_rotation.py`
   - 按真实三角形 cell centre 写旋转速度 `U`；
   - 按真实 cell centre 写复杂轮廓 `T`。
2. `foam_case.py`
   - 不再把旋转案例的 tri 分支直接拒绝；
   - 生成非周期外边界，而不是沿用正弦波的周期 patch。
3. `postprocess_case.py`
   - 使用真实 `C`、`Vc` 和三角形拓扑；
   - 绘制非结构三角形等值线。
4. `createPatchDict`
   - 旋转案例的 `xMin/xMax/yMin/yMax` 应是普通外边界；
   - `zMin/zMax` 仍为 `empty`。

所以当前正确的认识是：

```text
第一案例 tri
    -> 当前可以只改 JSON 运行

第二案例 tri
    -> JSON 形式可以先设计
    -> 但还必须先扩展代码
    -> 不能仅靠改 JSON 运行
```

## 8. 第二案例四边形和三角形的研究矩阵

题目要求可以整理成：

| 物理案例 | 网格 | 格式 | N |
|---|---|---|---|
| 正弦波平移 | quad | upwind | 10,20,40,80 |
| 正弦波平移 | quad | linearUpwind | 10,20,40,80 |
| 正弦波平移 | tri | upwind | 10,20,40,80 |
| 正弦波平移 | tri | linearUpwind | 10,20,40,80 |
| 刚体旋转 | quad | upwind | 100 |
| 刚体旋转 | quad | linearUpwind | 100 |
| 刚体旋转 | tri | upwind | 100 |
| 刚体旋转 | tri | linearUpwind | 100 |

当前已直接支持的行：

```text
前四行
刚体旋转 + quad + upwind
```

其余旋转组合需要完成相应代码扩展和最小验证。

## 9. 只准备、不运行

任何新组合都建议先：

```bash
python3 scripts/run_study.py \
    --config scripts/configs/01_advection_equation/<config>.json \
    --resolutions 50,100 \
    --prepare-only \
    --overwrite
```

检查：

```bash
cat cases/<caseName>/N50/metadata.json
cat cases/<caseName>/N50/system/controlDict
cat cases/<caseName>/N50/system/fvSchemes
```

对于 quad，还可以检查：

```bash
grep -n "hex" cases/<caseName>/N50/system/blockMeshDict
```

对于 tri，准备阶段不会生成最终三角形网格；必须运行 `Allrun` 后检查：

```bash
grep -n "cells:" cases/<caseName>/N50/log.checkMesh
```

## 10. 运行一个 N

```bash
source /opt/openfoam14/etc/bashrc

python3 scripts/run_case.py \
    --config scripts/configs/01_advection_equation/<config>.json \
    --N 100 \
    --overwrite
```

这个命令包含：

```text
prepare
    -> Allrun
        -> 生成网格
        -> checkMesh
        -> 求解器
    -> Python 后处理
```

## 11. 运行多个 N

```bash
source /opt/openfoam14/etc/bashrc

python3 scripts/run_study.py \
    --config scripts/configs/01_advection_equation/<config>.json \
    --resolutions 50,100,150 \
    --overwrite
```

对于第一案例，运行完成后还会自动产生：

```text
data/01_advection_equation/analysis/<caseName>/raw_results.csv
data/01_advection_equation/analysis/<caseName>/convergence_summary.csv
data/01_advection_equation/analysis/<caseName>/analysis.md
figures/01_advection_equation/analysis/<caseName>/convergence_errors.png
figures/01_advection_equation/analysis/<caseName>/convergence_order.png
```

对于当前第二案例，`run_study.py` 会完成每个 N 的单案例后处理，但不会调用
第一案例专用的 `collect/analyse/plot` 收敛分析。这是因为第二案例当前主要求
是最终等值线图，而不是正弦波解析解误差。

### 11.1 第二案例三角形网格

三角形旋转配置为：

```text
scripts/configs/01_advection_equation/04_solid_rotation_tri_upwind.json
```

当前三角形配置默认也测试 `N=50,100,200`：

```json
"problem": "solid_rotation_advection",
"meshType": "tri",
"resolutions": [50, 100, 200]
```

运行三角形的三个网格：

```bash
source /opt/openfoam14/etc/bashrc

python3 scripts/run_study.py \
    --config scripts/configs/01_advection_equation/04_solid_rotation_tri_upwind.json \
    --resolutions 50,100,200 \
    --overwrite
```

执行流程：

```text
Gmsh 生成三角形棱柱网格
    -> gmshToFoam
    -> createPatch
    -> xMin/xMax/yMin/yMax 生成普通外边界
    -> foamPostProcess 写 constant/C 和 constant/Vc
    -> 根据真实 cell centre 生成 U 和 T
    -> 运行到 t=2*pi
    -> 三角形后处理和绘图
```

结果位置：

```text
cases/01_advection_equation/04_solid_rotation_tri_upwind/N100/
cases/01_advection_equation/04_solid_rotation_tri_upwind/N200/
data/01_advection_equation/cases/04_solid_rotation_tri_upwind/N100/
data/01_advection_equation/cases/04_solid_rotation_tri_upwind/N200/
figures/01_advection_equation/cases/04_solid_rotation_tri_upwind/N100/
figures/01_advection_equation/cases/04_solid_rotation_tri_upwind/N200/
```

主要图片：

```text
figures/01_advection_equation/cases/04_solid_rotation_tri_upwind/N100/field_comparison.png
figures/01_advection_equation/cases/04_solid_rotation_tri_upwind/N100/contour_final.png
figures/01_advection_equation/cases/04_solid_rotation_tri_upwind/N100/cfl_history.png
```

三角形案例不能把场值强行 reshape 成 `N x N`。旋转三角形后处理使用：

```text
constant/C
    -> 真实 cell centre

constant/Vc
    -> 真实 cell 体积

mesh/mesh_geometry.json
    -> 三角形节点和连接关系
```

因此初始场、旋转速度、周期后的体积加权差异和图像都基于实际三角形网格。

## 12. 只做后处理

已有单个 N 的求解结果时：

```bash
python3 scripts/postprocess_case.py \
    --config scripts/configs/01_advection_equation/01_sine_wave_quad_upwind.json \
    --N 40
```

已有第一案例全部 N 结果时：

```bash
python3 scripts/collect_results.py \
    --config scripts/configs/01_advection_equation/01_sine_wave_quad_upwind.json \
    --resolutions 10,20,40,80

python3 scripts/analyze_study.py \
    --config scripts/configs/01_advection_equation/01_sine_wave_quad_upwind.json

python3 scripts/plot_study.py \
    --config scripts/configs/01_advection_equation/01_sine_wave_quad_upwind.json
```

## 13. 验收标准

### 13.1 所有案例通用

```text
checkMesh 输出 Mesh OK
solver log 没有 FOAM FATAL ERROR
solver log 出现 Stage 5 time loop completed
最终时间与 controlDict/endTime 一致
maxCo 不超过目标值
```

### 13.2 第一案例

```text
0/T 的单元数与网格 cell 数一致
最终场能计算 L1
N 加密时误差总体下降
convergence_summary.csv 中能看到 L1order
```

### 13.3 第二案例

```text
初值中出现圆盘、圆锥和光滑凸峰
最终时间为 2*pi
最终目录中存在 T
生成 contour_final.png
质量变化接近 0
场值没有 nan 或 inf
```

## 14. 推荐的实际工作顺序

```text
1. 复制或选择 JSON
2. 只修改 caseName、meshType、schemeName、divScheme、resolutions
3. --prepare-only
4. 检查 system/fvSchemes 和 system/controlDict
5. 先运行一个 N
6. 查看 log.checkMesh 和单案例 summary.json
7. 再运行全部 N
8. 最后做收敛分析和绘图
```

不要一开始就同时改变：

```text
物理问题 + 网格类型 + 空间格式 + N + 边界条件
```

否则出现错误时很难判断是哪个层次造成的。
