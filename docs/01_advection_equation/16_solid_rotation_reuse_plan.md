# 第一题第二案例的最大复用实现方案

本文依据：

```text
/home/a776/workdocuments/上交船舶/pdf/first_advection_problem.tex
```

目标是实现第一题 1.2.2：

```text
复杂轮廓的刚体旋转问题
```

本方案先只实现现有配置中已经预留的四边形 `N100` 案例：

```text
scripts/configs/01_advection_equation/04_solid_rotation_quad_upwind.json
```

三角形网格分支可以在四边形版本跑通后再接入。这样改动最小，也最适合教学推进。

## 1. 题目数学定义

控制方程仍然是：

$$\frac{\partial\phi}{\partial t}+\nabla\cdot(\boldsymbol{u}\phi)=0.$$

速度场为：

$$\boldsymbol{u}(x,y)=(0.5-y,\;x-0.5).$$

计算区域为：

$$[0,1]^2.$$

终止时间为一圈旋转：

$$t=2\pi.$$

初始标量场为：

$$\phi_0(x,y)=\phi_D(x,y)+\phi_C(x,y)+\phi_H(x,y).$$

三个轮廓半径均为：

$$r_0=0.15.$$

### 开槽圆盘

圆心：

$$
(x_D,y_D)=(0.5,0.75)
$$

半径函数：

$$
r_D(x,y)=\sqrt{(x-0.5)^2+(y-0.75)^2}
$$

初值：

$$
\phi_D(x,y)
=
\begin{cases}
1,
&
r_D(x,y)\leq0.15
\quad\text{且}\quad
\left(|x-0.5|\geq0.025\ \text{或}\ y\geq0.85\right),\\
0,
&\text{其他}.
\end{cases}
$$

### 圆锥

圆心：

$$
(x_C,y_C)=(0.5,0.25)
$$

半径函数：

$$
r_C(x,y)=\sqrt{(x-0.5)^2+(y-0.25)^2}
$$

初值：

$$
\phi_C(x,y)
=
\begin{cases}
1-r_C(x,y)/0.15,
&r_C(x,y)\leq0.15,\\
0,
&r_C(x,y)>0.15.
\end{cases}
$$

### 光滑凸峰

圆心：

$$
(x_H,y_H)=(0.25,0.5)
$$

半径函数：

$$
r_H(x,y)=\sqrt{(x-0.25)^2+(y-0.5)^2}
$$

初值：

$$
\phi_H(x,y)
=
\begin{cases}
\dfrac14\left[1+\cos\left(\pi r_H(x,y)/0.15\right)\right],
&r_H(x,y)\leq0.15,\\
0,
&r_H(x,y)>0.15.
\end{cases}
$$

## 2. 当前项目能复用什么

这个案例和正弦波案例属于同一个 PDE，所以最应该复用的是求解器主体：

```text
UDF/solver/01_advection_equation/explicitAdvectionFoamStudent/explicitAdvectionFoamStudent.C
```

它已经完成：

```text
读取 U 和 T
    -> fvc::flux(U) 生成面通量
        -> 根据 maxCo 计算时间步
            -> fvc::div(phi,T,"div(phi,T)") 计算显式残差
                -> 前向 Euler 更新
                    -> 时间循环写出 T
```

这些都不需要因为“正弦波”变成“复杂轮廓”而重写。原因是求解器只关心：

| 求解器需要 | 正弦波案例 | 刚体旋转案例 |
|---|---|---|
| `U` | 常速度场 | 空间变速度场 |
| `T` | 正弦初值 | 复杂轮廓初值 |
| `fvSchemes` | `Gauss upwind` 或 `linearUpwind` | 仍可用 `Gauss upwind` |
| `controlDict` | `endTime=1` | `endTime=2*pi` |

也就是说，求解器只要求磁盘上有合法的 `0/U` 和 `0/T`。至于这些场来自哪个公式，应由 Python 生成。

## 3. 为什么不能只改 JSON

现有 JSON 已经有：

```json
{
  "problem": "solid_rotation_advection",
  "velocity": [0.0, 0.0, 0.0],
  "endTime": 6.283185307179586,
  "resolutions": [100],
  "implemented": false
}
```

但它还不能直接运行，原因有三点。

第一，当前速度写入函数只支持常速度：

```text
scripts/common/foam_fields.py
    patch_uniform_vector_field()
```

它只能把：

```foam
internalField uniform (1 1 0);
```

替换成另一个常向量。刚体旋转需要：

$$\boldsymbol{u}(x,y)=(0.5-y,\;x-0.5),$$

所以每个 cell 的 `U` 都不同，必须写成 `nonuniform List<vector>`。

第二，当前初值分支只支持正弦波：

```text
scripts/common/foam_case.py
    _write_initial_field()
        -> sine_wave_advection
```

`scripts/common/advection_rotation.py` 目前还是占位：

```text
raise NotImplementedError
```

第三，当前后处理默认把 `quad` 案例当作正弦波平移，计算正弦波精确解。刚体旋转案例的验收目标是：

```text
t = 2*pi 后的标量场等值线图
```

所以后处理需要新分支，不能继续调用正弦波精确解。

结论是：JSON 可以也应该表达“函数和边界条件的选择”，但真正计算函数值的代码必须放在 Python 里。JSON 负责声明，Python 负责执行。

## 4. 推荐架构图

```mermaid
flowchart TD
    A[04_solid_rotation_quad_upwind.json] --> B[CaseConfig]
    B --> C[prepare_case.py]
    C --> D[foam_case.prepare_case]
    D --> E[blockMeshDict: 生成 N100 四边形网格]
    D --> F[advection_rotation.py: 写复杂轮廓 T]
    D --> G[foam_fields.py: 写旋转速度 U]
    D --> H[controlDict: endTime=2*pi, maxCo=0.2]
    D --> I[fvSchemes: div(phi,T)=Gauss upwind]
    E --> J[Allrun]
    F --> J
    G --> J
    H --> J
    I --> J
    J --> K[blockMesh]
    K --> L[checkMesh]
    L --> M[explicitAdvectionFoamStudent]
    M --> N[time directories: final T at t=2*pi]
    N --> O[postprocess_case.py solid_rotation branch]
    O --> P[figures/01_advection_equation/cases/04_solid_rotation_quad_upwind/N100/contour_final.png]
    O --> Q[data/01_advection_equation/cases/04_solid_rotation_quad_upwind/N100/summary.json]
```

如果 Mermaid 不能渲染，可以按下面的文本图理解：

```text
JSON 配置
    -> Python 解析 CaseConfig
        -> 生成 OpenFOAM case
            -> blockMesh 生成网格
            -> Python 写 U(x,y) 和 T0(x,y)
            -> 求解器推进到 2*pi
            -> Python 画最终等值线
```

## 5. 推荐 JSON 表达方式

当前 JSON 不建议继续只写：

```json
"velocity": [0.0, 0.0, 0.0]
```

因为这会误导读者，以为速度是零。更好的方式是保留可读、可扩展的声明：

```json
{
  "caseName": "04_solid_rotation_quad_upwind",
  "problem": "solid_rotation_advection",
  "meshType": "quad",
  "meshBackend": "blockMesh",
  "schemeName": "upwind",
  "divScheme": "Gauss upwind",
  "solver": "explicitAdvectionFoamStudent",
  "templateCaseName": "01_sine_wave_quad",
  "templateResolution": 20,
  "resolutions": [100],
  "domain": [0.0, 1.0, 0.0, 1.0],
  "resolutionMeaning": "cellsPerEdge",
  "velocityModel": "solidRotationAboutCenter",
  "rotationCenter": [0.5, 0.5],
  "angularVelocity": 1.0,
  "initialProfile": "slottedDiskConeCosineHump",
  "profileRadius": 0.15,
  "slotHalfWidth": 0.025,
  "slotTopY": 0.85,
  "endTime": 6.283185307179586,
  "maxCo": 0.2,
  "thickness": 0.1,
  "boundaryCondition": "zeroScalarAtOuterBoundary",
  "implemented": true
}
```

这里的关键原则是：

```text
JSON 写模型名字和参数，不写大量逐单元数值。
Python 根据模型名字和参数写 OpenFOAM 字段。
```

这样未来扩展第二题、第三题时，也能继续沿用：

```text
problem + model + parameters
```

而不是每个案例都临时硬编码。

## 6. 边界条件建议

正弦波平移使用周期边界，是因为题目明确要求：

$$\phi(0,y,t)=\phi(1,y,t),\qquad \phi(x,0,t)=\phi(x,1,t).$$

刚体旋转案例不同。它的三个轮廓都在区域内部旋转，中心绕 $(0.5,0.5)$ 转动，最大外缘不会碰到边界：

```text
中心轨道半径约 0.25
轮廓半径 0.15
最外距离约 0.40
离边界仍有约 0.10 的余量
```

因此建议四个外边界使用普通 `patch`，标量 `T` 使用零值边界：

```foam
T:
    xMin/xMax/yMin/yMax -> fixedValue uniform 0
    zMin/zMax           -> empty
```

速度 `U` 应在边界上也表达旋转速度场。最简单的可执行版本可以先让边界 `U` 使用 `zeroGradient`，因为被输运标量在边界附近为零，主目标是内部轮廓旋转图。但更严谨的版本应写入边界面中心处的旋转速度：

$$\boldsymbol{u}(x,y)=(0.5-y,\;x-0.5).$$

本项目第一版建议：

```text
先实现内部 cell 的非均匀 U；
T 边界使用 fixedValue 0；
U 边界使用 zeroGradient；
验证场不接触边界后，再考虑写精确边界 U。
```

这样实现最小，且对本题一圈后等值线图影响很小。

## 7. 分辨率语义

题面写的是：

```text
每条边布置 100 个顶点
```

现有项目中 `N` 的语义是：

```text
每条边 N 个单元
```

这两个概念相差 1：

```text
100 个顶点 -> 99 个单元
100 个单元 -> 101 个顶点
```

为了最大复用现有脚本，建议第一版采用：

```text
N100 = 每条边 100 个单元
```

原因是现有 `run_case.py --N 100`、`patch_block_mesh_resolution()`、输出目录 `N100` 都是按 cell 数理解的。如果严格按题面顶点数，应另开字段：

```json
"edgeVertices": 100
```

并把 blockMesh 单元数写成 `99 99 1`。这会让当前目录命名和脚本语义变复杂，不适合第一版。

## 8. 需要新增或改造的能力

### 8.1 `advection_rotation.py`

当前文件是占位。应扩展为专门负责刚体旋转案例的公式：

```text
solid_rotation_velocity(x,y)
solid_rotation_profile(x,y)
write_case_initial_field(case, nx, ny, params)
write_case_velocity_field(case, nx, ny, params)
```

其中：

$$u=0.5-y,\qquad v=x-0.5.$$

初值函数按 `first_advection_problem.tex` 的三段公式实现。

### 8.2 `foam_fields.py`

当前只有：

```text
write_scalar_field()
patch_uniform_vector_field()
```

需要补一个非均匀向量场写出函数：

```text
write_vector_field(path, values, boundary_types, object_name="U", location="0")
```

这样 `U` 可以写成：

```foam
internalField nonuniform List<vector>
...
```

### 8.3 `foam_case.py`

`_write_initial_field()` 当前只支持：

```text
sine_wave_advection
```

需要加入：

```text
solid_rotation_advection
```

同时 `_patch_velocity_field()` 不能再只按 `velocity` 常量写。推荐改成：

```text
if velocityModel exists:
    调用模型函数写 U
else:
    维持当前 uniform velocity 逻辑
```

这能保证已有正弦波案例不受影响。

### 8.4 `postprocess_case.py`

当前 `quad` 后处理默认计算正弦波精确解和误差。刚体旋转需要新增分支：

```text
if problem == "solid_rotation_advection":
    读取初始 T 和最终 T
    计算质量变化、极值、时间步、CFL
    绘制初始/最终/差值或最终等值线图
```

建议输出：

```text
figures/01_advection_equation/cases/04_solid_rotation_quad_upwind/N100/contour_final.png
figures/01_advection_equation/cases/04_solid_rotation_quad_upwind/N100/field_comparison.png
figures/01_advection_equation/cases/04_solid_rotation_quad_upwind/N100/cfl_history.png
data/01_advection_equation/cases/04_solid_rotation_quad_upwind/N100/summary.json
data/01_advection_equation/cases/04_solid_rotation_quad_upwind/N100/field_data.csv
```

由于一圈旋转后理论上图形回到初始位置，也可以额外计算：

$$
L_1^{cycle}
=
\frac{\sum_cV_c|T_c(t=2\pi)-T_c(0)|}{\sum_cV_c|T_c(0)|}.
$$

但题面只要求画等值线图，所以第一版不要把误差作为硬性验收。

## 9. 推荐实施顺序

第一步，只准备配置层：

```text
把 04 JSON 从 planned 改成 implemented=true
补充 velocityModel、initialProfile、boundaryCondition 等声明字段
```

第二步，写数学函数层：

```text
advection_rotation.py
    -> 速度函数
    -> 三个轮廓函数
    -> 初值写入函数
```

第三步，写 OpenFOAM 场输出层：

```text
foam_fields.py
    -> write_vector_field()
```

第四步，接入 case 生成：

```text
foam_case.py
    -> solid_rotation_advection 分支
    -> 非均匀 U 写入
    -> 非周期边界 patch 和 fixedValue 0 标量边界
```

第五步，接入后处理：

```text
postprocess_case.py
    -> solid_rotation_advection 分支
    -> 最终等值线图
    -> summary.json
```

第六步，最小验证：

```bash
source /opt/openfoam14/etc/bashrc
sh scripts/build_student_solver.sh
python3 scripts/run_case.py \
    --config scripts/configs/01_advection_equation/04_solid_rotation_quad_upwind.json \
    --N 100 \
    --overwrite
```

## 10. 验收标准

第一版四边形 `N100` 案例跑通后，应满足：

| 检查项 | 合格标准 |
|---|---|
| case 准备 | `cases/01_advection_equation/04_solid_rotation_quad_upwind/N100/` 生成成功 |
| 网格 | `log.checkMesh` 包含 `Mesh OK` |
| 求解器 | `log.explicitAdvectionFoamStudent` 包含 `Stage 5 time loop completed` 和 `End` |
| 时间 | 最终时间接近 `6.283185307179586` |
| CFL | `maxCo` 不超过目标值附近 |
| 初值 | 三个轮廓在正确位置出现 |
| 最终图 | 生成一圈后的标量场等值线图 |
| 数值稳定性 | 最终场没有 `nan`、`inf`，极值在合理范围内 |
| 质量 | 质量变化有记录，作为数值耗散/边界影响判断依据 |

## 11. 对五道题总目标的判断

你的终极目标是跑完 PDF 五道题。这个目标不能靠不断堆临时脚本完成，应该逐步形成统一结构：

```text
problem
    -> equation family
    -> initial condition model
    -> velocity model
    -> boundary condition model
    -> mesh type
    -> solver
    -> postprocess objective
```

第一题第二案例正好是一个分水岭：它逼着项目从“常速度 + 正弦初值”升级到“模型化速度 + 模型化初值 + 按问题分支后处理”。这一步做好后，后续第二题扩散、第三题对流扩散、第四题 Poisson、第五题 Navier-Stokes 都可以沿用同一种配置思想。

因此，当前最合理的路线不是马上重构整个项目，而是：

```text
在现有 JSON 驱动框架中加入少量模型字段；
在 Python 中实现具体公式；
保持 C++ 求解器只负责数值推进；
让后处理按 problem 分支输出题目要求的图或误差。
```

这能最大限度复用现有代码，也能为五道题保留清晰扩展路径。
