# 基础配置和公式的对应关系

这份文档面向刚开始学习 OpenFOAM 和有限体积法的同学。它不急着扩展新案例，而是先解释当前工程中已经跑通的“正弦波对流案例”为什么能工作，以及每一个基础配置文件在数学公式里对应什么。

你可以把本工程理解成三层：

```text
数学题目
    -> JSON 实验配置
        -> OpenFOAM case 文件和 Python/UDF 代码
```

当前第一题的控制方程是：

$$\frac{\partial\phi}{\partial t}+\nabla\cdot(\boldsymbol{u}\phi)=0.$$

在本工程里，数学中的 $\phi$ 没有直接叫 `phi`，而是叫 `T`。这是为了避免和 OpenFOAM 里常用的面通量字段 `phi` 混淆：

| 数学量 | 工程名字 | 文件或代码位置 | 含义 |
|---|---|---|---|
| $\phi$ | `T` | `0.orig/T`、`0/T`、时间目录里的 `T` | 被输运的标量 |
| $\boldsymbol{u}$ | `U` | `0.orig/U`、`0/U` | 给定速度场 |
| $F_f=\boldsymbol{u}_f\cdot\boldsymbol{S}_f$ | `phi` | C++ 中的 `surfaceScalarField phi` | 面体积通量 |
| $\Omega$ | `mesh` | `constant/polyMesh` | 计算区域和网格 |
| $t$ | `runTime` | C++ 中的 `Time` 对象 | 时间、输出和目录管理 |

## 1. JSON 是案例的第一层开关

当前项目不建议你直接手改每个 `cases/<caseName>/Nxx/` 目录。更合理的入口是：

```text
scripts/configs/<config>.json
```

以 `01_sine_wave_quad_upwind.json` 为例，它描述的是：

```json
{
  "caseName": "01_sine_wave_quad",
  "problem": "sine_wave_advection",
  "meshType": "quad",
  "schemeName": "upwind",
  "divScheme": "Gauss upwind",
  "solver": "explicitAdvectionFoamStudent",
  "resolutions": [10, 20, 40, 80],
  "velocity": [1.0, 1.0, 0.0],
  "endTime": 1.0,
  "maxCo": 0.2,
  "thickness": 0.1,
  "implemented": true
}
```

每个字段都应该和公式对应起来看：

| JSON 字段 | 公式或数值实验含义 | 最终影响 |
|---|---|---|
| `caseName` | 这组实验的名字 | 输出到 `cases/`、`data/`、`figures/` |
| `problem` | 选择哪一个数学问题 | 决定初值、精确解和后处理分支 |
| `meshType` | 区域 $\Omega$ 怎样离散 | `quad` 用 `blockMesh`，`tri` 用 Gmsh |
| `divScheme` | $\nabla\cdot(\boldsymbol{u}\phi)$ 的空间离散格式 | 写入 `system/fvSchemes` 的 `div(phi,T)` |
| `solver` | 用哪个求解器推进方程 | 写入 `controlDict.application` 和 `Allrun` |
| `resolutions` | 网格加密序列 | 生成 `N10/N20/N40/N80` |
| `velocity` | $\boldsymbol{u}$ | 写入 `0.orig/U`，也用于正弦波精确解 |
| `endTime` | 终止时间 $T_{\mathrm{end}}$ | 写入 `system/controlDict` |
| `maxCo` | 目标 CFL 数 | 写入 `system/controlDict`，求解器读取 |
| `thickness` | 二维问题的计算厚度 | 影响体积 $V_c$ 和通量量级 |
| `implemented` | 是否允许脚本执行 | `false` 时脚本拒绝运行 |

这里最重要的判断是：JSON 不是 OpenFOAM 原生格式，它是本项目自己的“实验说明书”。Python 脚本读 JSON，然后生成 OpenFOAM 能读懂的字典和字段文件。

## 2. `case_config.py` 把 JSON 变成 Python 对象

代码入口是：

```text
scripts/common/case_config.py
```

它里面的 `CaseConfig` 可以理解成“已经解析好的 JSON”。例如：

```python
case_name: str
problem: str
mesh_type: str
div_scheme: str
velocity: tuple[float, float, float]
end_time: float
max_co: float
```

数学上，这一步没有做离散，也没有求解。它只是把实验参数从文本变成程序可以访问的变量。

可以把它理解成：

```text
JSON 文件
    -> CaseConfig
        -> 后续生成 case 目录
```

如果 JSON 中写：

```json
"maxCo": 0.2
```

那么 `CaseConfig.max_co` 就是 `0.2`，后面会被写进：

```foam
maxCo           0.2;
```

再由 C++ 求解器读取，进入 CFL 公式：

$$\Delta t=\frac{2\,\mathrm{Co}_{\mathrm{target}}}{\max_c\left(\sum_f|F_{cf}|/V_c\right)}.$$

## 3. `prepare_case.py` 负责生成一个 OpenFOAM 算例

命令形态是：

```bash
python3 scripts/prepare_case.py \
    --config scripts/configs/01_sine_wave_quad_upwind.json \
    --N 40 \
    --overwrite
```

它的数学意义是：选择一个分辨率 $N$，把连续区域 $[0,1]^2$ 离散成 $N\times N$ 个控制体，并准备初值和运行字典。

主要调用链是：

```text
prepare_case.py
    -> load_config()
    -> prepare_case()
        -> patch_block_mesh_resolution()
        -> _patch_fv_schemes()
        -> _patch_control_dict()
        -> _patch_velocity_field()
        -> _write_initial_field()
        -> _write_case_scripts()
```

这些函数和公式的关系如下：

| 函数 | 数学或工程含义 |
|---|---|
| `patch_block_mesh_resolution()` | 设置 $N_x=N_y=N$ |
| `_patch_fv_schemes()` | 设置 $\nabla\cdot(\boldsymbol{u}\phi)$ 的离散格式 |
| `_patch_control_dict()` | 设置求解器、终止时间和 CFL |
| `_patch_velocity_field()` | 写入速度场 $\boldsymbol{u}$ |
| `_write_initial_field()` | 写入初始标量场 $\phi_0(x,y)$ |
| `_write_case_scripts()` | 写出 `Allrun` 和 `Allclean` |

所以 `prepare_case.py` 做的是“把数学实验翻译成 OpenFOAM 文件”，不是求解。

## 4. `blockMeshDict` 对应计算区域和网格

四边形案例中，网格由：

```text
system/blockMeshDict
```

生成。当前基础区域是：

$$[0,1]^2.$$

OpenFOAM 仍然用三维网格表示二维问题，所以实际顶点是：

```foam
(0 0 0)
(1 0 0)
(1 1 0)
(0 1 0)
(0 0 0.1)
(1 0 0.1)
(1 1 0.1)
(0 1 0.1)
```

这里的 `0.1` 是厚度。由于 `zMin` 和 `zMax` 是 `empty`，求解仍然是二维的。

网格单元数由这一行决定：

```foam
hex (0 1 2 3 4 5 6 7) (N N 1) simpleGrading (1 1 1)
```

数学上就是：

$$\Omega=[0,1]\times[0,1],\qquad \Omega=\bigcup_c\Omega_c.$$

当 `N=40` 时：

```text
单元数 = 40 * 40 * 1 = 1600
```

如果厚度是 $h_z=0.1$，则每个四边形柱状控制体的体积近似为：

$$V_c=\frac{1}{N}\cdot\frac{1}{N}\cdot 0.1.$$

## 5. 边界条件对应物理边界

正弦波平移案例的题目要求是：

```text
x 和 y 方向周期边界
```

所以 `blockMeshDict` 里有：

```foam
xMin
{
    type cyclic;
    neighbourPatch xMax;
}
xMax
{
    type cyclic;
    neighbourPatch xMin;
}
```

数学上，这对应：

$$\phi(0,y,t)=\phi(1,y,t),\qquad \phi(x,0,t)=\phi(x,1,t).$$

`0.orig/T` 和 `0.orig/U` 里也必须写同样的 patch 类型：

```foam
boundaryField
{
    xMin { type cyclic; }
    xMax { type cyclic; }
    yMin { type cyclic; }
    yMax { type cyclic; }
    zMin { type empty; }
    zMax { type empty; }
}
```

这说明：网格边界类型和字段边界类型必须互相匹配。网格说 `cyclic`，字段也要说 `cyclic`；网格说 `empty`，字段也要说 `empty`。

## 6. `0.orig/U` 对应速度场

正弦波案例的速度是：

$$\boldsymbol{u}=(1,1).$$

在 OpenFOAM 里写成三维向量：

```foam
internalField   uniform (1 1 0);
```

这说明每个 cell 的速度都相同：

$$\boldsymbol{U}_c=(1,1,0).$$

当前 `_patch_velocity_field()` 只能替换这种 `uniform` 速度。因此它适合正弦波平移案例，不足以直接表达刚体旋转速度：

$$\boldsymbol{u}(x,y)=(0.5-y,\;x-0.5).$$

刚体旋转速度每个 cell 都不同，后续需要写成：

```foam
internalField   nonuniform List<vector>
...
```

这就是为什么第二个案例不是“只改 JSON 数字”就能完全实现。

## 7. `0.orig/T` 对应初始标量场

正弦波案例的初值是：

$$\phi_0(x,y)=\sin 2\pi(x+y).$$

在代码中，这个公式由：

```text
scripts/common/advection_sine.py
scripts/common/advection_tools.py
```

写入 `0.orig/T`。

对四边形结构网格，单元中心是：

$$x_i=\frac{i+0.5}{N_x},\qquad y_j=\frac{j+0.5}{N_y}.$$

所以写入的每个数是：

$$T_{ij}^0=\sin 2\pi(x_i+y_j).$$

这和 `0.orig/T` 中的：

```foam
internalField   nonuniform List<scalar>
...
```

一一对应。`nonuniform List<scalar>` 的意思是：每个 cell 一个标量值。

## 8. `fvSchemes` 对应空间离散格式

控制方程中的核心空间项是：

$$\nabla\cdot(\boldsymbol{u}\phi).$$

有限体积积分后是：

$$V_c\frac{\mathrm d\phi_c}{\mathrm dt}+\sum_{f\in\partial\Omega_c}F_{cf}\phi_f=0.$$

其中：

$$F_{cf}=\boldsymbol{u}_f\cdot\boldsymbol{S}_{cf}.$$

在 OpenFOAM 里，这个离散格式由 `system/fvSchemes` 决定：

```foam
divSchemes
{
    default         none;
    div(phi,T)      Gauss upwind;
}
```

这里要特别注意两个 `phi`：

| 名字 | 含义 |
|---|---|
| 数学中的 $\phi$ | 被输运标量，本工程叫 `T` |
| OpenFOAM 中的 `phi` | 面通量 $F_f$ |

所以：

```foam
div(phi,T)
```

读作：

```text
用面通量 phi 去计算标量 T 的守恒型散度
```

如果使用二阶线性迎风，则 JSON 写：

```json
"divScheme": "Gauss linearUpwind grad(T)",
"gradTScheme": "Gauss linear"
```

最终字典写成：

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

## 9. `controlDict` 对应时间推进

`system/controlDict` 中最重要的是：

```foam
application     explicitAdvectionFoamStudent;
startTime       0;
endTime         1;
deltaT          0.01;
writeControl    timeStep;
writeInterval   10;
maxCo           0.2;
velocityField   U;
advectedField   T;
```

它们和数学关系如下：

| 条目 | 数学意义 |
|---|---|
| `application` | 选择求解器 |
| `startTime` | 初始时间 $t_0$ |
| `endTime` | 终止时间 $T_{\mathrm{end}}$ |
| `deltaT` | OpenFOAM 字典里的默认时间步 |
| `writeInterval` | 每隔多少步写一次结果 |
| `maxCo` | 目标 CFL 数 |
| `velocityField` | 速度场文件名 |
| `advectedField` | 被输运标量文件名 |

当前求解器会重新计算 CFL 时间步，所以 `deltaT 0.01` 主要是为了保持标准 `controlDict` 结构。真正用于推进的是 C++ 中根据 `maxCo` 算出的 `deltaT`。

## 10. C++ 求解器对应半离散公式

有限体积半离散形式是：

$$V_c\frac{\mathrm dT_c}{\mathrm dt}+\sum_fF_{cf}T_f=0.$$

除以 $V_c$：

$$\frac{\mathrm dT_c}{\mathrm dt}=-\frac{1}{V_c}\sum_fF_{cf}T_f.$$

定义残差：

$$R_c=\frac{1}{V_c}\sum_fF_{cf}T_f.$$

则：

$$\frac{\mathrm dT_c}{\mathrm dt}=-R_c.$$

C++ 中分三步实现。

第一步，读入速度和标量：

```cpp
volVectorField U(...);
volScalarField T(...);
```

第二步，计算面通量：

```cpp
surfaceScalarField phi(..., fvc::flux(U));
```

对应：

$$F_f=\boldsymbol{u}_f\cdot\boldsymbol{S}_f.$$

第三步，计算显式散度残差：

```cpp
tmp<volScalarField> tResidual
(
    fvc::div(phi, T, "div(phi,T)")
);
```

对应：

$$R_c^n=\frac{1}{V_c}\sum_fF_{cf}T_f^n.$$

第四步，前向 Euler 更新：

```cpp
T = T - deltaTDim*residual;
T.correctBoundaryConditions();
```

对应：

$$T_c^{n+1}=T_c^n-\Delta tR_c^n.$$

## 11. CFL 公式为什么有一个 `2`

当前求解器使用：

$$\Delta t=\frac{2\,\mathrm{maxCo}}{\max_c\left(\sum_f|F_{cf}|/V_c\right)}.$$

它来自：

$$\mathrm{Co}_c=\frac{\Delta t}{2}\frac{\sum_f|F_{cf}|}{V_c}.$$

对二维四边形网格和常速度 $(1,1)$，每个 cell 的东西南北四个面都会贡献通量大小。用 $\sum_f|F_{cf}|$ 统计所有出入通量时，会把两个方向的进出面都算进去，所以公式中使用 $\frac{1}{2}$ 与常见的单方向 CFL 表达保持一致。

C++ 里对应：

```cpp
scalarField sumPhi
(
    fvc::surfaceSum(mag(phi))().primitiveField()
);

scalarField rate
(
    sumPhi/mesh.V().primitiveField()
);

const scalar rateMax = gMax(rate);
const scalar deltaT = 2.0*maxCo/rateMax;
```

这一段的数学路线是：

```text
面通量 F_f
    -> |F_f|
        -> 每个 cell 求和
            -> 除以 V_c
                -> 全域最大值
                    -> 反推出 deltaT
```

## 12. `Allrun` 是把准备、网格、求解串起来

对于四边形案例，`Allrun` 的核心流程是：

```text
Allclean
    -> refresh 0.orig/T
    -> copy 0.orig to 0
    -> blockMesh
    -> checkMesh
    -> explicitAdvectionFoamStudent
```

数学上对应：

```text
清理旧结果
    -> 写初始条件
    -> 生成控制体集合
    -> 检查几何拓扑
    -> 推进离散方程
```

`Allrun` 不是数值方法本身，它是实验流程脚本。

## 13. 后处理对应误差定义

正弦波案例有精确解：

$$\phi(x,y,t)=\sin 2\pi(x+y-(u+v)t).$$

所以 `postprocess_case.py` 可以在最终时刻计算：

$$L_1=\frac{\sum_cV_c|T_c-T_c^{exact}|}{\sum_cV_c|T_c^{exact}|}.$$

这对应输出：

```text
data/cases/<caseName>/Nxx/summary.json
data/cases/<caseName>/Nxx/field_data.csv
figures/cases/<caseName>/Nxx/field_comparison.png
```

注意：这个后处理当前是为正弦波平移案例写的。刚体旋转第二案例的题目要求主要是画一圈后的等值线图，而不是给收敛阶。因此第二案例需要新的后处理分支，不能继续把它当作正弦波来计算误差。

## 14. 对第二案例最重要的启发

第一题第二个案例仍然解同一个对流方程：

$$\frac{\partial\phi}{\partial t}+\nabla\cdot(\boldsymbol{u}\phi)=0.$$

因此 C++ 求解器的核心离散过程可以复用：

```text
U, T
    -> fvc::flux(U)
        -> fvc::div(phi,T)
            -> forward Euler
                -> write T
```

但第二案例改变了三件事：

| 内容 | 正弦波案例 | 刚体旋转案例 |
|---|---|---|
| 初值 | $\sin 2\pi(x+y)$ | 开槽圆盘 + 圆锥 + 光滑凸峰 |
| 速度 | 常速度 $(1,1)$ | 空间变速度 $(0.5-y,x-0.5)$ |
| 边界 | 周期边界 | 紧支撑轮廓在区域内旋转，建议用零值标量边界 |

所以第二案例不是重写求解器，而是扩展“配置表达能力、场生成能力、后处理能力”。这就是最大复用现有代码的关键。
