# 学生版 UDF 开发记录

## 1. 开发边界

本目录只服务于第 1 题的守恒型线性对流方程：

$$
\frac{\partial T}{\partial t}
+
\nabla\cdot(\boldsymbol{U}T)=0.
$$

目标是让学生从有限体积公式出发，逐步写出一个 OpenFOAM 14 自定义显式求解器。

本轮不实现：

- 扩散方程；
- Poisson 方程；
- Navier--Stokes 方程；
- 高阶限制器；
- 并行分解；
- 复杂几何和工程船舶算例。

不得修改 `/opt/openfoam14/src`。所有学生代码放在本目录的 `UDF/` 下。

## 2. 参考对象

开发时参考以下本机对象：

- OpenFOAM 14 的 `applications/solvers` 求解器入口；
- OpenFOAM 14 的 `src/finiteVolume` 有限体积场和离散算子；
- 上级目录 `slover/from_scratch/` 中的完整教学参考实现。

完整参考实现只用于核对概念和 API，不作为第一步复制对象。

## 3. 数学模型

对控制体 $\Omega_c$ 积分：

$$
V_c\frac{\mathrm dT_c}{\mathrm dt}
+
\sum_{f\in\partial\Omega_c}F_{cf}T_f=0,
\qquad
F_{cf}=\boldsymbol{U}_f\cdot\boldsymbol{S}_{cf}.
$$

前向 Euler：

$$
T_c^{n+1}
=T_c^n
-\frac{\Delta t}{V_c}
\sum_fF_{cf}^nT_f^n.
$$

一阶迎风面值：

$$
T_f=
\begin{cases}
T_{\mathrm{owner}},&F_f\geq0,\\
T_{\mathrm{neighbour}},&F_f<0.
\end{cases}
$$

OpenFOAM 的 `fvSchemes` 负责声明面插值格式，学生求解器负责组织时间推进和
显式更新。这个边界划分会在 `docs/01_formula_to_code_map.md` 中逐项解释。

## 4. 文件职责

| 文件 | 职责 |
|---|---|
| `explicitAdvectionFoamStudent.C` | 求解器入口和逐阶段实现位置 |
| `Make/files` | 指定源文件和可执行文件名 |
| `Make/options` | 指定头文件目录和链接库 |
| `cases/01_sine_wave_quad/system/fvSchemes` | 选择 `Gauss upwind` |
| `cases/01_sine_wave_quad/system/controlDict` | 时间范围、输出频率和 `maxCo` |
| `cases/01_sine_wave_quad/system/blockMeshDict` | 生成二维周期方形网格 |
| `cases/01_sine_wave_quad/scripts/create_initial_fields.py` | 写入正弦波初始场 |

## 5. 预定实现阶段

### Stage 0：入口和字段读取

当前已搭好。程序读取 `fvMesh`、`U` 和 `T`，打印单元数、时间和字段范围。

### Stage 1：面通量

已实现。根据

$$
F_f=\boldsymbol{U}_f\cdot\boldsymbol{S}_f
$$

创建 `surfaceScalarField phi`，并输出所有面通量中的最小值和最大值。
实现位置为：

```text
UDF/solver/explicitAdvectionFoamStudent/explicitAdvectionFoamStudent.C
```

关键接口为：

```cpp
surfaceScalarField phi
(
    IOobject
    (
        "phi",
        runTime.name(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    fvc::flux(U)
);
```

`phi` 是面标量场，因为每个面对应一个标量体积通量；它由速度场和网格面面积向量计算得到，
不是从 `0/phi` 读取的初始场。

### Stage 2：CFL 时间步

教学设计和公式说明见：

```text
docs/05_stage2_cfl.md
```

本阶段代码已经写入学生版求解器。需要计算每个单元的

$$
\frac{\sum_f|F_{cf}|}{V_c}
$$

并根据 `maxCo` 得到 $\Delta t$：

$$
\Delta t
=\frac{2\,\mathrm{maxCo}}
{\max_c\left(\sum_f|F_{cf}|/V_c\right)}.
$$

OpenFOAM 对应的底层接口为：

```cpp
fvc::surfaceSum(mag(phi))
mesh.V()
gMax(rate)
```

对应源码位置：

```text
UDF/solver/explicitAdvectionFoamStudent/explicitAdvectionFoamStudent.C
```

核心代码为：

```cpp
const scalar maxCo
(
    controlDict.lookupOrDefault<scalar>("maxCo", 0.2)
);

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

数据类型和数学位置依次为：

```text
surfaceScalarField phi
    -> scalarField sumPhi
        -> scalarField rate
            -> scalar rateMax
                -> scalar deltaT
```

`fvc::surfaceSum(mag(phi))` 先把面通量绝对值累加到每个 cell；
`()` 访问 `tmp` 内部的临时场；`.primitiveField()` 取得 cell 内部数值；
`gMax(rate)` 取得全网格最大 rate。

当前案例的手算预期为：

```text
sumPhi per cell = 0.02
rateMax         = 80
maxCo           = 0.2
deltaT          = 0.005
```

这一阶段仍然不更新 `T`，不调用 `fvc::div`，也不推进时间。

### Stage 3：显式对流残差

Stage 3 的源码已经补入求解器，但它仍然只计算和打印残差，不执行时间推进。
目标是调用 OpenFOAM 的显式散度接口，并让 `fvSchemes` 中的
`div(phi,T) Gauss upwind;` 生效。

目标接口是：

```cpp
#include "fvcDiv.H"

tmp<volScalarField> tResidual
(
    fvc::div(phi, T, "div(phi,T)")
);

const volScalarField& residual = tResidual();
```

它对应的是已经除以控制体体积的显式对流残差：

$$
R_c^n=\frac{1}{V_c}\sum_fF_{cf}T_f^n.
$$

如果手动装配 owner/neighbour 面通量，可以把未除体积的中间累加量记为：

$$
\mathcal Q_c^n=\sum_fF_{cf}T_f^n.
$$

因此：

$$
R_c^n=\frac{\mathcal Q_c^n}{V_c}.
$$

`fvc::div` 直接对应最终残差 $R_c$，不需要再手动除以 $V_c$。

如果你自己写面循环，那么中间累加量应记为 $\mathcal Q_c$，更新时再用
$T_c^{n+1}=T_c^n-\frac{\Delta t}{V_c}\mathcal Q_c^n$。

详细解释见：

```text
docs/06_stage3_convection_residual.md
```

### Stage 4：前向 Euler 更新

Stage 4 已经写入求解器，但只做一个时间步，不进入完整时间循环。

源码使用：

```cpp
const dimensionedScalar deltaTDim("deltaT", dimTime, deltaT);

T = T - deltaTDim*residual;
T.correctBoundaryConditions();
```

它对应：

$$
T_c^{n+1}=T_c^n-\Delta tR_c^n,
\qquad
R_c^n=\frac{1}{V_c}\sum_fF_{cf}^nT_f^n.
$$

因为 `residual` 已经是体积归一化残差，所以 Stage 4 不再额外除以 `mesh.V()`。
`dimensionedScalar` 给 `deltaT` 补上时间量纲，使
`deltaTDim*residual` 能够和无量纲的 `T` 相减。
更新内部场后，`T.correctBoundaryConditions()` 同步周期等边界 patch。

当前末尾日志是：

```text
Stage 4 forward Euler update completed.
No time loop has been implemented yet.
```

### Stage 5：时间循环和输出

需要你控制终止时间、设置时间步、修正边界、输出字段和日志。

时间循环中的普通 `runTime.write()` 仍然遵守 `writeControl` 和 `writeInterval`。
因此循环结束后，求解器会额外调用：

```cpp
const bool finalWriteOK = runTime.writeNow();
```

这一步强制写出终止时间的最终场，避免最后一个缩短时间步没有恰好落在
`writeInterval` 上时，后处理只能读取上一个输出目录。案例脚本也使用高精度格式
写入 `controlDict/endTime`，并提高 `timePrecision`，因此像旋转案例的
$t=2\pi$ 不会在终止时间或时间目录名中被提前截断。

### Stage 6：守恒和误差检查

最后加入体积分、归一化质量误差和正弦波 L1 误差。

## 6. 编译和检查

```bash
source /opt/openfoam14/etc/bashrc
cd /home/a776/workdocuments/上交船舶/slover/student_project
sh scripts/build_student_solver.sh
build/bin/explicitAdvectionFoamStudent -help
```

算例检查：

```bash
source /opt/openfoam14/etc/bashrc
cd /home/a776/workdocuments/上交船舶/slover/student_project/cases/01_sine_wave_quad
sh Allrun
```

## 7. 当前状态和风险

当前源码已经写到 Stage 4 单步更新，但 Stage 4 本轮尚未由本代理重新编译和运行。
因此目前可以确认“代码和教学链条已经写入”，不能把本轮标记为“Stage 4 已验收”。

需要特别注意：

- `fvc::div` 的空间格式由 `fvSchemes` 决定；
- `fvm::div` 是隐式矩阵离散，不符合本题的全显式练习目标；
- `0.orig/` 是初始场源文件，`0/` 是运行时副本；
- `0.orig/U` 不能删除；
- `Allclean` 不能误删 `0.orig/`。
- `maxCo` 必须为正数；
- `rateMax` 过小或为零时不能直接计算 CFL 时间步；
- `fvc::surfaceSum` 在本阶段用于通量绝对值求和，不等同于 Stage 3 的有向散度。
- Stage 4 只做一次 `T` 更新，不会调用 `runTime++` 或 `runTime.write()`。

## 8. 当前代码状态

当前学生版源文件已经包含：

- Stage 0：读取 `runTime`、`fvMesh`、`U` 和 `T`；
- Stage 1：用 `fvc::flux(U)` 构造面体积通量 `phi`；
- Stage 2：根据 `maxCo`、`fvc::surfaceSum(mag(phi))` 和 `mesh.V()` 计算候选
  CFL 时间步。
- Stage 3：用 `fvc::div(phi, T, "div(phi,T)")` 计算体积归一化对流残差；
- Stage 4：用前向 Euler 对 `T` 做一次显式更新并修正边界条件。

Stage 2 的执行边界是：

```text
面通量
    -> 面通量绝对值
    -> 每个 cell 的通量绝对值总和
    -> 除以 cell 体积
    -> 全局最大 rate
    -> 候选 deltaT
```

Stage 2 这一个阶段本身仍然不会：

- 调用 `fvc::div(phi,T)`；
- 改变 `T`；
- 调用 `runTime++`；
- 写出新的时间目录。

因此“已经算出 `deltaT`”只说明 CFL 阶段完成；完整线性对流求解器仍需 Stage 5
时间循环和 Stage 6 验证。

### Stage 2 的代码阅读入口

详细的逐行教学说明见：

```text
docs/05_stage2_cfl.md
```

推荐按照以下顺序对照源码阅读：

1. `maxCo`：从 `controlDict` 读取全局目标；
2. `mag(phi)`：把带符号面通量变成通量大小；
3. `fvc::surfaceSum(...)`：从 face 数据累加到 cell 数据；
4. `()`：打开 `tmp` 临时对象；
5. `.primitiveField()`：取得内部 cell 数组；
6. `mesh.V()`：提供每个 cell 的体积；
7. `gMax(rate)`：得到全网格最坏 cell；
8. `2.0*maxCo/rateMax`：按 CFL 公式反解时间步。

Stage 2 读取 `maxCo` 的直接头文件依赖是：

```cpp
#include "dictionary.H"
```

这里没有 `maxCo.H`。`maxCo` 只是 `system/controlDict` 中的关键字，
真正执行读取的是 `dictionary::lookupOrDefault<scalar>()`。

## 9. 历史 Stage 2 构建与最小验证记录

执行命令：

```bash
source /opt/openfoam14/etc/bashrc
cd /home/a776/workdocuments/上交船舶/slover/student_project
sh scripts/build_student_solver.sh
build/bin/explicitAdvectionFoamStudent -help
cd cases/01_sine_wave_quad
../../build/bin/explicitAdvectionFoamStudent \
    > log.explicitAdvectionFoamStudent 2>&1
```

本轮没有执行 `Allrun`，因为它会先调用 `Allclean` 删除当前可再生成的网格、
运行目录和日志。本轮复用了已经存在并通过 `checkMesh` 的网格和 `0/` 初始场，
直接运行重新编译后的学生版求解器。

验证结果：

- 学生版求解器在 OpenFOAM 14 下重新编译成功；
- `explicitAdvectionFoamStudent -help` 正常输出；
- 已有网格包含 `400` 个单元，历史 `checkMesh` 日志包含 `Mesh OK`；
- 求解器成功读入 `U`、`T`；
- 初始标量场范围为 `min(T) = -1`、`max(T) = 1`；
- 面通量范围为 `min(phi) = -0.005`、`max(phi) = 0.005`；
- `rate max = 80`；
- `maxCo = 0.2`；
- `CFL deltaT = 0.005`；
- 日志输出 `Stage 2 CFL estimate completed.`；
- 日志输出 `No time update has been implemented yet.`；
- `0/T` 与 `0.orig/T` 完全一致，说明 Stage 2 没有更新标量场；
- 没有生成 `0.005`、`0.01` 等新时间目录，说明没有推进 `runTime`。

历史 Stage 2 记录：

- 已增加 `fvcSurfaceIntegrate.H`；
- 已写入 `maxCo`、`sumPhi`、`rate`、`rateMax` 和 `deltaT`；
- 已加入 `maxCo <= 0` 和 `rateMax <= SMALL` 的错误检查；
- 已显式加入 `dictionary.H`，用于 `lookupOrDefault<scalar>()`；
- 当时已完成重新编译和 Stage 2 直接运行验证。
- 后续已经把 Stage 3 的 `fvcDiv.H`、`fvc::div(...)` 和 Stage 4 单步更新写入学生源码。

当前源码末尾的日志已经改为：

```text
Stage 4 forward Euler update completed.
No time loop has been implemented yet.
```

### Stage 3--Stage 4 当前验收计划

本轮只完成代码和教学文档，不替你编译、不运行算例。由你执行：

```bash
source /opt/openfoam14/etc/bashrc
cd /home/a776/workdocuments/上交船舶/slover/student_project
sh scripts/build_student_solver.sh
build/bin/explicitAdvectionFoamStudent -help
cd cases/01_sine_wave_quad
sh Allrun
```

检查日志：

```bash
grep -E "cells|phi min|phi max|rate max|maxCo|CFL deltaT|residual dimensions|residual min|residual max|residual integral|updated T min|updated T max|Stage 4|No time loop" \
    log.explicitAdvectionFoamStudent
```

Stage 3 需要确认：

- `residual` 的最小值和最大值是有限数；
- 周期案例的 `residual integral` 接近零；
- `residual dimensions` 是时间倒数量纲，通常类似 `[0 0 -1 0 0 0 0]`。

Stage 4 需要确认：

- `updated T min` 和 `updated T max` 是有限数，不是 `nan` 或 `inf`；
- 日志包含 `Stage 4 forward Euler update completed.`；
- 日志包含 `No time loop has been implemented yet.`；
- 不应出现 `runTime++`、`Time = 0.005` 等完整时间循环日志；
- 不应生成新的时间目录；本阶段只是内存中的单步更新；
- `0/T` 是否被改变要按当前 `IOobject::AUTO_WRITE` 和是否调用
  `runTime.write()` 判断：当前没有 `runTime.write()`，因此磁盘上的 `0/T`
  不应被覆盖。

这两个通量范围值与当前网格一致：速度分量为 `1`，对应方向面面积约为
`0.05 * 0.1 = 0.005`；相反方向的面法向会产生相反符号。

曾遇到的编译问题：

- 初稿使用了 `runTime.timeName()`；
- 当前 OpenFOAM 14 中该接口不是无参数的实例方法；
- 查阅本地 `Time.H` 后改为 `runTime.name()`；
- 修改后重新编译和算例运行均成功。

初始化过程中出现的 `opal_ifinit` 和 ParaView 配置提示来自本机环境初始化，
没有导致 solver 编译或算例运行失败。
