# 有限体积公式到 OpenFOAM 工程的完整映射

本文回答一个核心问题：

> 题目中的偏微分方程、有限体积公式、求解器 C++ 代码、`system/` 字典、`constant/` 网格文件和 `0/` 场文件，分别是什么关系？

本文对应当前学生版求解器：

```text
UDF/solver/01_advection_equation/explicitAdvectionFoamStudent/explicitAdvectionFoamStudent.C
```

案例都放在：

```text
cases/<caseName>/N<resolution>/
```

## 0. 先建立一张“地图”

如果你是第一次接触 CFD，先不要把下面的文件看成一堆互不相干的配置。
它们其实是同一个问题在不同层次的表达：

```text
题目
    -> 连续偏微分方程
    -> 有限体积离散公式
    -> OpenFOAM 网格和场文件
    -> 学生版 C++ 求解器
    -> Python 准备、运行和后处理脚本
```

每一层回答的问题不同：

| 层次 | 它回答的问题 |
|---|---|
| 连续方程 | 现实中的物理量满足什么规律？ |
| 控制体积分 | 如何把连续方程变成一个小控制体上的守恒关系？ |
| 网格 | 这些控制体在空间中具体长什么样？ |
| 场文件 | 每个 cell 初始时有什么数值？ |
| `fvSchemes` | 面上的未知量如何从 cell 值构造？ |
| 求解器 | 每一个时间步按照什么顺序计算？ |
| 后处理 | 如何判断结果是否合理？ |

可以把一次计算想成“给每个小盒子记账”：

```text
每个 cell 都有一个标量 T
每个 face 都有一个体积通量 F
通量从某个面流入或流出 cell
根据净流出量改变 cell 中的 T
```

有限体积法最重要的思想不是“在很多点上套公式”，而是：

> 对每一个控制体做守恒记账。

### 0.1 当前项目中最重要的三个名字

```text
T
    被输运的 cell-centered 标量场

U
    给定的 cell-centered 速度场

phi
    由 U 和网格面面积计算出的 face-centered 体积通量
```

初学者最容易犯的错误是把数学题里的 $\phi$ 和 OpenFOAM 代码里的
`surfaceScalarField phi` 当成同一个对象。当前项目刻意使用：

```text
题目中的标量 phi       -> C++ 中的 T
OpenFOAM 习惯中的面通量 phi -> C++ 中的 phi
```

后文每次看到 `phi` 时，都要先问自己：

```text
它是题目的被输运标量？
还是 OpenFOAM 的面体积通量？
```

### 0.2 一次时间步的全景图

求解器每一步的核心顺序是：

```text
读取 U 和 T
    |
    v
根据 U 和网格计算 phi
    |
    v
根据 phi 和 cell 体积计算 CFL 时间步 deltaT
    |
    v
根据 phi、T 和 fvSchemes 计算 residual
    |
    v
T = T - deltaT*residual
    |
    v
修正边界条件
    |
    v
按 controlDict 写出结果
```

这张图就是整个 `explicitAdvectionFoamStudent.C` 的骨架。后面的每一节，
都在解释这张图中的一个箭头。

### 0.3 “输入”和“计算结果”如何区分

```text
长期维护的输入：
    scripts/configs/01_advection_equation/*.json
    cases/.../0.orig/
    cases/.../system/
    UDF/solver/...

运行时生成的内容：
    cases/.../0/
    cases/.../constant/polyMesh/
    cases/.../时间目录/
    data/
    figures/
```

`0.orig/` 是“原始初始场”，`0/` 是“这次运行实际读入的副本”。
`constant/polyMesh/` 是“网格数据库”，不是一个可有可无的图片文件。

## 1. 连续问题

第一题求解二维守恒型线性对流方程：

$$
\frac{\partial \phi}{\partial t}
+
\nabla\cdot(\boldsymbol{u}\phi)=0.
$$

上式中的 `+` 表示“时间变化项”和“空间输运项”相加后为零，不能省略。

在本项目中，题目里的标量 $\phi$ 为了避免和 OpenFOAM 约定的面通量字段
`phi` 重名，实际写成了标量场：

```text
数学标量 phi  ->  OpenFOAM 场 T
数学速度 u    ->  OpenFOAM 场 U
```

因此代码实际表达的是：

$$
\frac{\partial T}{\partial t}
+
\nabla\cdot(\boldsymbol{U}T)=0.
$$

这不是换了一个方程，而只是把同一个数学标量 $\phi$ 改名为代码字段 `T`，
把数学速度 $\boldsymbol{u}$ 改名为代码字段 `U`。

这里的 `T` 不是温度专用变量，只是本项目给被输运标量取的字段名。

## 2. 控制体积分

对一个控制体 $\Omega_c$ 积分：

$$
\int_{\Omega_c}
\frac{\partial T}{\partial t}\,\mathrm dV
+
\int_{\Omega_c}
\nabla\cdot(\boldsymbol{U}T)\,\mathrm dV
=0.
$$

有限体积法使用高斯定理把第二项变成边界面积分：

$$
\frac{\mathrm d}{\mathrm dt}
\int_{\Omega_c}T\,\mathrm dV
+
\int_{\partial\Omega_c}
(\boldsymbol{U}T)\cdot\boldsymbol{n}\,\mathrm dS
=0.
$$

如果把单元中心值 $T_c$ 作为单元平均值，则得到：

$$
V_c\frac{\mathrm dT_c}{\mathrm dt}
+
\sum_{f\in\partial\Omega_c}
F_{cf}T_f
=0.
$$

其中：

$$
V_c=|\Omega_c|,
\qquad
F_{cf}=\boldsymbol{U}_f\cdot\boldsymbol{S}_{cf}.
$$

符号含义：

| 数学对象 | 含义 | 工程中的对应 |
|---|---|---|
| $\Omega_c$ | 第 $c$ 个控制体 | `fvMesh` 中的一个 cell |
| $V_c$ | 控制体体积 | `mesh.V()` |
| $f$ | 控制体的一个面 | `mesh.Sf()`、owner/neighbour 拓扑 |
| $\boldsymbol{S}_{cf}$ | 带方向的面面积向量 | `mesh.Sf()` |
| $T_c$ | cell 中心标量 | `volScalarField T` |
| $\boldsymbol{U}_c$ | cell 中心速度 | `volVectorField U` |
| $F_{cf}$ | 有方向体积通量 | `surfaceScalarField phi` |

## 3. 时间离散

使用前向 Euler：

$$
\frac{\mathrm dT_c}{\mathrm dt}
\approx
\frac{T_c^{n+1}-T_c^n}{\Delta t}.
$$

代入控制体方程：

$$
V_c
\frac{T_c^{n+1}-T_c^n}{\Delta t}
+
\sum_fF_{cf}T_f^n
=0.
$$

整理得到：

$$
T_c^{n+1}
=
T_c^n
-
\frac{\Delta t}{V_c}
\sum_fF_{cf}T_f^n.
$$

本项目采用体积归一化残差：

$$
R_c^n
=
\frac{1}{V_c}
\sum_fF_{cf}T_f^n.
$$

所以更新式写成：

$$
T_c^{n+1}=T_c^n-\Delta t R_c^n.
$$

对应求解器代码：

```cpp
const dimensionedScalar deltaTDim
(
    "deltaT",
    dimTime,
    stepDeltaT
);

T = T - deltaTDim*residual;
```

这里 `residual` 已经包含 $1/V_c$，不能再额外除以 `mesh.V()`。

## 4. 面通量公式对应哪里

数学公式：

$$
F_f=\boldsymbol{U}_f\cdot\boldsymbol{S}_f.
$$

求解器代码：

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

其中：

```text
U
    -> volVectorField
    -> OpenFOAM 根据网格和插值规则得到面上的 U_f
    -> 与面面积向量 Sf 点乘
    -> surfaceScalarField phi
```

`phi` 的名字遵循 OpenFOAM 中对体积通量的常见命名。它和题目中的标量
$\phi$ 不是同一个数学对象：

```text
题目标量 phi       -> 代码 T
OpenFOAM 面通量 phi -> 代码 phi
```

### 4.1 `surfaceScalarField` 为什么是面标量场

$\boldsymbol{U}_f$ 是向量，$\boldsymbol{S}_f$ 也是向量：

$$
\boldsymbol{U}_f\cdot\boldsymbol{S}_f
$$

点乘结果是一个标量，而且每个面都有一个通量值，因此类型是：

```cpp
surfaceScalarField
```

不是：

```cpp
volScalarField
```

因为 `volScalarField` 表示每个 cell 一个值，而通量定义在 face 上。

### 4.2 `NO_READ` 和 `NO_WRITE`

`phi` 不是用户在 `0/phi` 中给定的初始场，而是运行时由 `U` 和网格计算出来的：

```text
0/U + constant/polyMesh
    -> fvc::flux(U)
    -> phi
```

所以：

```cpp
IOobject::NO_READ
```

表示不从磁盘读取 `phi`。

当前求解器只把 `phi` 作为内部计算量使用，因此：

```cpp
IOobject::NO_WRITE
```

表示不自动把它写到时间目录。

### 4.3 面通量代码逐行阅读

把当前代码完整写出来，并在旁边解释：

```cpp
// surfaceScalarField：
//   每一个网格面保存一个 scalar。
//   这里的 scalar 就是该面上的有方向体积通量 F_f。
surfaceScalarField phi
(
    // IOobject 不是通量数值本身。
    // 它只是告诉 OpenFOAM 如何管理这个对象。
    IOobject
    (
        "phi",                 // 字段名字；这里遵循 OpenFOAM 通量命名习惯
        runTime.name(),        // 所属时间目录，例如 "0"
        mesh,                  // 该字段属于当前 fvMesh
        IOobject::NO_READ,     // 不从 0/phi 读取
        IOobject::NO_WRITE     // 不把 phi 自动写入结果目录
    ),

    // 这是这个场的初始数值来源。
    // fvc::flux(U) 根据 U 和网格面几何计算每个面上的通量。
    fvc::flux(U)
);
```

最容易误读的地方是最后两个部分：

```cpp
IOobject(...),
fvc::flux(U)
```

它们分别是：

```text
IOobject(...) -> 对象的“身份证和读写说明”
fvc::flux(U) -> 对象的“实际初始数值”
```

可以类比为：

```text
先创建一个表格：
    名字叫 phi，属于 mesh，不从文件读

再把计算结果填入表格：
    fvc::flux(U)
```

`fvc::flux(U)` 不是一个自己凭空知道公式的 Python 函数。它会使用：

```text
U 的 cell 值
mesh 的 cell-face 拓扑
mesh 的面面积向量 Sf
OpenFOAM 的面插值规则
```

最终实现：

$$
F_f=\boldsymbol{U}_f\cdot\boldsymbol{S}_f.
$$

## 5. 面值 $T_f$ 和空间格式

有限体积通量中需要面值：

$$
T_f.
$$

不同空间格式的核心区别，就是如何由 cell 值构造 $T_f$。

### 5.1 一阶迎风

对于面通量方向已经确定的面：

$$
T_f=
\begin{cases}
T_{\mathrm{owner}},&F_f\geq0,\\
T_{\mathrm{neighbour}},&F_f<0.
\end{cases}
$$

字典：

```foam
divSchemes
{
    default         none;
    div(phi,T)      Gauss upwind;
}
```

JSON：

```json
"schemeName": "upwind",
"divScheme": "Gauss upwind"
```

代码的作用是把 JSON 中的 `divScheme` 写入：

```text
cases/<caseName>/Nxx/system/fvSchemes
```

求解器调用：

```cpp
fvc::div(phi, T, "div(phi,T)")
```

`"div(phi,T)"` 是查找 `fvSchemes` 中离散格式的关键字。

### 5.2 `linearUpwind`

当前项目使用：

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

JSON：

```json
"schemeName": "linearUpwind",
"divScheme": "Gauss linearUpwind grad(T)",
"gradTScheme": "Gauss linear"
```

它仍然是上游偏置格式，但面值中加入了梯度重构。这里的 `grad(T)` 不是
另一个独立方程，而是为 `linearUpwind` 提供面值重构所需的梯度。

当前项目没有单独名为 `"secondOrderUpwind"` 的实现。若题目把“二阶迎风”
作为泛称，当前对应的是 `linearUpwind`；严格的格式名称仍应以
`divSchemes` 中的 OpenFOAM 语法为准。

## 6. 残差公式对应哪里

数学定义：

$$
R_c^n
=
\frac{1}{V_c}
\sum_fF_{cf}T_f^n.
$$

求解器代码：

```cpp
tmp<volScalarField> tResidual
(
    fvc::div(phi, T, "div(phi,T)")
);

const volScalarField& residual = tResidual();
```

这里：

```text
phi
    -> 面通量 F_f
T
    -> cell 标量
fvSchemes
    -> 面值 T_f 的计算规则
fvc::div(...)
    -> 组装面通量并除以 V_c
residual
    -> 每个 cell 一个 R_c
```

`fvc::div` 返回的是 `volScalarField`，因为残差最终定义在 cell 上。

如果把未除以体积的量单独记为：

$$
\mathcal Q_c^n=\sum_fF_{cf}T_f^n,
$$

则：

$$
R_c^n=\frac{\mathcal Q_c^n}{V_c}.
$$

当前项目采用前一种更适配 OpenFOAM 的表示：直接使用 `fvc::div` 得到
已经除以体积的 `residual`。

### 6.1 残差代码逐行阅读

```cpp
// tmp<volScalarField>：
//   fvc::div 返回一个临时的 volScalarField。
//   tmp 是 OpenFOAM 用来管理临时对象生命周期的包装器。
//
//   这里不要把 tmp 理解成“另一个数学量”。
//   数学量仍然是 R_c；
//   tmp 只是 C++ 内存管理方式。
tmp<volScalarField> tResidual
(
    // 第一个参数 phi：
    //   面体积通量 F_f。
    //
    // 第二个参数 T：
    //   被输运标量的 cell 场。
    //
    // 第三个参数 "div(phi,T)":
    //   离散算子的名字，必须和 fvSchemes 中的条目一致。
    fvc::div(phi, T, "div(phi,T)")
);
```

把函数括号拆开：

```cpp
fvc::div
(
    phi,              // 面通量
    T,                // cell 标量
    "div(phi,T)"      // fvSchemes 查找名称
)
```

它的数学过程是：

```text
phi
    -> 得到每个面上的 F_f

T + fvSchemes
    -> 得到每个面上的 T_f

F_f*T_f
    -> 对每个 cell 周围的面求和

除以 V_c
    -> 得到 R_c
```

然后：

```cpp
// tResidual()：
//   从 tmp 包装器中取得真正的 volScalarField。
//
// const volScalarField&：
//   用只读引用访问它，不复制整张残差场。
const volScalarField& residual = tResidual();
```

这里的两个括号要区分：

```text
tResidual
    -> 一个 tmp 对象

tResidual()
    -> 调用 tmp 的访问操作，取得内部场
```

残差积分检查：

```cpp
const scalar residualIntegral
(
    // residual.primitiveField()
    //     -> [R_0, R_1, R_2, ...]
    //
    // mesh.V().primitiveField()
    //     -> [V_0, V_1, V_2, ...]
    //
    // * 是逐 cell 相乘：
    //     [R_0*V_0, R_1*V_1, ...]
    //
    // gSum 是全网格求和。
    gSum
    (
        residual.primitiveField()
        *
        mesh.V().primitiveField()
    )
);
```

数学对应：

$$
\texttt{residualIntegral}
=
\sum_cR_cV_c
=
\sum_c\sum_fF_{cf}T_f.
$$

周期边界和守恒通量正确时，这个值应该接近零。

## 7. CFL 公式对应哪里

当前求解器使用：

$$
\mathrm{Co}_{\max}
=
\frac{1}{2}
\Delta t
\max_c
\left(
\frac{\sum_f|F_{cf}|}{V_c}
\right).
$$

令：

$$
\lambda_c
=
\frac{\sum_f|F_{cf}|}{V_c},
$$

则：

$$
\Delta t
=
\frac{2\,\mathrm{Co}_{\mathrm{target}}}
{\max_c\lambda_c}.
$$

求解器代码：

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

逐项对应：

| 数学量 | 代码 |
|---|---|
| $\mathrm{Co}_{\mathrm{target}}$ | `maxCo` |
| $|F_f|$ | `mag(phi)` |
| $\sum_f|F_{cf}|$ | `fvc::surfaceSum(mag(phi))` |
| $V_c$ | `mesh.V()` |
| $\max_c$ | `gMax(rate)` |
| $\Delta t$ | `deltaT` |

JSON：

```json
"maxCo": 0.2
```

脚本把它写到：

```text
cases/<caseName>/Nxx/system/controlDict
```

求解器再从 `controlDict` 读取它。

### 7.1 CFL 代码逐行阅读

```cpp
// const：
//   这个变量在创建后不允许再次赋值。
//
// scalar：
//   OpenFOAM 的浮点标量类型。
//
// maxCo：
//   目标 Courant 数，不是某个 cell 的局部结果，
//   而是整个计算案例希望控制的上限。
const scalar maxCo
(
    // lookupOrDefault<scalar>：
    //   从 dictionary 中读取一个 scalar。
    controlDict.lookupOrDefault<scalar>
    (
        "maxCo",       // controlDict 中要查找的关键字
        0.2            // 没有写 maxCo 时的默认值
    )
);
```

这里的括号不是多余的格式，而是 C++ 的构造和函数调用：

```text
const scalar maxCo( 某个数值 );
    -> 创建一个名叫 maxCo 的 scalar

lookupOrDefault<scalar>( "maxCo", 0.2 );
    -> 查字典，没有就返回 0.2
```

接下来：

```cpp
// mag(phi)：
//   把有方向的 F_f 变成大小 |F_f|。
//   CFL 需要知道一个 cell 周围“总共能通过多少通量”，
//   所以这里不能保留正负号。
//
// fvc::surfaceSum(...)：
//   输入是面场，输出是每个 cell 的面值总和。
//
// ()：
//   fvc::surfaceSum 返回的是 tmp 临时对象，
//   这里的 () 用来访问其中的场对象。
//
// primitiveField()：
//   取出内部 cell 数值数组，便于构造 scalarField。
scalarField sumPhi
(
    fvc::surfaceSum(mag(phi))().primitiveField()
);
```

数学对应：

$$
\texttt{sumPhi}[c]
=
\sum_{f\in\partial\Omega_c}|F_{cf}|.
$$

注意 `sumPhi` 已经从“每个面一个值”变成了“每个 cell 一个值”。

然后：

```cpp
// mesh.V() 返回每个 cell 的体积场。
// primitiveField() 取出：
//
//     [V_0, V_1, V_2, ...]
//
// 两个 scalarField 做逐元素除法：
//
//     rate[c] = sumPhi[c] / V[c]
scalarField rate
(
    sumPhi/mesh.V().primitiveField()
);
```

数学对应：

$$
\texttt{rate}[c]
=
\frac{\sum_f|F_{cf}|}{V_c}.
$$

最后：

```cpp
// gMax(rate)：
//   在所有 cell 中寻找最大的 rate。
//   g 表示 global；并行计算时还会跨进程归约。
const scalar rateMax = gMax(rate);

// 根据：
//
//     Co_target = 0.5*deltaT*rateMax
//
// 反解本步候选时间步。
const scalar deltaT = 2.0*maxCo/rateMax;
```

整个 CFL 计算链可以写成：

```text
phi
    -> mag(phi)
    -> surfaceSum(...)
    -> sumPhi[c]
    -> sumPhi[c]/V[c]
    -> rate[c]
    -> gMax(rate)
    -> rateMax
    -> 2*maxCo/rateMax
    -> deltaT
```

## 8. 时间循环和最终时间

终止时间来自：

```foam
endTime 1;
```

或第二案例中的：

```foam
endTime 6.2831853071795862;
```

求解器循环：

```cpp
while (endTime - runTime.value() > timeTolerance)
{
    const scalar remainingTime = endTime - runTime.value();
    const scalar stepDeltaT = min(deltaT, remainingTime);
    runTime.setDeltaT(stepDeltaT);
    runTime++;
    ...
    T = T - deltaTDim*residual;
    T.correctBoundaryConditions();
    runTime.write();
}
```

最后一次时间步使用：

$$
\Delta t^n
=
\min(\Delta t_{\mathrm{CFL}},
t_{\mathrm{end}}-t^n),
$$

因此不会超过终止时间。

循环结束后还会调用：

```cpp
runTime.writeNow();
```

因为普通的 `runTime.write()` 遵守 `writeInterval`，不一定写出最后一个时间层。

### 8.1 时间循环代码逐行阅读

```cpp
// 从 OpenFOAM 的 Time 对象中读取终止时间。
// .value() 把带量纲的时间对象转换成可参与计算的 scalar。
const scalar endTime = runTime.endTime().value();

// 给时间比较设置容差，避免浮点数误差导致：
//
//     6.283185307179586 != 6.2831853071795862
//
// 从而多进入一次或少进入一次循环。
const scalar timeTolerance
(
    1.0e-12*max(1.0, mag(endTime))
);

// label 是 OpenFOAM 的整数类型，用来记录时间步编号。
label step = 0;
```

循环条件：

```cpp
while (endTime - runTime.value() > timeTolerance)
```

括号里的表达式是：

```text
endTime - runTime.value()
    -> 还剩多少物理时间

> timeTolerance
    -> 剩余时间是否还大到值得推进一步
```

进入循环后：

```cpp
// 当前物理时间是 runTime.value()。
// endTime - 当前时间就是剩余时间。
const scalar remainingTime
(
    endTime - runTime.value()
);

// 正常使用 CFL 时间步；
// 最后一步如果不足一个完整 CFL 步，就只走剩余时间。
const scalar stepDeltaT
(
    min(deltaT, remainingTime)
);
```

这对应：

$$
\Delta t^n
=
\min\left(
\Delta t_{\mathrm{CFL}},
t_{\mathrm{end}}-t^n
\right).
$$

然后：

```cpp
// 把本次真正采用的步长写入 Time 对象。
// 这一步只设置 deltaT，还没有推进当前时间。
runTime.setDeltaT(stepDeltaT);

// ++ 操作推进 OpenFOAM 的时间标签：
//
//     t^n -> t^(n+1)
//
// 此时 T 仍然是更新前的 T^n。
runTime++;

// 时间步编号加一。
++step;
```

为什么先 `runTime++` 再计算残差？

```text
runTime++       -> 告诉 OpenFOAM 当前正在写 t^(n+1)
T                -> 仍保存旧场 T^n
fvc::div(phi,T)  -> 用旧场计算 R^n
T = T-deltaT*R   -> 得到新场 T^(n+1)
runTime.write()  -> 把新场写入 t^(n+1)
```

这是“时间标签”和“场变量更新”在代码中的分工。

### 8.2 每一步中残差和更新的详细顺序

```cpp
// 当前 T 仍然是 T^n。
tmp<volScalarField> tResidual
(
    fvc::div(phi, T, "div(phi,T)")
);

// 取出只读残差场 R^n。
const volScalarField& residual = tResidual();

// 给普通数值 stepDeltaT 加上时间量纲。
const dimensionedScalar deltaTDim
(
    "deltaT",       // 这个量在 OpenFOAM 中的名字
    dimTime,        // 时间量纲
    stepDeltaT      // 实际数值
);

// 显式前向 Euler 更新：
//
//     T^(n+1) = T^n - deltaT*R^n
T = T - deltaTDim*residual;

// 让 T 的边界 patch 满足当前网格边界类型。
T.correctBoundaryConditions();
```

这里 `dimensionedScalar` 很重要：

```text
stepDeltaT
    -> 普通 C++ 数值，没有 OpenFOAM 物理量纲

deltaTDim
    -> 带有 dimTime 的时间量

residual
    -> 通常具有 1/time 的量纲

deltaTDim*residual
    -> 无量纲，和 T 可以相减
```

### 8.3 写出结果

```cpp
// T 在创建时使用 AUTO_WRITE。
// runTime.write() 会根据 controlDict 的 writeControl 和 writeInterval
// 决定当前时间层是否真正写入磁盘。
runTime.write();
```

最后：

```cpp
// 不管最后一步是否刚好满足 writeInterval，
// writeNow() 都强制把当前 T 写入当前时间目录。
const bool finalWriteOK = runTime.writeNow();
```

这里返回 `bool`：

```text
true  -> 写出成功
false -> 写出失败
```

所以代码继续检查：

```cpp
if (!finalWriteOK)
{
    FatalErrorInFunction
        << "Failed to write the final field at time "
        << runTime.name()
        << exit(FatalError);
}
```

如果最终写出失败，程序不会安静地结束，而是明确报错。

## 9. `system/` 中每个文件的职责

### 9.1 `system/controlDict`

它描述运行控制参数：

```foam
application     explicitAdvectionFoamStudent;
startFrom       startTime;
startTime       0;
stopAt          endTime;
endTime         1;
deltaT          0.01;
writeControl    timeStep;
writeInterval   10;
maxCo           0.2;
velocityField   U;
advectedField   T;
```

映射：

| `controlDict` 项 | 数学/程序含义 |
|---|---|
| `application` | 使用哪个可执行求解器 |
| `startTime` | 初始时间 $t^0$ |
| `endTime` | 终止时间 |
| `deltaT` | 初始占位时间步，实际步长由 CFL 重算 |
| `writeControl` | 结果写出控制方式 |
| `writeInterval` | 写出频率 |
| `maxCo` | CFL 目标值 |
| `velocityField` | 速度字段名 |
| `advectedField` | 被输运字段名 |

### 9.2 `system/fvSchemes`

它描述空间和时间离散格式：

```foam
ddtSchemes
{
    default Euler;
}

gradSchemes
{
    default Gauss linear;
}

divSchemes
{
    default none;
    div(phi,T) Gauss upwind;
}
```

本项目求解器不调用 `fvm::div` 或线性方程组求解器，因此
`system/fvSolution` 不承担矩阵求解任务。

### 9.3 `system/fvSolution`

当前案例是完全显式更新：

```cpp
T = T - deltaTDim*residual;
```

不组装线性方程组，因此 `fvSolution` 只保留标准空结构，方便 OpenFOAM
案例目录完整。

### 9.4 `system/blockMeshDict`

它只用于四边形案例。核心内容：

```foam
hex (0 1 2 3 4 5 6 7) (N N 1) simpleGrading (1 1 1)
```

含义：

```text
x 方向 N 个 cell
y 方向 N 个 cell
z 方向 1 个 cell
```

所以四边形案例的总 cell 数是：

$$
N_{\mathrm{cell}}=N\times N\times1.
$$

### 9.5 `system/createPatchDict`

它只在三角形案例运行时发挥作用。Gmsh 先生成：

```text
xMinSource, xMaxSource,
yMinSource, yMaxSource,
zMinSource, zMaxSource
```

然后 `createPatch` 把它们改造成 OpenFOAM 最终 patch：

```text
xMin <-> xMax  cyclic
yMin <-> yMax  cyclic
zMin, zMax      empty
```

## 10. `0.orig/`、`0/` 和场文件

### 10.1 `0.orig/U`

它是速度场输入：

```foam
class       volVectorField;
object      U;
dimensions  [0 1 -1 0 0 0 0];
```

第一案例：

```foam
internalField uniform (1 1 0);
```

对应：

$$
\boldsymbol{u}=(1,1).
$$

第二案例：

```foam
internalField nonuniform List<vector>
...
```

每个 cell 的速度由：

$$
\boldsymbol{u}(x,y)
=(0.5-y,\;x-0.5)
$$

计算得到。

### 10.2 `0.orig/T`

它是初始标量场输入：

第一案例：

$$
T(x,y,0)=\sin(2\pi(x+y)).
$$

四边形网格中，脚本按规则 cell centre 生成：

```text
x_i = (i+0.5)/N
y_j = (j+0.5)/N
```

三角形网格中，脚本先生成网格，再读取真实 cell centre，最后生成：

```text
T(xc, yc, 0)
```

第二案例中，`T` 是三个函数之和：

$$
T_0(x,y)=T_D(x,y)+T_C(x,y)+T_H(x,y).
$$

### 10.3 为什么要有 `0.orig/` 和 `0/`

```text
0.orig/
    长期维护的初始场源文件

0/
    本次运行实际读入的初始场副本
```

`Allrun` 会清理并重新复制：

```text
0.orig -> 0
```

这样运行产生的修改不会污染原始输入。

## 11. `constant/` 中到底有什么

### 11.1 `constant/polyMesh`

这是 OpenFOAM 的真实网格数据库，包含：

```text
points
faces
owner
neighbour
boundary
```

它们共同定义：

```text
cell 几何
face 几何
cell-face 邻接关系
boundary patch
```

求解器通过：

```cpp
#include "createMesh.H"
```

得到：

```cpp
fvMesh mesh;
```

后续 `mesh.V()`、`mesh.Sf()`、`mesh.nCells()` 和 `fvc::flux(U)` 都依赖它。

### 11.2 `constant/C` 和 `constant/Vc`

三角形案例运行：

```bash
foamPostProcess -constant -func writeCellCentres
foamPostProcess -constant -func writeCellVolumes
```

后得到：

```text
constant/C
constant/Vc
```

它们不是求解器必须读取的物理输入，而是三角形案例的脚本和后处理辅助数据：

| 文件 | 含义 |
|---|---|
| `constant/C` | 每个 cell 的真实中心坐标 |
| `constant/Vc` | 每个 cell 的真实体积 |

求解器本身仍直接从 `fvMesh` 读取几何；Python 后处理用 `C` 和 `Vc` 计算
三角形案例的初值、精确解和体积加权误差。

### 11.3 为什么通常没有 `transportProperties`

当前题目是给定速度的线性对流：

$$
\boldsymbol{U}
\quad\text{是输入场，而不是待求解的动量方程未知量。}
$$

因此没有密度、黏度、扩散系数等物理模型，不需要 `transportProperties`
来完成本题。

## 12. 网格在公式中扮演什么角色

网格并不是“只负责画图”。它直接决定：

```text
V_c       -> 控制体体积
S_f       -> 面面积向量
owner     -> 面的拥有单元
neighbour -> 内部面的另一侧单元
boundary  -> 边界面和边界类型
```

因此网格参与：

1. 面通量 $F_f=\boldsymbol{U}_f\cdot\boldsymbol{S}_f$；
2. 残差 $R_c=(1/V_c)\sum_fF_{cf}T_f$；
3. CFL 时间步；
4. 边界条件；
5. 误差的体积加权。

## 13. 四边形与三角形的区别

### 四边形案例

```text
JSON meshType=quad
    -> 复制四边形模板
    -> 修改 blockMeshDict 的 (N N 1)
    -> blockMesh
    -> constant/polyMesh
```

每个二维控制体是一个四边形柱体，z 方向只有一层。

### 三角形案例

```text
JSON meshType=tri
    -> Gmsh 生成二维三角形
    -> 沿 z 方向挤出一层
    -> mesh.msh
    -> gmshToFoam
    -> createPatch
    -> constant/polyMesh
```

当前 Gmsh 脚本通过 `setTransfiniteCurve` 和 `setTransfiniteSurface`
让 `N` 表示每条边的区间数；但二维面被三角形剖分，所以：

```text
N=20 不代表 20*20 个三角形
```

实际三角形数量由剖分结果决定，必须以 `checkMesh` 和导入后的网格为准。

## 14. 周期边界对应哪里

周期边界的数学含义是：

$$
T(0,y,t)=T(1,y,t),
\qquad
T(x,0,t)=T(x,1,t).
$$

四边形案例在 `system/blockMeshDict` 中直接写：

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
yMin
{
    type cyclic;
    neighbourPatch yMax;
}
yMax
{
    type cyclic;
    neighbourPatch yMin;
}
```

场文件 `0.orig/U` 和 `0.orig/T` 的 `boundaryField` 也要匹配：

```foam
xMin { type cyclic; }
xMax { type cyclic; }
yMin { type cyclic; }
yMax { type cyclic; }
zMin { type empty; }
zMax { type empty; }
```

三角形案例中：

```text
Gmsh physical group
    -> createPatchDict
    -> xMin/xMax/yMin/yMax cyclic
```

求解器每一步调用：

```cpp
T.correctBoundaryConditions();
```

这使 `T` 的边界 patch 与周期拓扑保持一致。`fvc::flux` 和 `fvc::div`
也会使用网格和字段的边界信息处理跨周期面。

## 15. 质量守恒和残差积分

周期边界下，所有内部面和周期面上的有向通量应相互抵消，因此：

$$
\sum_cR_cV_c
=
\sum_c\sum_fF_{cf}T_f
\approx0.
$$

代码：

```cpp
const scalar residualIntegral
(
    gSum(residual.primitiveField()*mesh.V().primitiveField())
);
```

注意这个量是“残差的体积分”，不是误差，也不是把所有 `T` 相加。

质量守恒检查使用：

$$
M^n=\sum_cT_c^nV_c.
$$

后处理中的：

```text
initialMass
finalMass
massChange
normalizedMassError
```

用于检查输运过程中总量是否被保持。

## 16. 从公式到程序的总链条

```text
连续方程
    -> 控制体积分
    -> F_cf = U_f dot S_cf
    -> phi = fvc::flux(U)
    -> 面值 T_f 由 fvSchemes 决定
    -> residual = fvc::div(phi,T)
    -> residual 已除以 V_c
    -> T = T - deltaT*residual
    -> maxCo 控制 deltaT
    -> correctBoundaryConditions()
    -> runTime.write()/writeNow()
```

文件层次对应：

```text
JSON
    -> scripts/common/case_config.py
    -> scripts/common/foam_case.py
    -> system/controlDict
    -> system/fvSchemes
    -> system/blockMeshDict 或 Gmsh
    -> constant/polyMesh
    -> 0.orig/U、0.orig/T
    -> UDF/solver/01_advection_equation/explicitAdvectionFoamStudent.C
    -> data/ 和 figures/
```

## 17. 从 `main()` 开始读源码

如果你打开 C++ 源码，建议不要从第一行开始逐字背诵，而是先找到：

```cpp
int main(int argc, char *argv[])
```

它是程序真正的入口。

### 17.1 `argc` 和 `argv`

```cpp
int main(int argc, char *argv[])
```

括号内：

```text
argc -> 命令行参数的数量
argv -> 命令行参数的字符串数组
```

例如：

```bash
explicitAdvectionFoamStudent -case cases/01_advection_equation/01_sine_wave_quad/N20
```

程序收到的参数大致是：

```text
argv[0] = explicitAdvectionFoamStudent
argv[1] = -case
argv[2] = cases/01_advection_equation/01_sine_wave_quad/N20
```

OpenFOAM 的 `argList.H`、`setRootCase.H` 会把这些普通字符串解析成
OpenFOAM 能理解的 case 路径和选项。

### 17.2 为什么源码中大量使用 `#include`

你会看到：

```cpp
#include "setRootCase.H"
#include "createTime.H"
#include "createMesh.H"
```

对于普通 C++ 初学者，容易误以为它们只是“读取三个头文件”。在 OpenFOAM
求解器中，这些文件很多是带有执行代码的 include 片段：

```text
setRootCase.H
    -> 解析 case 和命令行

createTime.H
    -> 创建 runTime

createMesh.H
    -> 创建 mesh
```

可以把它们理解成 OpenFOAM 提供的“标准启动步骤”。

### 17.3 `using namespace Foam`

```cpp
using namespace Foam;
```

OpenFOAM 的类型通常位于 `Foam` 命名空间中。没有这句时，需要写：

```cpp
Foam::scalar
Foam::volScalarField
Foam::Info
```

有了这句后，可以直接写：

```cpp
scalar
volScalarField
Info
```

它不改变数学意义，只是省略命名空间前缀。

## 18. OpenFOAM 类型和数学对象对照表

### 18.1 基础类型

```cpp
scalar
```

表示一个浮点数，例如：

```cpp
const scalar maxCo = 0.2;
```

对应数学中的一个实数。

```cpp
label
```

表示整数，例如：

```cpp
label step = 0;
```

用于时间步编号、cell 编号等离散索引。

```cpp
word
```

表示 OpenFOAM 中的一个合法名称，例如：

```cpp
word velocityName("U");
```

对应字段名字，而不是字段数值。

### 18.2 场类型

```cpp
volScalarField
```

含义：

```text
vol -> 定义在 cell/控制体上
Scalar -> 每个 cell 一个标量
Field -> 一整张场
```

当前：

```cpp
volScalarField T;
```

对应：

$$
\{T_0,T_1,T_2,\ldots,T_{N_{\mathrm{cell}}-1}\}.
$$

```cpp
volVectorField
```

对应每个 cell 一个三维向量：

```text
U_0=(u_0,v_0,w_0)
U_1=(u_1,v_1,w_1)
...
```

```cpp
surfaceScalarField
```

对应每个 face 一个标量：

```text
phi_0, phi_1, phi_2, ...
```

当前这些标量就是有方向体积通量。

### 18.3 `scalarField` 和 `volScalarField` 的区别

```cpp
scalarField rate;
```

更像一个单纯的数值数组：

```text
[rate_0, rate_1, rate_2, ...]
```

它没有完整的 OpenFOAM 字段对象信息，不负责：

```text
边界 patch
读写文件
网格注册
物理量纲
```

而：

```cpp
volScalarField T;
```

不仅保存数值，还知道：

```text
属于哪个 mesh
有哪些边界 patch
如何读写
具有什么量纲
```

## 19. `system/`、`constant/`、`0/` 的学习顺序

对于一个具体案例，例如：

```text
cases/01_advection_equation/01_sine_wave_quad/N20/
```

建议按下面的顺序阅读。

### 第一步：看 `system/controlDict`

先回答：

```text
使用哪个 solver？
从什么时候开始？
什么时候结束？
maxCo 是多少？
多久写一次结果？
字段叫 U 还是别的名字？
字段叫 T 还是别的名字？
```

重点：

```foam
application     explicitAdvectionFoamStudent;
endTime         1;
maxCo           0.2;
velocityField   U;
advectedField   T;
```

### 第二步：看 `system/fvSchemes`

回答：

```text
时间格式是什么？
梯度格式是什么？
散度项 div(phi,T) 使用什么格式？
```

对于一阶迎风：

```foam
div(phi,T)      Gauss upwind;
```

对于 `linearUpwind`：

```foam
grad(T)         Gauss linear;
div(phi,T)      Gauss linearUpwind grad(T);
```

### 第三步：看 `system/blockMeshDict` 或 `mesh/`

四边形：

```text
system/blockMeshDict
    -> vertices
    -> blocks
    -> boundary
```

三角形：

```text
mesh/mesh.msh
mesh/mesh_geometry.json
constant/polyMesh/
```

### 第四步：看 `0.orig/U` 和 `0.orig/T`

回答：

```text
速度场是什么？
初始标量场是什么？
边界类型是什么？
```

### 第五步：看 C++ 求解器

按下面顺序定位：

```text
createMesh
    -> U/T
    -> fvc::flux(U)
    -> CFL
    -> fvc::div(phi,T)
    -> T = T-deltaT*residual
    -> correctBoundaryConditions
    -> write
```

### 第六步：看日志和后处理

日志重点：

```text
cells
phi min/max
rate max
maxCo
residual integral
T min/max
final time
```

后处理重点：

```text
summary.json
time_history.csv
field_data.csv
convergence_summary.csv
```

## 20. 一个 cell 的手算示意

假设某个 cell 有四个面，体积为：

$$
V_c=0.01.
$$

四个有方向通量和面值为：

| 面 | $F_{cf}$ | $T_f$ | $F_{cf}T_f$ |
|---|---:|---:|---:|
| 左 | $-0.10$ | $0.8$ | $-0.08$ |
| 右 | $0.10$ | $0.3$ | $0.03$ |
| 下 | $-0.05$ | $0.6$ | $-0.03$ |
| 上 | $0.05$ | $0.4$ | $0.02$ |

未归一化通量和：

$$
\mathcal Q_c
=-0.08+0.03-0.03+0.02
=-0.06.
$$

体积归一化残差：

$$
R_c=\frac{\mathcal Q_c}{V_c}
=\frac{-0.06}{0.01}
=-6.
$$

若：

$$
\Delta t=0.001,
\qquad
T_c^n=0.5,
$$

则：

$$
T_c^{n+1}
=0.5-0.001(-6)
=0.506.
$$

代码中：

```cpp
T = T - deltaTDim*residual;
```

就是把每个 cell 都做同样的计算，只不过 `residual` 是整张场，不是一个
手算中的单个数。

## 21. 初学者最容易混淆的六件事

### 21.1 `T` 不是温度

当前 `T` 只是被输运标量的名字。它可以代表：

```text
浓度
染料浓度
温度
示踪物
数学题中的 phi
```

### 21.2 `phi` 有两个语境

```text
数学题 phi       -> 被输运标量
OpenFOAM phi     -> 面体积通量
```

本项目用 `T` 避免了二者冲突。

### 21.3 `fvc` 不是 `fvm`

```text
fvc -> finite-volume calculus，显式计算场
fvm -> finite-volume method，组装隐式矩阵
```

本题要求显式，所以使用：

```cpp
fvc::flux(...)
fvc::surfaceSum(...)
fvc::div(...)
```

而不是：

```cpp
fvm::div(...)
fvm::ddt(...)
```

### 21.4 `residual` 不是“误差”

`residual` 是方程右端的离散输运项：

$$
R_c=\frac{1}{V_c}\sum_fF_{cf}T_f.
$$

误差是数值解与解析解的差：

$$
e_c=T_c-T_c^{\mathrm{exact}}.
$$

二者用途完全不同。

### 21.5 `checkMesh` 不是精度验证

```text
checkMesh
    -> 检查网格几何和拓扑是否可用

L1 误差/收敛阶
    -> 检查数值方法是否准确
```

网格 `Mesh OK` 不等于数值结果一定准确。

### 21.6 `constant/C` 不是所有案例都必须有

四边形案例的规则中心可以由 $N$ 直接计算，因此后处理可以使用规则坐标。

三角形案例的 cell centre 不规则，所以额外写出：

```text
constant/C
constant/Vc
```

它们主要服务于三角形初值生成和误差后处理。

## 22. 教授式源码阅读练习

建议你每次只回答一个问题：

### 练习 1

```cpp
fvc::flux(U)
```

问自己：

```text
输入 U 的位置是什么？
U 是 cell 场还是 face 场？
输出 phi 的位置是什么？
为什么输出是 scalar 而不是 vector？
```

答案：

```text
U 在 cell 上，是 volVectorField；
OpenFOAM 将它评价到 face；
与 Sf 点乘；
输出每个 face 一个 scalar；
所以是 surfaceScalarField。
```

### 练习 2

```cpp
fvc::surfaceSum(mag(phi))().primitiveField()
```

按顺序解释：

```text
mag(phi)
    -> 面通量绝对值

surfaceSum(...)
    -> 面数据累加到 cell

()
    -> 取得 tmp 内部对象

primitiveField()
    -> 取得 cell 数值数组
```

### 练习 3

```cpp
fvc::div(phi, T, "div(phi,T)")
```

问自己：

```text
phi 提供什么？
T 提供什么？
"div(phi,T)" 去哪里找？
为什么返回 volScalarField？
```

答案：

```text
phi 提供面通量 F_f；
T 提供 cell 标量；
"div(phi,T)" 在 system/fvSchemes 中查找；
散度结果属于每个 cell，所以是 volScalarField。
```

## 23. 最终记忆版

不要试图先背 OpenFOAM API。先记住下面这组数学关系：

$$
F_f=\boldsymbol{U}_f\cdot\boldsymbol{S}_f,
$$

$$
R_c=\frac{1}{V_c}\sum_fF_{cf}T_f,
$$

$$
T_c^{n+1}=T_c^n-\Delta tR_c^n.
$$

然后再记住三个代码接口：

```cpp
fvc::flux(U)
    -> F_f

fvc::div(phi, T, "div(phi,T)")
    -> R_c

T = T - deltaTDim*residual
    -> T^(n+1)
```

最后记住三个文件层：

```text
system/
    -> 怎么离散、什么时候运行、什么时候写出

constant/
    -> 网格几何和拓扑

0/
    -> 初始时刻的 U 和 T
```

这样你就不是在背“某个函数怎么写”，而是在理解：

```text
数学公式
    -> 数据结构
    -> OpenFOAM API
    -> 具体文件
    -> 数值结果
```

## 24. 继续追到 OpenFOAM 源码内部

前面的内容已经说明了学生版求解器“调用了什么”。这一节再向下走一层，
回答一个更严格的问题：

> `fvc::div(phi, T, "div(phi,T)")` 到底怎样从一个 C++ 函数调用，
> 变成了 $\frac{1}{V_c}\sum_f F_{cf}T_f$？

你现在不需要一次读懂 OpenFOAM 的全部模板和继承体系。先只跟踪下面
这条最短调用链：

```text
学生版求解器
    |
    | fvc::div(phi, T, "div(phi,T)")
    v
finiteVolume/finiteVolume/fvc/fvcDiv.C
    |
    | fv::convectionScheme<Type>::New(...)
    v
具体对流格式，例如 Gauss upwind
    |
    | fvcDiv(flux, T)
    v
面通量组装
    |
    | fvc::surfaceIntegrate(...)
    v
每个 cell 的面通量求和，再除以 V_c
    |
    v
volScalarField residual
```

### 24.1 学生版代码所在的位置

当前学生版代码中的主调用位于：

```text
UDF/solver/01_advection_equation/explicitAdvectionFoamStudent/explicitAdvectionFoamStudent.C
```

关键位置可以用下面的命令查找：

```bash
rg -n \
  "fvc::flux|surfaceSum|fvc::div|correctBoundaryConditions|writeNow" \
  UDF/solver/01_advection_equation/explicitAdvectionFoamStudent/explicitAdvectionFoamStudent.C
```

输出中的含义是：

```text
fvc::flux
    -> 由速度场 U 生成面体积通量 phi

surfaceSum
    -> 计算 CFL 所需的 sum_f(|F_cf|)

fvc::div
    -> 计算离散散度 residual

correctBoundaryConditions
    -> 更新 T 的边界 patch 数值

writeNow
    -> 强制写出最后时间层
```

不要把 `rg -n` 当成数值计算命令。它只是搜索文件并显示行号，方便你
从源码跳转到真正的实现。

### 24.2 `fvcDiv.C` 中最关键的函数

在 OpenFOAM 14 中，相关源码位于：

```text
/opt/openfoam14/src/finiteVolume/finiteVolume/fvc/fvcDiv.C
```

与本项目匹配的重载可以概括为下面的形式：

```cpp
template<class Type>
tmp<VolField<Type>>
div
(
    const surfaceScalarField& flux,
    const VolField<Type>& vf,
    const word& name
)
{
    return fv::convectionScheme<Type>::New
    (
        vf.mesh(),                         // 当前网格
        flux,                              // 面体积通量 phi
        vf.mesh().schemes().div(name)      // 从 fvSchemes 查找格式
    ).ref().fvcDiv(flux, vf);              // 调用具体格式的显式散度
}
```

把它和学生版调用逐个对照：

```cpp
fvc::div
(
    phi,              // 对应源码参数 flux
    T,                // 对应源码参数 vf
    "div(phi,T)"      // 对应源码参数 name
)
```

对应关系：

| 学生版实参 | 源码形参 | 数学含义 |
|---|---|---|
| `phi` | `flux` | 面上的 $F_f$ |
| `T` | `vf` | cell 上的 $T_c$ |
| `"div(phi,T)"` | `name` | 离散格式查找名 |
| 返回值 | `tmp<VolField<Type>>` | 每个 cell 一个散度值 |

这里的 `template<class Type>` 表示：这套函数可以处理不同类型的场。
本题的 `T` 是标量，所以最终的 `Type` 是 `scalar`，返回的就是
`volScalarField`。

### 24.3 `"div(phi,T)"` 为什么不是数学变量

下面三处名字看起来相同，但身份不同：

```text
phi
    -> C++ 中的 surfaceScalarField 变量

T
    -> C++ 中的 volScalarField 变量

"div(phi,T)"
    -> 一个 word 字符串，用来索引 fvSchemes
```

`"div(phi,T)"` 不会自己计算任何数值。它的作用类似字典中的键：

```foam
divSchemes
{
    div(phi,T)      Gauss upwind;
}
```

当源码执行：

```cpp
vf.mesh().schemes().div("div(phi,T)")
```

OpenFOAM 才会从 `system/fvSchemes` 取出：

```text
Gauss upwind
```

然后根据这个配置创建对应的 `convectionScheme` 对象。

所以完整逻辑是：

```text
"div(phi,T)"
    -> 查到 "Gauss upwind"
    -> 创建 Gauss 对流格式对象
    -> 对 T 做面值重构
    -> 生成 F_f*T_f
    -> 求散度
```

### 24.4 `fvcDiv` 为什么最后会除以体积

继续查看：

```text
/opt/openfoam14/src/finiteVolume/finiteVolume/fvc/fvcSurfaceIntegrate.C
```

其中的核心操作可以抽象成：

```cpp
template<class Type>
void surfaceIntegrate
(
    Field<Type>& ivf,                  // 要写入的 cell 内部数组
    const SurfaceField<Type>& ssf      // 输入的面通量场
)
{
    const fvMesh& mesh = ssf.mesh()();

    const labelUList& owner = mesh.owner();
    const labelUList& neighbour = mesh.neighbour();

    forAll(owner, facei)
    {
        ivf[owner[facei]] += ssf[facei];
        ivf[neighbour[facei]] -= ssf[facei];
    }

    forAll(mesh.boundary(), patchi)
    {
        const labelUList& pFaceCells =
            mesh.boundary()[patchi].faceCells();

        const fvsPatchField<Type>& pssf =
            ssf.boundaryField()[patchi];

        forAll(mesh.boundary()[patchi], facei)
        {
            ivf[pFaceCells[facei]] += pssf[facei];
        }
    }

    ivf /= mesh.Vsc();                 // 关键：除以每个 cell 的体积
}
```

上面是教学用的带注释版本，不是建议你复制替换 OpenFOAM 源码。
逐段解释如下。

#### `owner` 和 `neighbour`

对于每个内部面，OpenFOAM 保存两个 cell 编号：

```text
owner[facei]
    -> 这个面所属的第一侧 cell

neighbour[facei]
    -> 这个面另一侧的 cell
```

假设某个内部面上的有向通量为 $q_f$，OpenFOAM 的守恒加减就是：

```cpp
ivf[owner[facei]]     += q_f;
ivf[neighbour[facei]] -= q_f;
```

数学上：

$$
Q_{\mathrm{owner}} \mathrel{+}= q_f,
\qquad
Q_{\mathrm{neighbour}} \mathrel{-}= q_f.
$$

这正是“同一个内部面，对两侧 cell 贡献符号相反”。因此内部面的总贡献
在全域求和时会抵消。

#### 为什么边界面没有 `neighbour`

边界面只属于一个内部 cell，外面不是本计算域中的另一个 cell。因此边界
循环只有：

```cpp
ivf[pFaceCells[facei]] += pssf[facei];
```

边界的第二侧由边界条件提供，而不是由 `neighbour` 数组提供。对周期 patch，
OpenFOAM 通过周期 patch 的拓扑和边界场接口，把对应的另一侧信息接起来。

#### 最后一行就是 $1/V_c$

累加完成后，`ivf[c]` 还是：

$$
\mathcal Q_c=\sum_fF_{cf}T_f.
$$

源码最后执行：

```cpp
ivf /= mesh.Vsc();
```

这里是逐 cell 除法，因此得到：

$$
ivf[c]
=
\frac{\mathcal Q_c}{V_c}
=
\frac{1}{V_c}\sum_fF_{cf}T_f.
$$

这就是文档前面所说的：

```text
fvc::div(...)
    -> 不是未归一化通量和
    -> 而是已经除以 V_c 的散度
```

这也是为什么求解器中不能再写：

```cpp
T = T - deltaTDim*(residual/mesh.V());
```

那会把体积除两次，数学上变成：

$$
T_c^{n+1}
=
T_c^n
-
\Delta t
\frac{1}{V_c^2}
\sum_fF_{cf}T_f.
$$

这不是本题推导得到的离散格式。

### 24.5 `Gauss upwind` 在调用链中的位置

`fvc::div` 负责找到对流格式，但“面值 $T_f$ 如何计算”由具体的
对流格式完成。

以一阶迎风为例，逻辑可以写成：

```text
面通量 F_f 的符号
    |
    +-- F_f >= 0 -> 选择 owner 侧 T
    |
    +-- F_f < 0  -> 选择 neighbour 侧 T
    v
得到面值 T_f
    |
    v
形成数值面通量 F_f*T_f
    |
    v
surfaceIntegrate
    |
    v
除以 V_c
```

因此，`fvc::div` 不是“只把 `phi` 和 `T` 做普通乘法”。它还必须读取
`fvSchemes`，因为 $T_f$ 不是原始输入，而是由空间离散格式决定的。

对于：

```foam
div(phi,T) Gauss linearUpwind grad(T);
```

调用链多了一步梯度重构：

```text
T_c
    -> 根据 grad(T) 计算梯度
    -> 按 linearUpwind 重构 T_f
    -> 与 F_f 相乘
    -> 面积分
    -> 除以 V_c
```

所以改变 `fvSchemes` 可以改变空间离散方法，而不必重写学生版求解器
的时间推进代码。

## 25. `0.orig/U` 和 `0.orig/T` 不是“写了几个数”这么简单

初学者看到 `internalField` 里的长列表，容易产生一个误解：

> 这些数是不是求解器自己推导出来的？

对于当前案例，答案是：不是。它们是脚本根据题目给定的解析函数，在
真实网格 cell centre 上生成的初始离散数据。

### 25.1 `U` 文件逐项解释

第一案例中的核心部分是：

```foam
dimensions      [0 1 -1 0 0 0 0];
internalField   uniform (1 1 0);

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

逐行理解：

```text
dimensions
    -> 声明速度的物理量纲

uniform (1 1 0)
    -> 所有 cell 的速度都为 (1, 1, 0)

xMin/xMax/yMin/yMax cyclic
    -> 四个二维方向边界互相周期连接

zMin/zMax empty
    -> 这是二维计算，厚度方向不求解三维变化
```

速度的数学表达是：

$$
\boldsymbol U(x,y)=(1,1).
$$

但注意：文件中的 `uniform` 并不是说网格只有一个速度值。它表示
“用同一个向量初始化所有 cell”，所以从数据结构角度仍然是一个
`volVectorField`。

### 25.2 `T` 文件逐项解释

核心结构是：

```foam
dimensions      [0 0 0 0 0 0 0];
internalField   nonuniform List<scalar>
100
(
    ...
);
```

逐项解释：

```text
dimensions
    -> T 在当前案例中取无量纲标量

nonuniform
    -> 不同 cell 可以有不同数值

List<scalar>
    -> 后面是一串 scalar 数值

100
    -> 数值个数，必须等于网格 cell 总数

( ... )
    -> 按 OpenFOAM 的 cell 内部编号顺序逐个填写
```

如果 `N=10` 的四边形案例有 $10\times10=100$ 个 cell，那么这里必须有
100 个数。它们不应该按你眼睛看到的二维行列顺序猜测，而应该按照
OpenFOAM 网格的 cell 编号顺序写入。

脚本的数学步骤是：

$$
x_c=x_i,
\qquad
y_c=y_j,
\qquad
T_c^0=\sin\bigl(2\pi(x_c+y_c)\bigr).
$$

三角形网格不能简单使用规则的 $(i+0.5)/N$，因为三角形中心由真实网格
几何决定。正确顺序是：

```text
生成三角形网格
    -> 读取真实 cell centre C
    -> 对每个 (x_c,y_c) 计算 T_c^0
    -> 按 OpenFOAM cell 编号写入 T
```

这就是为什么三角形案例通常需要 `constant/C` 辅助文件，而四边形规则
案例可以直接由 $N$ 推出中心坐标。

## 26. 周期边界为什么能实现“跨边界搬运”

在数学上，周期边界不是把边界值简单设成一个常数，而是：

$$
T(0,y,t)=T(1,y,t),
\qquad
T(x,0,t)=T(x,1,t).
$$

对速度也必须满足相容的周期关系。当前案例中 `U` 是常量，因此这一点
自然满足。

从代码角度，周期边界通过三层共同完成：

```text
网格拓扑
    -> boundary 中声明成成对 cyclic patch

场文件
    -> U/T 的 boundaryField 声明 type cyclic

时间推进
    -> T.correctBoundaryConditions()
```

请注意最后一项的含义。它不是“重新求解 PDE”，而是让场对象根据其
边界 patch 类型，把边界数据更新到与当前内部场一致的状态。

真正的跨周期输运还依赖：

```text
mesh 的 owner/neighbour 和 patch 拓扑
fvc::flux 的面通量
fvc::div 的边界面处理
fvSchemes 的面值重构
```

因此只在 `T` 文件里写 `type cyclic` 还不够；网格本身也必须存在正确的
周期 patch 配对。

## 27. 为什么这个求解器没有 `fvSolution` 里的线性求解器

你可能会在标准 OpenFOAM 案例里看到：

```foam
solvers
{
    T
    {
        solver          smoothSolver;
        tolerance       1e-08;
        relTol          0;
    }
}
```

这些条目用于求解由 `fvm::ddt`、`fvm::div` 等接口组装出的线性方程组。

当前学生版求解器使用的是：

```cpp
T = T - deltaTDim*residual;
```

这条语句直接对每个 cell 做显式更新，没有形成：

```text
矩阵 A
未知向量 T
线性系统 A*T=b
```

所以当前 `fvSolution` 不负责线性求解。它存在的主要原因是让案例保留
标准 OpenFOAM 目录结构，而不是因为这里真的调用了一个迭代线性求解器。

把两类写法对照起来：

```text
显式：
    fvc::div(...)
    T = T - deltaT*residual
    不需要线性方程组求解器

隐式：
    fvm::ddt(...)
    fvm::div(...)
    solve(...)
    需要 fvSolution 中的线性求解器设置
```

这是“显式算法”和“隐式算法”的结构性区别，不只是函数名字不同。

## 28. 一次完整时间步的“纸笔版”

为了把所有对象连起来，假设你正在处理 cell $c$。一次时间步可以在纸上
写成下面五步。

### 第一步：确定网格几何

从 `constant/polyMesh` 得到：

$$
V_c,\qquad \boldsymbol S_{cf},\qquad
\text{owner/neighbour}.
$$

代码层面由：

```cpp
#include "createMesh.H"
```

创建的 `fvMesh mesh` 保存这些信息。

### 第二步：得到面通量

速度文件提供 cell 上的 $\boldsymbol U_c$。OpenFOAM 根据网格和面评价得到
$\boldsymbol U_f$，然后：

$$
F_{cf}=\boldsymbol U_f\cdot\boldsymbol S_{cf}.
$$

代码：

```cpp
surfaceScalarField phi
(
    IOobject(...),
    fvc::flux(U)
);
```

### 第三步：根据格式得到面标量

由 `T` 的 cell 值和 `fvSchemes`：

$$
T_f=
\begin{cases}
T_{\mathrm{upwind}}, & \text{upwind},\\
T_{\mathrm{linearUpwind}}, & \text{linearUpwind}.
\end{cases}
$$

代码入口：

```cpp
fvc::div(phi, T, "div(phi,T)")
```

### 第四步：组装并归一化

先形成：

$$
\mathcal Q_c=\sum_fF_{cf}T_f.
$$

再除以控制体积：

$$
R_c=\frac{\mathcal Q_c}{V_c}.
$$

### 第五步：前向 Euler 更新

$$
T_c^{n+1}=T_c^n-\Delta tR_c^n.
$$

代码：

```cpp
const dimensionedScalar deltaTDim
(
    "deltaT",
    dimTime,
    stepDeltaT
);

T = T - deltaTDim*residual;
```

如果这五步你能够逐步写出来，就已经真正理解了这个求解器的数学骨架。

## 29. 教授给初学者的阅读方法

学习时建议每次只看一个对象，不要同时追所有源码。

### 第一轮：只看物理量位置

回答：

```text
T 在 cell 还是 face？
U 在 cell 还是 face？
phi 在 cell 还是 face？
residual 在 cell 还是 face？
```

标准答案：

```text
T         -> cell
U         -> cell
phi       -> face
residual  -> cell
```

### 第二轮：只看数据流

```text
U + mesh
    -> phi

phi + T + fvSchemes
    -> residual

residual + deltaT
    -> new T
```

### 第三轮：只看量纲

对于本题：

```text
U
    -> m/s

Sf
    -> m^2

phi = U dot Sf
    -> m^3/s

phi*T/V
    -> 1/s * T

deltaT*residual
    -> T
```

量纲检查是非常有力量的自检方法。如果某一步的量纲对不上，通常说明
公式、字段类型或体积因子处理错了。

### 第四轮：只看离散格式

只改变这一行：

```foam
div(phi,T)      Gauss upwind;
```

变为：

```foam
div(phi,T)      Gauss linearUpwind grad(T);
```

然后问：

```text
求解器的主循环有没有改变？
phi 的定义有没有改变？
改变的是哪一步？
```

答案是：主要改变面值 $T_f$ 的重构方式，时间推进骨架不变。

## 30. 最终自测题

在你认为自己已经理解前，建议不看答案，先自己回答：

1. 为什么 `T` 是 `volScalarField` 而不是 `surfaceScalarField`？
2. 为什么 `phi` 是 `surfaceScalarField`？
3. `fvc::flux(U)` 中为什么不传入 `T`？
4. `"div(phi,T)"` 是数值还是字典键？
5. `fvSchemes` 中哪一行决定 $T_f$ 的构造？
6. `fvc::div` 结果为什么已经除以 $V_c$？
7. `owner` 和 `neighbour` 如何保证内部面守恒？
8. 为什么周期边界不是简单的 `fixedValue`？
9. 为什么当前求解器不需要 `solve()`？
10. `residual` 和 L1 error 的区别是什么？

参考答案：

```text
1. T 是每个 cell 的被输运标量。
2. phi 是每个 face 的有向体积通量。
3. fvc::flux(U) 只负责生成速度通量，T 在散度算子中参与面值重构。
4. "div(phi,T)" 是查找离散格式的字符串键。
5. divSchemes 中的 div(phi,T) 条目。
6. fvc::surfaceIntegrate 的最后一步执行逐 cell 的体积归一化。
7. 同一个内部面给 owner 加、给 neighbour 减。
8. 周期边界需要把两侧拓扑和场值配成相同的周期接口。
9. 当前算法直接进行前向 Euler 显式更新，没有组装线性系统。
10. residual 是离散 PDE 的输运项；L1 error 是数值解与精确解的差。
```

当你能不看代码回答这些问题时，再去读 OpenFOAM 源码会轻松很多。此时
源码中的模板、`tmp`、patch 和 field 类不再是孤立的 API，而是已经有了
明确的数学位置。
