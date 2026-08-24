# 公式到代码的映射

这份文档专门解释第 1 题线性对流方程中，数学对象怎样变成 OpenFOAM
代码对象。重点不是背类名，而是理解三个层次之间的关系：

```text
连续数学对象  ->  有限体积离散对象  ->  OpenFOAM C++ 类型
```

本工程目前处在教学骨架阶段。Stage 0 已经能读取 `mesh`、`U` 和 `T`，
Stage 1 已经加入面通量，Stage 2 已经加入 CFL 候选时间步，
Stage 3 已经加入体积归一化残差，Stage 4 已经加入单步前向 Euler 更新。
当前仍未实现完整时间循环和误差后处理。

## 1. 总问题

题目中的守恒型线性对流方程是：

$$\frac{\partial \phi}{\partial t}+\nabla\cdot(\boldsymbol{u}\phi)=0.$$

本工程中使用 OpenFOAM 字段名：

| 数学记号 | 工程变量名 | 含义 |
|---|---|---|
| $\boldsymbol{u}$ | `U` | 给定速度场 |
| $\phi$ | `T` | 被输运的标量场 |
| $\Omega$ | `mesh` | 计算区域及其离散网格 |
| $t$ | `runTime` | 时间对象 |

这里的“守恒型”很重要。我们离散的是：

$$\nabla\cdot(\boldsymbol{u}\phi).$$

不是直接离散：

$$\boldsymbol{u}\cdot\nabla\phi.$$

当速度场不是严格无散时，这两个形式一般不等价。有限体积法天然适合守恒型方程，
因为它直接计算控制体边界上的通量。

## 2. 从连续区域到有限体积网格

连续数学问题定义在区域 $\Omega$ 上。有限体积法先把区域拆成很多小控制体：

$$\Omega=\bigcup_c\Omega_c.$$

每个控制体 $\Omega_c$ 有：

| 数学对象 | 含义 | OpenFOAM 中在哪里 |
|---|---|---|
| $\Omega_c$ | 第 $c$ 个控制体 | `fvMesh mesh` 的一个 cell |
| $V_c$ | 单元体积 | `mesh.V()` |
| $\partial\Omega_c$ | 单元边界 | cell 周围的 faces |
| $f$ | 单元的某个面 | `mesh` 中的 face |
| $\boldsymbol{S}_{cf}$ | 面面积向量 | `mesh.Sf()`，方向和 owner/neighbour 约定有关 |

`fvMesh` 是 OpenFOAM 有限体积网格类型。它不是单纯的点坐标表，而是一个完整的
有限体积拓扑对象，包含：

- 单元数量；
- 面数量；
- 点坐标；
- 单元体积；
- 面面积向量；
- owner/neighbour 关系；
- 边界 patch 类型。

所以当 C++ 里出现：

```cpp
#include "createMesh.H"
```

它实际是在创建：

```cpp
fvMesh mesh;
```

数学上，这一步就是完成：

$$\Omega\longrightarrow\{\Omega_c\}.$$

Stage 0 只需要确认这件事成功。`checkMesh` 输出 `Mesh OK`，说明这个离散区域
在几何和拓扑上可以用。

## 3. OpenFOAM 类型总览

这一节先把后面会遇到的类型全部放在一张表里。

| OpenFOAM 类型 | 数学角色 | 定义位置 | 例子 |
|---|---|---|---|
| `scalar` | 一个实数 | 单个数值 | `0.2`、`deltaT` |
| `label` | 一个整数编号 | 索引 | cell 编号、face 编号 |
| `word` | 一个短字符串 | 名字 | `"U"`、`"T"` |
| `vector` | 三维向量 | 单个点或单个单元上的向量 | `(1 1 0)` |
| `scalarField` | 一列实数 | 内存数组 | 所有单元体积 |
| `vectorField` | 一列向量 | 内存数组 | 所有面面积向量 |
| `dimensionedScalar` | 带量纲的实数 | 时间步等有物理单位的数 | $\Delta t$ |
| `dictionary` | 字典文件 | 运行参数 | `controlDict` |
| `IOobject` | 读写说明 | 字段文件和网格关联 | 读取 `0/U` |
| `Time` | 时间管理器 | 全局运行对象 | `runTime` |
| `fvMesh` | 有限体积网格 | cell/face/patch 拓扑 | `mesh` |
| `volScalarField` | 单元标量场 | 每个 cell 一个标量 | `T` |
| `volVectorField` | 单元向量场 | 每个 cell 一个向量 | `U` |
| `surfaceScalarField` | 面标量场 | 每个 face 一个标量 | `phi` |
| `surfaceVectorField` | 面向量场 | 每个 face 一个向量 | 概念上的 $\boldsymbol{U}_f$ |
| `tmp<volScalarField>` | 临时单元标量场 | 算子返回的临时结果 | `fvc::div(...)` |

这些类型不是随便起的名字。它们基本由两个信息组合而成：

```text
位置 + 值的种类
```

例如：

| 类型 | 位置 | 值的种类 |
|---|---|---|
| `volScalarField` | volume/cell 上 | scalar |
| `volVectorField` | volume/cell 上 | vector |
| `surfaceScalarField` | face/surface 上 | scalar |
| `surfaceVectorField` | face/surface 上 | vector |

这是理解 OpenFOAM 的第一把钥匙。

## 4. `scalar`、`label`、`word`：最基础的三类量

### 4.1 `scalar`

`scalar` 可以先理解为双精度实数。数学中一个普通实数，比如：

$$\Delta t,\quad \mathrm{Co}_{\max},\quad \lambda_{\max}.$$

在 OpenFOAM 中通常写成：

```cpp
scalar deltaT;
scalar maxCo;
scalar rate;
```

如果某个量只是一个全局数，而不是每个单元一个值、每个面一个值，它通常就是
`scalar`。

### 4.2 `label`

`label` 是整数编号。数学中我们常写：

$$c=1,2,\ldots,N_c,\qquad f=1,2,\ldots,N_f.$$

这里的 $c$、$f$ 若出现在 OpenFOAM 索引里，通常就是 `label`。

本求解器目前没有显式写循环访问 cell/face，所以你暂时不需要手写 `label`。
但当你以后手写 owner/neighbour 上风格式时，就会用到它。

### 4.3 `word`

`word` 是 OpenFOAM 的短字符串类型。它常用于字段名和字典关键字。

例如：

```cpp
const word velocityName("U");
const word advectedName("T");
```

在本工程中：

```cpp
controlDict.lookupOrDefault<word>("velocityField", "U")
```

含义是：从 `controlDict` 读一个名字。如果用户没有写 `velocityField`，默认字段名
就是 `U`。

数学上，`word` 本身没有公式含义。它是把数学对象和文件名连接起来的标签。

## 5. `dictionary`：运行参数的数学意义

OpenFOAM 的 `dictionary` 对应一个字典文件，例如：

```cpp
const dictionary& controlDict = runTime.controlDict();
```

它读取的是：

```text
cases/01_sine_wave_quad/N20/system/controlDict
```

数学上，`controlDict` 不是未知函数，也不是算子；它是计算任务的参数集合。比如：

| `controlDict` 条目 | 数学或数值意义 |
|---|---|
| `startTime 0` | 初始时间 $t_0=0$ |
| `endTime 1` | 终止时间 $T_{\mathrm{end}}=1$ |
| `deltaT 0.01` | 默认时间步 |
| `maxCo 0.2` | 目标 Courant 数 |
| `velocityField U` | 速度场名字 |
| `advectedField T` | 被输运标量名字 |

后续 Stage 2 中，`maxCo` 会进入 CFL 公式：

$$\Delta t=\frac{2\,\mathrm{Co}_{\mathrm{target}}}{\max_c\left(\sum_f|F_{cf}|/V_c\right)}.$$

所以 `dictionary` 的作用是：把“数值实验参数”从代码中分离出来，放在 case 文件里。

## 6. `Time runTime`：时间变量和输出控制

数学中，时间变量是：

$$t.$$

时间层是：

$$t^n,\quad t^{n+1}=t^n+\Delta t.$$

OpenFOAM 中的 `Time runTime` 管理这些东西：

- 当前时间；
- 当前时间目录名；
- 终止时间；
- 时间步；
- 什么时候写文件；
- 当前 case 的根目录；
- `system/controlDict`。

在当前 C++ 骨架中，`runTime` 来自：

```cpp
#include "createTime.H"
```

常见用法：

| 代码 | 含义 |
|---|---|
| `runTime.name()` | 当前时间目录名，例如 `0` |
| `runTime.value()` | 当前物理时间数值 |
| `runTime.endTime().value()` | 终止时间 |
| `runTime.setDeltaT(deltaT)` | 设置下一步时间步长 |
| `runTime++` | 时间推进到下一层 |
| `runTime.write()` | 按输出规则写场文件 |

Stage 0 只用到了 `runTime.name()`。它告诉 `IOobject`：现在应该从哪个时间目录读字段。

例如：

```cpp
IOobject
(
    "T",
    runTime.name(),
    mesh,
    IOobject::MUST_READ,
    IOobject::AUTO_WRITE
)
```

如果当前 `runTime.name()` 是 `0`，那么它会读：

```text
0/T
```

## 7. `fvMesh mesh`：控制体、面和边界

有限体积法的基本公式是对每个控制体写的：

$$\int_{\Omega_c}\frac{\partial\phi}{\partial t}\,\mathrm d\Omega+\int_{\Omega_c}\nabla\cdot(\boldsymbol{u}\phi)\,\mathrm d\Omega=0.$$

这里所有几何对象都来自 `fvMesh mesh`。

| 公式对象 | OpenFOAM 表达 | 类型理解 |
|---|---|---|
| $\Omega_c$ | 第 $c$ 个 cell | `mesh` 的一个控制体 |
| $V_c$ | `mesh.V()` | 所有单元体积组成的 `scalarField` |
| $f\in\partial\Omega_c$ | cell 周围的 faces | 网格拓扑关系 |
| $\boldsymbol{S}_f$ | `mesh.Sf()` | 所有面面积向量组成的面向量场 |
| 边界 $\partial\Omega$ | `mesh.boundary()` | patch 列表 |

`blockMesh` 根据 `system/blockMeshDict` 生成：

```text
constant/polyMesh/
```

`createMesh.H` 再把这个目录读进 C++，形成 `fvMesh mesh`。

因此文件和代码的关系是：

```text
system/blockMeshDict
    -> blockMesh
        -> constant/polyMesh/*
            -> createMesh.H
                -> fvMesh mesh
```

在数学上就是：

$$\Omega\rightarrow\Omega_h\rightarrow\{\Omega_c,V_c,\boldsymbol{S}_f\}.$$

## 8. `volScalarField T`：单元中心标量场

`volScalarField` 的意思是：

```text
volume + scalar + field
```

也就是每个控制体一个标量值。

数学上：

$$T(\boldsymbol{x},t)\rightarrow\{T_c^n\}_{c=1}^{N_c}.$$

在本题里，`T` 是被输运的标量，等价于题目里的 $\phi$。

Stage 0 中读取 `T` 的代码是：

```cpp
volScalarField T
(
    IOobject
    (
        advectedName,
        runTime.name(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);
```

逐项解释：

| 代码 | 含义 |
|---|---|
| `volScalarField T` | 创建一个单元标量场 |
| `advectedName` | 字段文件名，默认 `T` |
| `runTime.name()` | 时间目录，当前为 `0` |
| `mesh` | 这个场定义在哪个网格上 |
| `MUST_READ` | 必须从磁盘读取 |
| `AUTO_WRITE` | 后续时间推进后可以自动写出 |

它读取的文件是：

```text
cases/01_sine_wave_quad/N20/0/T
```

文件中 `internalField` 存的是所有内部单元的 $T_c^0$。

所以：

| 文件内容 | 数学含义 |
|---|---|
| `internalField` | $\{T_c^0\}$ |
| `boundaryField` | 边界上的 $T$ 如何处理 |
| `dimensions [0 0 0 0 0 0 0]` | $T$ 是无量纲量 |

## 9. `volVectorField U`：单元中心速度场

`volVectorField` 的意思是：

```text
volume + vector + field
```

也就是每个控制体一个向量值。

数学上：

$$\boldsymbol{u}(\boldsymbol{x},t)\rightarrow\{\boldsymbol{U}_c^n\}_{c=1}^{N_c}.$$

本算例中速度是常量：

$$\boldsymbol{U}=(1,1,0).$$

对应文件：

```text
cases/01_sine_wave_quad/N20/0/U
```

里面写：

```foam
dimensions      [0 1 -1 0 0 0 0];
internalField   uniform (1 1 0);
```

这里的量纲 `[0 1 -1 0 0 0 0]` 表示：

$$L\,T^{-1}.$$

也就是速度单位。

`U` 为什么是 `volVectorField` 而不是 `surfaceVectorField`？因为输入文件给的是单元速度。
后续计算通量时，OpenFOAM 会把单元速度插值到面上，再与面面积向量点乘。

## 10. `surfaceScalarField phi`：面通量

这是 Stage 1 最重要的类型。

有限体积法中，守恒通量不是定义在单元中心，而是定义在控制体边界面上：

$$F_{cf}=\boldsymbol{u}_f\cdot\boldsymbol{S}_{cf}.$$

这里：

- $\boldsymbol{u}_f$ 是面上的速度；
- $\boldsymbol{S}_{cf}$ 是面面积向量；
- 点乘以后得到一个标量；
- 每个面都有一个这样的标量。

所以代码类型是：

```cpp
surfaceScalarField phi;
```

名称 `phi` 在 OpenFOAM 里常用来表示面体积通量。不要和题目里的未知量 $\phi$
混淆。在本工程里：

| 符号 | 含义 |
|---|---|
| 题目 $\phi$ | 被输运标量，本工程叫 `T` |
| OpenFOAM `phi` | 面通量，本工程后续用于 $F_f$ |

后续 Stage 1 会写：

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

它对应公式：

$$F_f=\boldsymbol{U}_f\cdot\boldsymbol{S}_f.$$

`NO_READ` 表示它不从磁盘读取，而是由代码计算出来。
`NO_WRITE` 表示第一版不把它写出到时间目录。

## 11. `surfaceVectorField`：为什么文档里说它是“概念值”

公式里有：

$$\boldsymbol{u}_f.$$

这是面上的速度，是向量，并且每个面一个值。

如果我们显式保存它，它会接近：

```cpp
surfaceVectorField Uf;
```

但在本求解器第一版里，我们不单独保存 `Uf`。OpenFOAM 的：

```cpp
fvc::flux(U)
```

内部会完成“把 `volVectorField U` 插值到面上，然后与面面积向量点乘”这件事。

所以文档中说：

```text
surfaceVectorField 是概念上的 U_f
```

意思是数学上它存在，但代码里第一版不需要单独声明。

## 12. `IOobject`：字段怎样从文件变成 C++ 对象

OpenFOAM 的一个字段对象不仅仅是一组数值，还必须知道：

- 自己叫什么名字；
- 从哪个时间目录读取；
- 属于哪个网格；
- 是否必须读取；
- 是否需要写出。

这些信息由 `IOobject` 提供。

例如读取 `T`：

```cpp
IOobject
(
    advectedName,
    runTime.name(),
    mesh,
    IOobject::MUST_READ,
    IOobject::AUTO_WRITE
)
```

对应关系：

| `IOobject` 参数 | 实际作用 |
|---|---|
| `advectedName` | 文件对象名，例如 `T` |
| `runTime.name()` | 时间目录，例如 `0` |
| `mesh` | 注册到哪个网格对象 |
| `MUST_READ` | 文件必须存在 |
| `AUTO_WRITE` | 后续写时间目录时自动输出 |

如果读取失败，通常不是数学错，而是文件组织错，比如：

- `0/T` 不存在；
- 文件头 `class` 写错；
- 边界 patch 名称和网格不一致；
- 量纲不匹配。

## 13. `dimensionedScalar`：带单位的数

普通 `scalar` 只是一个数，比如：

```cpp
scalar deltaT = 0.01;
```

但 OpenFOAM 的场运算会检查量纲。例如：

$$T^{n+1}=T^n-\Delta t\,R^n.$$

其中：

- $T$ 是无量纲；
- $R$ 的量纲是 $1/s$；
- $\Delta t$ 的量纲是 $s$；
- 所以 $\Delta t R$ 是无量纲，可以和 $T$ 相减。

为了让 OpenFOAM 知道 $\Delta t$ 带时间量纲，后续会使用：

```cpp
dimensionedScalar deltaTDim("deltaT", dimTime, deltaT);
```

对应数学量：

$$\Delta t.$$

`dimensionedScalar` 的作用是防止你把时间、长度、速度、无量纲量混在一起误算。
这对 CFD 程序很重要，因为很多错误编译能过，但量纲上根本不对。

## 14. `tmp<volScalarField>`：临时场

OpenFOAM 的很多显式算子会返回临时对象。例如：

```cpp
tmp<volScalarField> tResidual
(
    fvc::div(phi, T, "div(phi,T)")
);
```

它对应体积归一化后的对流残差：

$$R_c^n=\frac{1}{V_c}\sum_fF_{cf}^nT_f^n.$$

如果手动写 owner/neighbour 面循环，可以把未除体积的面通量和临时记为：

$$\mathcal Q_c^n=\sum_fF_{cf}^nT_f^n.$$

因此：

$$R_c^n=\frac{\mathcal Q_c^n}{V_c}.$$

为什么是 `volScalarField`？

- `fvc::div(...)` 计算散度；
- 散度结果回到每个控制体；
- 被输运量 `T` 是标量；
- 所以每个单元得到一个体积归一化后的标量散度。

为什么外面有 `tmp<...>`？

因为这个空间散度场通常只是中间结果，用完就可以释放。`tmp` 是 OpenFOAM 的临时对象管理器，
用于减少不必要的复制和内存开销。

你可以先把它理解成：

```text
tmp<volScalarField> = 一个临时的 volScalarField
```

等你熟悉以后，再深入理解 OpenFOAM 的内存管理。

## 15. `fvc` 和 `fvm`：显式算子与隐式算子

这两个名字非常关键。

| 前缀 | 含义 | 返回结果 | 本题是否使用 |
|---|---|---|---|
| `fvc` | finite volume calculus | 直接算出一个场 | 使用 |
| `fvm` | finite volume method matrix | 组装线性方程矩阵 | 不使用 |

本题要求显式离散，因此体积归一化后的空间散度可使用：

```cpp
fvc::div(phi, T, "div(phi,T)")
```

它直接给出已经除以控制体体积的 $R^n$。

对应的源码调用链是：

```text
fvcDiv.C
    -> 选择 fvSchemes 里的离散格式
    -> gaussConvectionScheme.C
    -> fvcSurfaceIntegrate.C
    -> ivf /= mesh.Vsc()
```

不要写成：

```cpp
fvm::div(phi, T)
```

因为 `fvm::div` 会进入隐式矩阵形式，不再是你要手写的前向 Euler 显式更新。

## 16. 控制体积分公式中的类型关系

从连续方程开始：

$$\frac{\partial T}{\partial t}+\nabla\cdot(\boldsymbol{U}T)=0.$$

对控制体积分：

$$\int_{\Omega_c}\frac{\partial T}{\partial t}\,\mathrm d\Omega+\int_{\Omega_c}\nabla\cdot(\boldsymbol{U}T)\,\mathrm d\Omega=0.$$

第一项：

$$\int_{\Omega_c}\frac{\partial T}{\partial t}\,\mathrm d\Omega\approx V_c\frac{\mathrm dT_c}{\mathrm dt}.$$

这里的类型映射是：

| 公式 | 类型 |
|---|---|
| $T_c$ | `volScalarField T` 的一个单元值 |
| $V_c$ | `mesh.V()` 的一个元素 |
| $\mathrm dT_c/\mathrm dt$ | 后续由时间更新近似 |

第二项用高斯散度定理：

$$\int_{\Omega_c}\nabla\cdot(\boldsymbol{U}T)\,\mathrm d\Omega=\sum_{f\in\partial\Omega_c}F_{cf}T_f.$$

类型映射是：

| 公式 | 类型 |
|---|---|
| $F_{cf}$ | `surfaceScalarField phi` |
| $T_f$ | `fvSchemes` 根据 `T` 和 `phi` 选择的面值 |
| $\sum_f$ | `fvc::div` 内部根据 mesh 拓扑完成 |
| 除以 $V_c$ | `fvc::div` 返回的是体积归一化后的残差 $R_c$ |

因此半离散形式：

$$V_c\frac{\mathrm dT_c}{\mathrm dt}+\sum_fF_{cf}T_f=0.$$

在 OpenFOAM 中更接近：

```cpp
fvc::div(phi, T, "div(phi,T)")
```

因为 `fvc::div` 返回的是：

$$\frac{1}{V_c}\sum_fF_{cf}T_f.$$

## 17. 一阶迎风公式中的类型关系

一阶迎风面值：

$$\phi_f=\begin{cases}\phi_{\mathrm{owner}},&F_f\geq0,\\\phi_{\mathrm{neighbour}},&F_f<0.\end{cases}$$

在本工程把题目中的 $\phi$ 命名为 `T`，所以代码中的对应关系是：

| 公式对象 | OpenFOAM 对象 |
|---|---|
| $\phi_{\mathrm{owner}}$ | `T` 在 owner 单元的值 |
| $\phi_{\mathrm{neighbour}}$ | `T` 在 neighbour 单元的值 |
| $F_f$ | `surfaceScalarField phi` 在该面的值 |

本工程第一版不手写 owner/neighbour，而是在 `fvSchemes` 中声明：

```foam
div(phi,T)      Gauss upwind;
```

然后在 C++ 中调用同名离散项：

```cpp
fvc::div(phi, T, "div(phi,T)")
```

关系是：

```text
surfaceScalarField phi 的正负
    -> 决定上游方向
        -> fvSchemes 的 upwind 选取 T_f
            -> fvc::div 求每个 cell 的散度
```

## 18. 前向 Euler 公式中的类型关系

按照新的 `.tex` 记号，半离散后：

$$\frac{\mathrm dT_c}{\mathrm dt}=-R_c,\qquad R_c=\frac{1}{V_c}\sum_fF_{cf}T_f.$$

前向 Euler：

$$T_c^{n+1}=T_c^n-\Delta tR_c^n.$$

类型对应：

| 公式对象 | OpenFOAM 类型 |
|---|---|
| $T_c^n$ | `volScalarField T` 当前 internal field |
| $R_c^n$ | `tmp<volScalarField>`，由 `fvc::div` 返回，已经除以体积 |
| $\mathcal Q_c^n=V_cR_c^n$ | 手动面循环时的未归一化通量和，OpenFOAM 路线中通常不显式保存 |
| $\Delta t$ | `dimensionedScalar` 或 `scalar` |
| $T_c^{n+1}$ | 修改后的 `T.internalFieldRef()` |

Stage 4 现在在源码里就是这样写的：

```cpp
const dimensionedScalar deltaTDim("deltaT", dimTime, deltaT);

T = T - deltaTDim*residual;
T.correctBoundaryConditions();
```

这段代码的数学顺序是：

```text
先算 R^n = (1/V) sum_f F*T_f
再算 deltaT*R^n
再从 T^n 中减掉
最后更新边界条件
```

显式方法的关键是右端项必须来自旧时间层 $n$。

### 18.1 如果改成未除体积的原始面通量和

如果你想自己写面循环，而不是直接调用 `fvc::div`，那更像这样：

```cpp
scalarField fluxSum(mesh.nCells(), 0.0);

forAll(owner, facei)
{
    const label own = owner[facei];
    const label nei = neighbour[facei];
    const scalar F = phi[facei];

    const scalar Tf = (F >= 0 ? T[own] : T[nei]);
    const scalar H  = F * Tf;

    fluxSum[own] += H;
    fluxSum[nei] -= H;
}
```

这时 `fluxSum` 对应的是未除体积的原始面通量和：

$$\mathcal Q_c=\sum_{f\in\partial\Omega_c}F_{cf}T_f.$$

如果要得到最终残差，仍然要再除一次体积：

$$R_c=\frac{\mathcal Q_c}{V_c}.$$

所以两种写法的区别是：

| 写法 | 结果 | 是否已除体积 |
|---|---|---|
| `fvc::div(phi, T, "div(phi,T)")` | `tResidual` / `R_c` | 是 |
| 手写 `fluxSum` | $\mathcal Q_c$ | 否 |

如果你用 `fluxSum` 直接更新，就必须写成：

$$T_c^{n+1}=T_c^n-\frac{\Delta t}{V_c}\mathcal Q_c^n.$$

如果你先做归一化，再更新，就回到当前主文档的写法：

$$T_c^{n+1}=T_c^n-\Delta tR_c^n.$$

## 19. CFL 公式中的类型关系

CFL 控制时间步：

$$\mathrm{Co}_{\max}=\frac{\Delta t}{2}\max_c\left(\frac{\sum_f|F_{cf}|}{V_c}\right).$$

令：

$$\lambda_c=\frac{\sum_f|F_{cf}|}{V_c}.$$

则：

$$\Delta t=\frac{2\,\mathrm{Co}_{\mathrm{target}}}{\max_c\lambda_c}.$$

类型对应：

| 公式对象 | OpenFOAM 类型 |
|---|---|
| $F_{cf}$ | `surfaceScalarField phi` |
| $|F_{cf}|$ | `mag(phi)` |
| $\sum_f|F_{cf}|$ | `fvc::surfaceSum(mag(phi))` |
| $V_c$ | `mesh.V()` |
| $\lambda_c$ | 临时的 cell 上 `scalarField` |
| $\max_c\lambda_c$ | `scalar` |
| $\Delta t$ | `scalar`，用于 `runTime.setDeltaT(deltaT)` |

这里有一个重要的空间位置变化：

```text
面上的 phi
    -> 对每个 cell 周围的面求和
        -> 得到 cell 上的 lambda_c
            -> 全局最大值得到一个 scalar
```

这就是 `surfaceScalarField`、`volScalarField` 和 `scalar` 之间的关系。

## 20. 运行文件和 C++ 类型的关系

OpenFOAM 的字段不是只在 C++ 里存在。它们也对应磁盘文件。

| 磁盘文件 | C++ 类型 | 数学对象 |
|---|---|---|
| `0/U` | `volVectorField U` | $\boldsymbol{U}_c$ |
| `0/T` | `volScalarField T` | $T_c^0$ |
| `constant/polyMesh/*` | `fvMesh mesh` | $\Omega_h$ |
| `system/controlDict` | `dictionary` | 时间和运行参数 |
| `system/fvSchemes` | `dictionary` | 离散格式规则 |
| `system/fvSolution` | `dictionary` | 线性求解器规则 |

Stage 0 的读取链条是：

```text
system/controlDict
    -> createTime.H
        -> Time runTime

constant/polyMesh/*
    -> createMesh.H
        -> fvMesh mesh

0/U
    -> IOobject + mesh
        -> volVectorField U

0/T
    -> IOobject + mesh
        -> volScalarField T
```

## 21. 类型之间的整体关系图

第 1 题从输入到更新的关系可以这样看：

```text
Time runTime
    |
    |-- 读取 controlDict: endTime, maxCo, field names
    |
fvMesh mesh
    |
    |-- mesh.V()  -> cell volume V_c
    |-- mesh.Sf() -> face area vector S_f
    |
volVectorField U
    |
    |-- fvc::flux(U)
    v
surfaceScalarField phi
    |
    |-- fvc::div(phi, T, "div(phi,T)")
    v
tmp<volScalarField> R
    |
volScalarField T
    |
    |-- T = T - deltaT*R
    v
new volScalarField T
```

这张图可以作为你读 C++ 的主线。

## 22. 目前已经做到哪里

当前 `explicitAdvectionFoamStudent.C` 已完成：

```text
runTime
mesh
U
T
phi
maxCo
sumPhi
rate
rateMax
deltaT
residual
deltaTDim
updated T
```

也就是：

$$
t,\quad
\Omega_h,\quad
\{\boldsymbol{U}_c\},\quad
\{T_c^0\},\quad
\{F_f\},\quad
\left\{\frac{\sum_f|F_{cf}|}{V_c}\right\},\quad
\Delta t_{\mathrm{CFL}},\quad
\{R_c^n\},\quad
\{T_c^{n+1}\}.
$$

它已经创建并执行了：

```text
tmp<volScalarField> R
dimensionedScalar deltaTDim
T = T - deltaT*R
```

但还没有执行完整时间循环：

```text
runTime++
runTime.write()
```

Stage 4 的单步公式已经进入源码：

$$
T_c^{n+1}=T_c^n-\Delta tR_c^n.
$$

当前仍然没有完整时间循环，因此还没有重复执行

```text
R^n -> T^(n+1) -> runTime++ -> write
```

并推进到 `endTime`。

## 23. 学习时的判断口诀

读一个 OpenFOAM 类型时，先问四个问题：

1. 它是一个数，还是一列数，还是一个场？
2. 它定义在 cell 上，还是 face 上？
3. 它是 scalar，还是 vector？
4. 它是从文件读入，还是由公式计算出来？

例如：

| 对象 | 判断 |
|---|---|
| `T` | cell 上的 scalar field，从 `0/T` 读取 |
| `U` | cell 上的 vector field，从 `0/U` 读取 |
| `phi` | face 上的 scalar field，由 `fvc::flux(U)` 计算 |
| `mesh.V()` | cell 上的 scalar list，由网格几何给出 |
| `deltaT` | 一个全局 scalar，由 CFL 公式给出 |
| `tResidual` | cell 上的临时 scalar field，由 `fvc::div` 给出 |

只要这四个问题答清楚，公式到代码的映射基本就不会乱。
