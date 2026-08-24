# Stage 0：代码、文件与数学对象

本文档解释学生版求解器当前已经完成的 Stage 0。

Stage 0 的目标不是求解线性对流方程，而是完成计算问题的初始化：

```text
读取运行参数
    -> 读取时间对象
        -> 读取网格
            -> 读取速度场 U
                -> 读取标量场 T
                    -> 检查量纲和数值范围
```

只有这些对象全部准备好，后面才有可能计算面通量、CFL 和显式更新。

## 1. Stage 0 不做什么

当前 Stage 0 还没有计算：

$$F_f=\boldsymbol{U}_f\cdot\boldsymbol{S}_{cf}.$$

也没有计算：

$$R_c^n=\frac{1}{V_c}\sum_{f\in\partial\Omega_c}F_{cf}T_f^n.$$

更没有进行时间更新：

$$T_c^{n+1}=T_c^n-\Delta tR_c^n.$$

当前程序只做：

- 找到 case；
- 读取 `controlDict`；
- 创建时间对象；
- 创建有限体积网格；
- 读取 `U`；
- 读取 `T`；
- 检查字段的量纲和范围；
- 正常退出。

因此日志中出现：

```text
No numerical update has been implemented yet.
```

不是错误，而是说明当前教学阶段有意停在这里。

## 2. 连续问题与离散问题

连续方程是：

$$\frac{\partial\phi}{\partial t}
+\nabla\cdot(\boldsymbol{u}\phi)=0.$$

连续问题中的对象是：

- 连续区域 $\Omega$；
- 连续时间 $t$；
- 连续速度场 $\boldsymbol{u}(\boldsymbol{x},t)$；
- 连续标量场 $\phi(\boldsymbol{x},t)$。

有限体积法先把区域划分成控制体：

$$\Omega=\bigcup_{c=1}^{N_c}\Omega_c.$$

然后把连续场用单元值表示：

$$
\boldsymbol{u}(\boldsymbol{x},t)
\longrightarrow
\{\boldsymbol{U}_c^n\}_{c=1}^{N_c},
$$

$$
\phi(\boldsymbol{x},t)
\longrightarrow
\{T_c^n\}_{c=1}^{N_c}.
$$

本工程中将题目中的 $\phi$ 命名为 `T`，因此：

```text
数学中的 phi      -> 代码中的 T
数学中的 u        -> 代码中的 U
```

## 3. 求解器的头文件

当前源码的开头是：

```cpp
#include "argList.H"
#include "volFields.H"
```

### 3.1 `argList.H`

```cpp
#include "argList.H"
```

它提供 OpenFOAM 程序的命令行参数处理能力。

例如用户可以运行：

```bash
explicitAdvectionFoamStudent
```

也可以运行：

```bash
explicitAdvectionFoamStudent -case /path/to/case
```

还可以运行：

```bash
explicitAdvectionFoamStudent -help
```

`argList.H` 本身不对应某一条 PDE 公式。它属于程序基础设施，作用是让求解器知道：

- 当前使用哪个 case；
- 是否并行运行；
- 是否需要显示帮助；
- 是否指定了其他运行参数。

### 3.2 `volFields.H`

```cpp
#include "volFields.H"
```

它提供体场类型，例如：

```cpp
volScalarField
volVectorField
```

这里的 `vol` 表示 volume，意思是场定义在控制体上。

所以：

```text
volScalarField = 每个控制体一个标量
volVectorField = 每个控制体一个向量
```

本阶段：

```cpp
volVectorField U;
volScalarField T;
```

正是通过这个头文件获得类型定义。

注意，当前源码还没有包含：

```cpp
#include "fvcFlux.H"
#include "fvcDiv.H"
#include "surfaceFields.H"
```

这是有意的，因为 Stage 0 还没有创建面通量，也没有计算散度。

## 4. `using namespace Foam`

代码中有：

```cpp
using namespace Foam;
```

OpenFOAM 的类型都放在 `Foam` 命名空间中。严格写法可以是：

```cpp
Foam::volScalarField T;
Foam::volVectorField U;
Foam::Info << "message";
```

有了：

```cpp
using namespace Foam;
```

以后可以直接写：

```cpp
volScalarField T;
volVectorField U;
Info << "message";
```

它只是 C++ 命名空间规则，不改变数学意义。

## 5. `main`：程序入口

```cpp
int main(int argc, char *argv[])
{
    ...
}
```

这是标准 C++ 程序入口。

其中：

- `argc` 是命令行参数个数；
- `argv` 是命令行参数内容；
- `int` 是程序退出状态。

程序最后：

```cpp
return 0;
```

表示程序正常结束。

Stage 0 的完整逻辑都放在 `main` 中。今后随着求解器变复杂，通常会把重复的数学操作
提取成函数，但当前阶段放在入口中最容易学习。

## 6. `setRootCase.H`

代码：

```cpp
#include "setRootCase.H"
```

它的作用是定位当前 case。

一个 OpenFOAM case 通常包含：

```text
case/
├── 0/
├── constant/
└── system/
```

如果当前工作目录就是 case，那么：

```bash
explicitAdvectionFoamStudent
```

程序会使用当前目录。

如果通过 `-case` 指定路径：

```bash
explicitAdvectionFoamStudent -case /some/case
```

程序会使用指定路径。

`setRootCase.H` 不是一个普通函数调用，而是一个被编译器插入当前位置的代码片段。
OpenFOAM 中大量使用这种 `.H` 初始化文件。

## 7. `createTime.H`

代码：

```cpp
#include "createTime.H"
```

它创建一个时间对象，通常叫：

```cpp
runTime
```

数学上对应独立变量：

$$t.$$

但 `runTime` 不只是一个实数。它还负责：

- 当前时间值；
- 当前时间目录名；
- 终止时间；
- 时间步长；
- 输出时间；
- 当前 case 路径；
- `system/controlDict`。

常见成员函数与数学含义：

| 代码 | 数学或工程含义 |
|---|---|
| `runTime.name()` | 当前时间目录名，例如 `0` |
| `runTime.value()` | 当前时间数值 $t^n$ |
| `runTime.endTime()` | 终止时间 |
| `runTime.setDeltaT(deltaT)` | 设置 $\Delta t$ |
| `runTime++` | $t^n\rightarrow t^{n+1}$ |
| `runTime.write()` | 输出当前时间层的场 |

Stage 0 只使用：

```cpp
runTime.name()
```

因为读取字段时需要知道当前时间目录。

例如当前时间为初始时刻：

```cpp
runTime.name() == "0"
```

于是字段路径就是：

```text
0/U
0/T
```

## 8. `createMesh.H`

代码：

```cpp
#include "createMesh.H"
```

它读取：

```text
constant/polyMesh/
```

并创建：

```cpp
fvMesh mesh;
```

`fvMesh` 是有限体积网格对象。它包含：

- 点；
- 面；
- 单元；
- 单元体积；
- 面面积向量；
- owner/neighbour 关系；
- 边界 patch；
- 网格的几何质量信息。

数学上：

$$
\Omega
\longrightarrow
\Omega_h
\longrightarrow
\{\Omega_c,V_c,\boldsymbol{S}_f\}.
$$

其中：

| 数学对象 | `fvMesh` 访问方式 |
|---|---|
| 单元数量 | `mesh.nCells()` |
| 单元体积 $V_c$ | `mesh.V()` |
| 面面积向量 | `mesh.Sf()` |
| 边界 patch | `mesh.boundary()` |

`createMesh.H` 只负责把网格读入内存，不计算对流方程。

## 9. `controlDict`

代码：

```cpp
const dictionary& controlDict = runTime.controlDict();
```

读取的文件是：

```text
system/controlDict
```

当前关键内容：

```foam
application     explicitAdvectionFoamStudent;
startTime       0;
endTime         1;
deltaT          0.01;
maxCo           0.2;
velocityField   U;
advectedField   T;
```

这些项目的含义：

| 条目 | 含义 |
|---|---|
| `application` | 当前 case 使用的求解器 |
| `startTime` | 初始时间 |
| `endTime` | 终止时间 |
| `deltaT` | 默认时间步 |
| `maxCo` | 目标 CFL 数 |
| `velocityField` | 速度字段名 |
| `advectedField` | 被输运字段名 |

Stage 0 中真正读取的是：

```cpp
velocityName
advectedName
```

Stage 2 才会读取 `maxCo`，Stage 5 才会使用时间循环。

## 10. `word` 和 `lookupOrDefault`

代码：

```cpp
const word velocityName
(
    controlDict.lookupOrDefault<word>("velocityField", "U")
);
```

`word` 是 OpenFOAM 的字符串类型。

这句话相当于：

```cpp
std::string velocityName = ...
```

但它遵守 OpenFOAM 的字段命名规则。

`lookupOrDefault` 的逻辑是：

```text
读取关键字 velocityField；
如果存在，就使用用户给的值；
如果不存在，就使用 U。
```

因此求解器不必把字段名永远写死为 `U` 和 `T`。

数学对象和字段名字的关系是：

```text
velocityName = "U"
    -> 读取文件 0/U
        -> 创建 volVectorField U
```

```text
advectedName = "T"
    -> 读取文件 0/T
        -> 创建 volScalarField T
```

## 11. `Info` 输出

代码：

```cpp
Info<< "Student checkpoint: reading fields" << nl
    << "  case       = " << runTime.caseName() << nl
    << "  time       = " << runTime.name() << nl
    << "  cells      = " << mesh.nCells() << nl
    << "  velocity   = " << velocityName << nl
    << "  advected   = " << advectedName << nl
    << endl;
```

`Info` 可以理解为 OpenFOAM 版本的标准输出流，功能类似：

```cpp
std::cout
```

但它和 OpenFOAM 的日志、并行运行机制配合得更好。

输出内容的目的不是求解，而是确认：

- 当前使用了哪个 case；
- 当前时间是多少；
- 网格有多少个单元；
- 速度字段叫什么；
- 标量字段叫什么。

这是一种很重要的工程习惯：在真正计算前打印关键状态。

## 12. 读取 `volVectorField U`

代码：

```cpp
volVectorField U
(
    IOobject
    (
        velocityName,
        runTime.name(),
        mesh,
        IOobject::MUST_READ,
        IOobject::NO_WRITE
    ),
    mesh
);
```

### 12.1 类型含义

```text
volVectorField
```

表示：

```text
每个 cell 一个 vector
```

数学上：

$$
\boldsymbol{u}(\boldsymbol{x},t)
\longrightarrow
\{\boldsymbol{U}_c^n\}_{c=1}^{N_c}.
$$

在本案例中：

$$
\boldsymbol{U}_c=(1,1,0).
$$

### 12.2 `IOobject`

内部的 `IOobject` 告诉字段如何与文件关联：

| 参数 | 当前意义 |
|---|---|
| `velocityName` | 字段名字 `U` |
| `runTime.name()` | 当前时间目录 `0` |
| `mesh` | 字段属于这个网格 |
| `MUST_READ` | 文件必须存在 |
| `NO_WRITE` | 不自动写出 |

所以程序会读取：

```text
0/U
```

### 12.3 `0/U` 的 `class`

```foam
class       volVectorField;
```

它必须与 C++ 的字段类型匹配：

```cpp
volVectorField U;
```

这是文件格式和内存对象之间的类型契约。

### 12.4 `0/U` 的量纲

```foam
dimensions      [0 1 -1 0 0 0 0];
```

表示：

$$
[\boldsymbol{U}]=L\,T^{-1}.
$$

也就是速度量纲。

### 12.5 `0/U` 的内部值

```foam
internalField   uniform (1 1 0);
```

表示所有控制体中的速度都是：

$$
\boldsymbol{U}_c=(1,1,0),
\qquad c=1,\ldots,400.
$$

如果改成 `nonuniform List<vector>`，则可以为每个 cell 指定不同速度。

## 13. 读取 `volScalarField T`

代码：

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

### 13.1 类型含义

```text
volScalarField
```

表示：

```text
每个 cell 一个 scalar
```

数学上：

$$
\phi(\boldsymbol{x},t)
\longrightarrow
\{T_c^n\}_{c=1}^{N_c}.
$$

本工程约定：

```text
题目中的 phi -> 代码中的 T
```

所以控制方程在代码变量名下是：

$$
\frac{\partial T}{\partial t}
+\nabla\cdot(\boldsymbol{U}T)=0.
$$

### 13.2 `T` 的读写规则

`T` 使用：

```cpp
IOobject::MUST_READ
IOobject::AUTO_WRITE
```

含义是：

- 初始时必须从 `0/T` 读取；
- 后续调用 `runTime.write()` 时，自动写出 `T`。

因为 `T` 是要求解的未知场，所以它会随着时间改变：

$$
T^n\longrightarrow T^{n+1}.
$$

### 13.3 `0/T` 的内部场

```foam
internalField   nonuniform List<scalar>
400
(
    ...
);
```

解释：

| 部分 | 含义 |
|---|---|
| `nonuniform` | 不同单元可以有不同值 |
| `List<scalar>` | 这是一个标量列表 |
| `400` | 必须提供 400 个值 |

因为网格有：

$$
20\times20\times1=400
$$

个控制体。

这些数值对应初始函数：

$$
T(x,y,0)=\sin(2\pi(x+y)).
$$

对于第 $c$ 个单元：

$$
T_c^0=\sin\left(2\pi(x_c+y_c)\right),
$$

其中 $(x_c,y_c)$ 是该单元中心。

## 14. `0/T` 的量纲

```foam
dimensions      [0 0 0 0 0 0 0];
```

表示 `T` 是无量纲量：

$$
[T]=1.
$$

这与正弦函数的输出相容：

$$
T=\sin(\cdot).
$$

正弦函数输出本身没有长度、时间或质量单位。

## 15. `boundaryField`

`U` 和 `T` 都有：

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

`boundaryField` 不只是文件格式要求，它对应 PDE 的边界条件。

### 15.1 周期边界

```text
xMin <-> xMax
yMin <-> yMax
```

数学上表示：

$$
T(0,y,t)=T(1,y,t),
$$

$$
T(x,0,t)=T(x,1,t).
$$

计算上，一个 patch 的面会和另一个 patch 的面配对。

### 15.2 `empty`

```text
zMin、zMax -> empty
```

这表示 z 方向没有真正求解的变化，当前问题是二维问题：

$$
\frac{\partial}{\partial z}=0.
$$

所以虽然网格形式上有一个 z 方向厚度，但物理和数值问题只在 x、y 平面上发生。

## 16. 为什么字段必须绑定 `mesh`

代码中的字段构造最后都传入：

```cpp
mesh
```

原因是 `U` 和 `T` 不是普通数组。

它们不仅包含数值，还必须知道：

- 这些值对应哪些单元；
- 共有多少个单元；
- 边界有多少个 patch；
- 每个 patch 有多少个面；
- 单元和面的拓扑关系是什么。

数学上，离散场不是单纯的向量：

$$
(T_1,T_2,\ldots,T_{400}).
$$

更准确地说，它是定义在离散区域上的函数：

$$
T_h:\Omega_h\rightarrow\mathbb{R}.
$$

其中：

- `mesh` 提供定义域 $\Omega_h$；
- `T` 提供函数值 $T_c$。

所以 `volScalarField T` 必须和 `fvMesh mesh` 绑定。

## 17. 量纲与范围检查

代码：

```cpp
Info<< "  U dimensions = " << U.dimensions() << nl
    << "  T dimensions = " << T.dimensions() << nl
    << "  T min       = " << min(T).value() << nl
    << "  T max       = " << max(T).value() << nl
    << endl;
```

### 17.1 `dimensions()`

```cpp
U.dimensions()
T.dimensions()
```

读取字段文件中声明的量纲。

预期：

```text
U dimensions = [0 1 -1 0 0 0 0]
T dimensions = []
```

这里 `[]` 表示所有量纲指数为 0。

### 17.2 `min(T)`

数学上：

$$\min_cT_c.$$

它检查所有单元的最小标量值。

### 17.3 `max(T)`

数学上：

$$\max_cT_c.$$

它检查所有单元的最大标量值。

对于正弦初值：

$$
-1\leq T_c\leq1.
$$

因此日志中出现：

```text
T min = -1
T max = 1
```

说明输入场基本正确。

## 18. Stage 0 的完整数据流

```text
system/controlDict
    |
    | runTime.controlDict()
    v
dictionary controlDict
    |
    | createTime.H
    v
Time runTime

system/blockMeshDict
    |
    | blockMesh
    v
constant/polyMesh/*
    |
    | createMesh.H
    v
fvMesh mesh

0/U
    |
    | IOobject + mesh
    v
volVectorField U

0/T
    |
    | IOobject + mesh
    v
volScalarField T
```

数学对应：

$$
t,\quad\Omega_h,\quad
\{\boldsymbol{U}_c^0\},\quad
\{T_c^0\}.
$$

截至 Stage 0，程序还没有创建：

```cpp
surfaceScalarField phi;
tmp<volScalarField> tResidual;
dimensionedScalar deltaTDim;
```

因此它还没有进入：

$$
F_f=\boldsymbol{U}_f\cdot\boldsymbol{S}_{cf},
$$

也还没有进入：

$$
T_c^{n+1}=T_c^n-\Delta tR_c^n.
$$

## 19. Stage 0 的判断标准

Stage 0 通过，需要同时满足：

- `wmake` 编译成功；
- 求解器能运行；
- `blockMesh` 成功；
- `checkMesh` 输出 `Mesh OK`；
- `0/U` 能读取；
- `0/T` 能读取；
- `U` 的量纲正确；
- `T` 的量纲正确；
- `T` 的范围符合初始条件；
- 程序在没有数值更新的情况下正常退出。

Stage 0 通过不代表方程已经求解，只代表：

```text
离散区域和初始数据已经准备好。
```

## 20. 复习问题

学习完本阶段后，应能回答：

1. `volVectorField U` 为什么是 cell 上的 vector？
2. `volScalarField T` 为什么是 cell 上的 scalar？
3. 为什么 `0/T` 要有 400 个值？
4. `IOobject::MUST_READ` 和 `IOobject::AUTO_WRITE` 分别是什么意思？
5. 为什么 `U` 使用 `NO_WRITE`，而 `T` 使用 `AUTO_WRITE`？
6. `mesh` 在数学上对应什么？
7. `boundaryField` 在 PDE 中对应什么？
8. 为什么 `surfaceScalarField phi` 还没有在 Stage 0 中出现？

如果这些问题能够回答清楚，就可以进入 Stage 1：

$$
F_f=\boldsymbol{U}_f\cdot\boldsymbol{S}_{cf}.
$$
