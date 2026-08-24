# 周期边界的代码体现

本文解释当前正弦波线性对流算例中，PDF 所说的：

```text
periodic boundary condition imposed in x-y directions
```

在 OpenFOAM 工程中到底由哪些文件、哪些语句体现。

当前案例目录是：

```text
/home/a776/workdocuments/上交船舶/slover/student_project/cases/01_sine_wave_quad
```

## 1. 数学上的周期边界是什么意思

本算例的二维区域是：

$$\Omega=[0,1]\times[0,1]$$

x 方向周期表示：

$$T(0,y,t)=T(1,y,t)$$

y 方向周期表示：

$$T(x,0,t)=T(x,1,t)$$

它的物理含义是：

```text
从右边界流出的标量，会从左边界同一位置进入；
从左边界流出的标量，会从右边界同一位置进入；
从上边界流出的标量，会从下边界同一位置进入；
从下边界流出的标量，会从上边界同一位置进入。
```

所以这个区域不是有墙的方腔，而是一个首尾相接的周期区域。

## 2. 顶点编号和定义域

周期边界首先依赖网格几何。当前定义域在：

```text
cases/01_sine_wave_quad/N20/system/blockMeshDict
```

顶点为：

```foam
vertices
(
    (0 0 0)      // vertex 0
    (1 0 0)      // vertex 1
    (1 1 0)      // vertex 2
    (0 1 0)      // vertex 3
    (0 0 0.1)    // vertex 4
    (1 0 0.1)    // vertex 5
    (1 1 0.1)    // vertex 6
    (0 1 0.1)    // vertex 7
);
```

这些顶点构成一个很薄的三维盒子。z 方向厚度是 `0.1`，但 z 方向只有一层单元，并且使用 `empty` 边界，所以物理上是二维问题。

可以把顶点位置想象成：

```text
z = 0.1 顶面:

4 -------- 5
|          |
|          |
7 -------- 6

z = 0 底面:

0 -------- 1
|          |
|          |
3 -------- 2
```

注意：

```text
顶点 0, 1, 4, 5 的 y 坐标都是 0
顶点 3, 2, 7, 6 的 y 坐标都是 1
顶点 0, 3, 4, 7 的 x 坐标都是 0
顶点 1, 2, 5, 6 的 x 坐标都是 1
```

这就是后面 `faces (...)` 能代表边界面的原因。

## 3. y 方向周期边界如何写出来

在 `blockMeshDict` 的 `boundary` 中有：

```foam
yMin
{
    type cyclic;
    neighbourPatch yMax;
    faces ((0 1 5 4));
}
yMax
{
    type cyclic;
    neighbourPatch yMin;
    faces ((3 7 6 2));
}
```

### 3.1 `yMin` 为什么是 `faces ((0 1 5 4))`

`faces ((0 1 5 4))` 的意思是：用顶点 `0, 1, 5, 4` 围成一个四边形面。

这四个点是：

```text
0 = (0, 0, 0)
1 = (1, 0, 0)
5 = (1, 0, 0.1)
4 = (0, 0, 0.1)
```

它们的共同特点是：

$$y=0$$

因此这个面就是区域下边界：

```text
yMin = y 最小的一侧 = y = 0
```

### 3.2 `yMax` 为什么是 `faces ((3 7 6 2))`

`faces ((3 7 6 2))` 的意思是：用顶点 `3, 7, 6, 2` 围成一个四边形面。

这四个点是：

```text
3 = (0, 1, 0)
7 = (0, 1, 0.1)
6 = (1, 1, 0.1)
2 = (1, 1, 0)
```

它们的共同特点是：

$$y=1$$

因此这个面就是区域上边界：

```text
yMax = y 最大的一侧 = y = 1
```

### 3.3 `type cyclic` 是什么

```foam
type cyclic;
```

表示这个 patch 不是墙面、不是入口、不是出口，而是周期边界。

如果它是墙面，标量不能穿过去；如果它是周期边界，OpenFOAM 会把这个边界和另一个 patch 连接起来。

### 3.4 `neighbourPatch` 是什么

```foam
yMin
{
    type cyclic;
    neighbourPatch yMax;
}
```

表示：

```text
yMin 的周期邻居是 yMax
```

也就是：

```text
y = 0 这一侧连接到 y = 1 这一侧
```

反过来：

```foam
yMax
{
    type cyclic;
    neighbourPatch yMin;
}
```

表示：

```text
yMax 的周期邻居是 yMin
```

也就是：

```text
y = 1 这一侧连接到 y = 0 这一侧
```

合起来就是数学条件：

$$T(x,0,t)=T(x,1,t)$$

## 4. x 方向周期边界如何写出来

x 方向同理：

```foam
xMin
{
    type cyclic;
    neighbourPatch xMax;
    faces ((0 4 7 3));
}
xMax
{
    type cyclic;
    neighbourPatch xMin;
    faces ((1 2 6 5));
}
```

`faces ((0 4 7 3))` 中的四个点为：

```text
0 = (0, 0, 0)
4 = (0, 0, 0.1)
7 = (0, 1, 0.1)
3 = (0, 1, 0)
```

共同特点是：

$$x=0$$

所以它是 `xMin`。

`faces ((1 2 6 5))` 中的四个点为：

```text
1 = (1, 0, 0)
2 = (1, 1, 0)
6 = (1, 1, 0.1)
5 = (1, 0, 0.1)
```

共同特点是：

$$x=1$$

所以它是 `xMax`。

`xMin` 和 `xMax` 互相作为 `neighbourPatch`，对应数学条件：

$$T(0,y,t)=T(1,y,t)$$

## 5. 为什么 `faces` 中顶点顺序不能乱写

OpenFOAM 用 `faces ((...))` 里的顶点顺序确定面的方向，也就是面面积向量的方向。

例如 `yMin` 写成：

```foam
faces ((0 1 5 4));
```

这个顺序让 `yMin` 的外法向指向负 y 方向。

`yMax` 写成：

```foam
faces ((3 7 6 2));
```

这个顺序让 `yMax` 的外法向指向正 y 方向。

周期面配对时，两个面的法向应该方向相反。这样 OpenFOAM 才能正确理解：

```text
一个 patch 上流出的通量，对应另一个 patch 上流入的通量。
```

因此 `faces ((0 1 5 4))` 不只是列出四个角点，它同时告诉 OpenFOAM：

```text
这个面在哪里；
这个面朝哪边；
这个面属于哪个 patch。
```

## 6. 为什么只有一个 `faces`，但 checkMesh 看到 20 个面

`blockMeshDict` 中的块定义是：

```foam
hex (0 1 2 3 4 5 6 7) (20 20 1) simpleGrading (1 1 1)
```

这表示整个大六面体会被切成：

```text
x 方向 20 个单元
y 方向 20 个单元
z 方向 1 个单元
```

因此 `faces ((0 1 5 4))` 先定义的是 `y=0` 的大边界面。`blockMesh` 生成网格时，会把这个大面切成 20 个小 patch face。

所以 `checkMesh` 里会看到：

```text
yMin    20 faces
yMax    20 faces
```

这不是矛盾，而是：

```text
blockMeshDict 中写大面；
blockMesh 根据网格分辨率把大面切成多个小边界面。
```

## 7. 场文件里也必须写 `cyclic`

只在 `blockMeshDict` 里写周期边界还不够。每个场文件的 `boundaryField` 也要和网格 patch 类型一致。

速度场 `U` 在：

```text
cases/01_sine_wave_quad/N20/0.orig/U
```

其中写了：

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

标量场 `T` 是由脚本生成的，脚本在：

```text
scripts/01_sine_wave_quad/create_initial_fields.py
```

生成的 `T` 文件中同样有：

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

这表示：

```text
网格说 x/y 是 cyclic；
U 场也承认 x/y 是 cyclic；
T 场也承认 x/y 是 cyclic。
```

如果三者不一致，OpenFOAM 通常会报 patch 类型不匹配，或者边界行为不符合预期。

## 8. 求解器里为什么没有手写周期逻辑

求解器里没有写：

```cpp
if (y < 0) y = 1;
if (y > 1) y = 0;
```

这是因为 OpenFOAM 的有限体积算子会根据 mesh patch 类型自动处理边界。

当前求解器先计算面通量：

```cpp
surfaceScalarField phi
(
    IOobject(...),
    fvc::flux(U)
);
```

然后计算对流散度：

```cpp
tmp<volScalarField> tResidual
(
    fvc::div(phi, T, "div(phi,T)")
);
```

这里的 `fvc::div(phi,T,...)` 会读取 `fvSchemes` 中的：

```foam
div(phi,T) Gauss upwind;
```

当某个面是 `cyclic` 边界面时，OpenFOAM 会把它的 `neighbourPatch` 上对应的面当成周期邻居来处理。也就是说，边界外侧不是空的，也不是墙，而是另一个 patch 的相应单元。

更新后，求解器调用：

```cpp
T.correctBoundaryConditions();
```

这一步会让 `T` 的边界值重新满足场文件中定义的边界条件。对于 `cyclic` patch，就是同步配对 patch 的周期关系。

## 9. checkMesh 如何验证周期配对

运行 `Allrun` 时会执行：

```bash
checkMesh
```

在当前日志：

```text
cases/01_sine_wave_quad/N20/log.checkMesh
```

可以看到类似信息：

```text
Boundary definition OK.
xMin                20       42       ok
xMax                20       42       ok
yMin                20       42       ok
yMax                20       42       ok
Coupled point location match (average 0) OK.
Mesh OK.
```

其中 `Coupled point location match` 很关键。它说明 OpenFOAM 检查了这些耦合边界，也就是 cyclic patch 的点位匹配关系。

如果周期面数量不一致、几何位置无法配对、点顺序不合理，就可能在这一步报错。

## 10. 初值为什么适合周期边界

当前初值是：

$$T(x,y,0)=\sin(2\pi(x+y))$$

检查 x 方向周期性：

$$T(x+1,y,0)=\sin(2\pi(x+1+y))$$

$$T(x+1,y,0)=\sin(2\pi(x+y)+2\pi)=\sin(2\pi(x+y))$$

所以：

$$T(0,y,0)=T(1,y,0)$$

检查 y 方向周期性：

$$T(x,y+1,0)=\sin(2\pi(x+y+1))$$

$$T(x,y+1,0)=\sin(2\pi(x+y)+2\pi)=\sin(2\pi(x+y))$$

所以：

$$T(x,0,0)=T(x,1,0)$$

这就是为什么正弦波平移问题适合用周期边界。标量平移穿过边界时不会产生跳跃，解析解也能在周期区域里自然延续。

## 11. 一句话总结

当前案例的周期边界由三层共同体现：

```text
1. blockMeshDict:
   xMin <-> xMax, yMin <-> yMax，通过 type cyclic 和 neighbourPatch 配对。

2. U 和 T 的 boundaryField:
   每个场都声明 x/y patch 是 cyclic。

3. 求解器:
   fvc::flux、fvc::div 和 correctBoundaryConditions 根据 mesh/field 的 cyclic 信息自动跨周期边界计算。
```

所以 `yMin` 和 `yMax` 那段代码的核心意思是：

```text
把 y = 0 的边界面和 y = 1 的边界面接起来，让它们成为同一个周期方向的两端。
```
