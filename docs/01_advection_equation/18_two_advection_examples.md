# 第一题两个案例的函数、速度和边界对应关系

第一题的两个案例解的是同一个守恒型线性对流方程：

$$
\frac{\partial \phi}{\partial t}
+
\nabla\cdot(\boldsymbol{u}\phi)=0.
$$

区别主要在于：

```text
速度场不同
初始标量场不同
边界条件不同
最终检查目标不同
```

求解器主体仍然可以复用：

```text
fvc::flux(U)
fvc::div(phi,T)
T = T - deltaT*residual
```

## 1. 案例一：正弦波平移

题目参数：

$$
\boldsymbol{u}=(1,1),
\qquad
\Omega=[0,1]^2,
$$

$$
\phi(x,y,0)=\sin(2\pi(x+y)).
$$

解析解：

$$
\phi(x,y,t)
=
\sin\left(2\pi(x+y-2t)\right).
$$

题目要求：

```text
x 和 y 方向周期边界
CFL=0.2
t=1.0 时的 L1 误差
网格收敛阶
```

### 1.1 JSON

四边形 + 一阶迎风：

```json
{
  "caseName": "01_sine_wave_quad",
  "problem": "sine_wave_advection",
  "meshType": "quad",
  "schemeName": "upwind",
  "divScheme": "Gauss upwind",
  "velocity": [1.0, 1.0, 0.0],
  "endTime": 1.0,
  "maxCo": 0.2
}
```

四边形 + `linearUpwind`：

```json
{
  "caseName": "02_sine_wave_quad_linearUpwind",
  "problem": "sine_wave_advection",
  "meshType": "quad",
  "schemeName": "linearUpwind",
  "divScheme": "Gauss linearUpwind grad(T)",
  "gradTScheme": "Gauss linear",
  "velocity": [1.0, 1.0, 0.0],
  "endTime": 1.0,
  "maxCo": 0.2
}
```

### 1.2 速度场在哪里

JSON：

```json
"velocity": [1.0, 1.0, 0.0]
```

准备脚本：

```text
scripts/common/foam_case.py
    -> _patch_velocity_field()
```

最终写入：

```text
cases/01_advection_equation/01_sine_wave_quad/N20/0.orig/U
```

内容：

```foam
internalField uniform (1 1 0);
```

求解器读取：

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

这段代码中每个括号里的内容都对应一个具体对象：

```cpp
volVectorField U( ..., mesh );
```

最外层构造函数的两个参数是：

| 参数 | 实际内容 | 作用 |
|---|---|---|
| 第一个参数 | `IOobject(...)` | 描述字段名字、时间目录和读写方式 |
| 第二个参数 | `mesh` | 指定这个体场属于哪个有限体积网格 |

内部 `IOobject(...)` 的参数依次是：

```cpp
IOobject
(
    velocityName,             // 字段名字，默认是 "U"
    runTime.name(),           // 当前时间目录，初始时通常是 "0"
    mesh,                     // 该字段注册到哪个网格对象
    IOobject::MUST_READ,      // 必须从磁盘读取
    IOobject::NO_WRITE        // 求解器不自动写 U
)
```

字段名由下面的调用得到：

```cpp
const word velocityName
(
    controlDict.lookupOrDefault<word>("velocityField", "U")
);
```

`lookupOrDefault<word>(...)` 括号内的两个参数是：

```cpp
controlDict.lookupOrDefault<word>
(
    "velocityField",  // 要查找的 controlDict 关键字
    "U"               // 关键字不存在时使用的默认字段名
)
```

因此完整链条是：

```text
system/controlDict:
    velocityField U;

lookupOrDefault<word>("velocityField", "U")
    -> velocityName = "U"

IOobject("U", runTime.name(), mesh, MUST_READ, NO_WRITE)
    -> 读取 0/U

volVectorField(..., mesh)
    -> 得到 cell-centered 速度场 U
```

第一案例中，`0.orig/U` 的内容为：

```foam
internalField uniform (1 1 0);
```

所以：

```text
U 的 x 分量 = 1
U 的 y 分量 = 1
U 的 z 分量 = 0
```

### 1.3 初始函数在哪里

数学函数：

$$
\phi_0(x,y)=\sin(2\pi(x+y)).
$$

代码：

```text
scripts/common/advection_sine.py
    -> write_case_initial_field()
    -> write_case_initial_field_from_centres()
```

四边形初值函数的完整签名是：

```python
def write_case_initial_field(
    case: Path,
    nx: int,
    ny: int,
    velocity: tuple[float, float, float] = (1.0, 1.0, 0.0),
) -> Path:
```

括号内参数含义：

| 参数 | 类型 | 作用 |
|---|---|---|
| `case` | `Path` | 当前案例目录，例如 `cases/01_advection_equation/01_sine_wave_quad/N20` |
| `nx` | `int` | x 方向单元数 |
| `ny` | `int` | y 方向单元数 |
| `velocity` | 三元组 | 速度 `(1.0,1.0,0.0)`，用于保持统一接口 |

它内部调用：

```python
write_initial_field(case, nx, ny, velocity)
```

这四个位置参数依次传入：

```text
case      -> 输出目录
nx, ny    -> 生成多少个 cell 值
velocity  -> 当前案例速度
```

底层在每个规则网格中心计算：

```python
x = (i + 0.5) / nx
y = (j + 0.5) / ny
value = math.sin(2.0 * math.pi * (x + y))
```

因此：

$$
T(x,y,0)=\sin(2\pi(x+y)).
$$

三角形初值函数的完整签名是：

```python
def write_case_initial_field_from_centres(
    case: Path,
    centres: list[tuple[float, float, float]],
    velocity: tuple[float, float, float] = (1.0, 1.0, 0.0),
) -> Path:
```

调用：

```python
values = exact_values_at_centres(centres, 0.0, velocity)
```

括号中的参数分别是：

```python
exact_values_at_centres(
    centres,   # 从 constant/C 读取的真实 cell centre
    0.0,       # 初始时间 t=0
    velocity,  # (1,1,0)
)
```

随后调用：

```python
write_scalar_field(
    case / "0.orig" / "T",
    values,
    {
        "xMin": "cyclic",
        "xMax": "cyclic",
        "yMin": "cyclic",
        "yMax": "cyclic",
        "zMin": "empty",
        "zMax": "empty",
    },
    object_name="T",
    location="0",
)
```

参数含义：

| 参数 | 作用 |
|---|---|
| `case / "0.orig" / "T"` | 输出文件路径 |
| `values` | 按 cell 顺序排列的标量初值 |
| 第三个字典 | 每个 patch 的边界类型 |
| `object_name="T"` | OpenFOAM 字段名 |
| `location="0"` | 字段属于 0 时间目录 |

所以三角形初值的实际流程是：

```text
constant/C
    -> centres
    -> exact_values_at_centres(centres, 0.0, velocity)
    -> values
    -> write_scalar_field(...)
    -> 0.orig/T
```

四边形使用规则 cell centre：

$$
x_i=\frac{i+1/2}{N},
\qquad
y_j=\frac{j+1/2}{N}.
$$

三角形使用真实 cell centre：

```text
constant/C
    -> read_cell_geometry()
    -> exact_values_at_centres()
    -> 写入 0.orig/T
```

最终文件：

```text
cases/<caseName>/Nxx/0.orig/T
```

### 1.4 周期边界在哪里

四边形：

```text
cases/01_advection_equation/01_sine_wave_quad/N20/system/blockMeshDict
cases/01_advection_equation/01_sine_wave_quad/N20/0.orig/U
cases/01_advection_equation/01_sine_wave_quad/N20/0.orig/T
```

四个方向：

```foam
xMin <-> xMax cyclic
yMin <-> yMax cyclic
zMin, zMax empty
```

三角形：

```text
scripts/common/gmsh_tri_mesh.py
scripts/common/foam_case.py::_write_create_patch_dict()
```

三角形的 `createPatchDict` 把 Gmsh 的 source patch 变成：

```text
xMin/xMax cyclic
yMin/yMax cyclic
zMin/zMax empty
```

### 1.5 误差和收敛阶在哪里

解析解：

```text
scripts/common/advection_tools.py
    -> exact_values()
    -> exact_values_at_centres()
```

误差：

$$
L_1
=
\frac{
\sum_cV_c|T_c-T_c^{\mathrm{exact}}|
}{
\sum_cV_c|T_c^{\mathrm{exact}}|
}.
$$

代码：

```text
scripts/common/metrics.py
scripts/common/postprocess_case.py
```

收敛阶：

$$
p
=
\frac{\log(E_N/E_{2N})}{\log(2)}.
$$

代码：

```text
scripts/common/study_analysis.py
    -> observed_order()
    -> analyse()
```

## 2. 案例二：复杂轮廓刚体旋转

题面速度场：

$$
\boldsymbol{u}(x,y)
=(0.5-y,\;x-0.5).
$$

这表示绕：

$$
(x_c,y_c)=(0.5,0.5)
$$

以单位角速度逆时针旋转。

终止时间：

$$
t=2\pi.
$$

### 2.1 旋转速度的代码

当前实现：

```text
scripts/common/advection_rotation.py
    -> solid_rotation_velocity()
```

完整函数签名：

```python
def solid_rotation_velocity(
    x: float,
    y: float,
    center: tuple[float, float] = (0.5, 0.5),
    angular_velocity: float = 1.0,
) -> tuple[float, float, float]:
```

括号内参数含义：

| 参数 | 数学含义 | 当前值 |
|---|---|---|
| `x` | 当前 cell centre 的 x 坐标 | 例如 `0.505` |
| `y` | 当前 cell centre 的 y 坐标 | 例如 `0.745` |
| `center` | 旋转中心 $(x_c,y_c)$ | `(0.5,0.5)` |
| `angular_velocity` | 角速度 $\omega$ | `1.0` |

函数内部：

```python
xc, yc = center
```

把：

```python
center=(0.5, 0.5)
```

拆成：

```text
xc=0.5
yc=0.5
```

返回表达式：

```python
return (
    angular_velocity * (yc - y),
    angular_velocity * (x - xc),
    0.0,
)
```

对应数学公式：

$$
u=\omega(y_c-y),
\qquad
v=\omega(x-x_c),
\qquad
w=0.
$$

实际调用：

```python
solid_rotation_velocity(
    x,
    y,
    center=(0.5, 0.5),
    angular_velocity=1.0,
)
```

返回：

```python
(0.5-y, x-0.5, 0.0)
```

代码公式：

```python
(
    angular_velocity * (yc - y),
    angular_velocity * (x - xc),
    0.0,
)
```

当：

```text
center=(0.5,0.5)
angularVelocity=1.0
```

得到：

$$
(u,v)=(0.5-y,\;x-0.5).
$$

### 2.2 旋转速度写入哪里

JSON：

```json
"velocity": [0.0, 0.0, 0.0],
"velocityModel": {
  "type": "solidRotation",
  "center": [0.5, 0.5],
  "angularVelocity": 1.0
}
```

这里 `velocity` 只是保留字段，旋转案例真正使用的是 `velocityModel`。

准备流程：

```text
foam_case.py::_patch_velocity_field()
    -> advection_rotation.py::write_case_velocity_field()
    -> foam_fields.py::write_vector_field()
    -> 0.orig/U
```

具体的写入函数签名：

```python
def write_case_velocity_field(
    case: Path,
    nx: int,
    ny: int,
    domain: tuple[float, float, float, float] = (0.0, 1.0, 0.0, 1.0),
    velocity_config: dict[str, Any] | None = None,
) -> Path:
```

参数含义：

| 参数 | 内容 | 用途 |
|---|---|---|
| `case` | 当前案例目录 | 写出 `0.orig/U` |
| `nx` | x 方向 cell 数 | 生成 cell centre |
| `ny` | y 方向 cell 数 | 生成 cell centre |
| `domain` | `(xmin,xmax,ymin,ymax)` | 定义计算区域 |
| `velocity_config` | JSON 的 `velocityModel` | 读取旋转中心和角速度 |

函数先执行：

```python
velocity_config = velocity_config or {}
```

意思是：没有配置时使用空字典。

然后：

```python
center_value = velocity_config.get("center", [0.5, 0.5])
angular_velocity = float(
    velocity_config.get("angularVelocity", 1.0)
)
```

例如 JSON：

```json
"velocityModel": {
  "type": "solidRotation",
  "center": [0.5, 0.5],
  "angularVelocity": 1.0
}
```

会变成：

```python
center = (0.5, 0.5)
angular_velocity = 1.0
```

接下来：

```python
values = [
    solid_rotation_velocity(
        x,
        y,
        center,
        angular_velocity,
    )
    for x, y in _cell_centres(nx, ny, domain)
]
```

也就是说，每个 cell centre 都调用一次：

```python
solid_rotation_velocity(x, y, center, angular_velocity)
```

返回一个：

```text
(u_x,u_y,u_z)
```

最后写出：

```python
write_vector_field(
    case / "0.orig" / "U",
    values,
    VELOCITY_BOUNDARY,
)
```

三个参数分别是：

```text
输出文件
所有 cell 的速度列表
速度场边界条件
```

最终 `U` 是非均匀 cell-centered 场：

```foam
internalField nonuniform List<vector>
...
```

求解器仍然不需要知道“这是旋转”。它只读取 `U`，然后统一执行：

```cpp
fvc::flux(U)
```

这就是案例二复用求解器主体的关键。

### 2.3 开槽圆盘

圆心：

$$
(x_D,y_D)=(0.5,0.75),
\qquad r_0=0.15.
$$

代码：

```text
advection_rotation.py::slotted_disk()
```

完整函数签名：

```python
def slotted_disk(
    x: float,
    y: float,
    radius: float = 0.15,
    center: tuple[float, float] = (0.5, 0.75),
    slot_half_width: float = 0.025,
    slot_top_y: float = 0.85,
) -> float:
```

参数含义：

| 参数 | 含义 |
|---|---|
| `x`, `y` | 当前评价点坐标 |
| `radius` | 圆盘半径 |
| `center` | 圆盘中心 |
| `slot_half_width` | 槽宽的一半 |
| `slot_top_y` | 槽的顶部高度 |

函数内部：

```python
xc, yc = center
distance = math.hypot(x - xc, y - yc)
outside_slot = abs(x - xc) >= slot_half_width or y >= slot_top_y
```

括号对应：

```python
math.hypot(x - xc, y - yc)
```

即：

$$
\sqrt{(x-x_c)^2+(y-y_c)^2}.
$$

```python
abs(x - xc) >= slot_half_width
```

判断点是否在槽的水平范围之外。

```python
y >= slot_top_y
```

判断点是否在槽顶部以上。

最后：

```python
return 1.0 if distance <= radius and outside_slot else 0.0
```

等价于：

```python
if distance <= radius and outside_slot:
    return 1.0
else:
    return 0.0
```

对应：

$$
r_D=\sqrt{(x-0.5)^2+(y-0.75)^2}.
$$

当：

$$
r_D\leq0.15
\quad\text{且}\quad
\left(|x-0.5|\geq0.025\ \text{或}\ y\geq0.85\right),
$$

取：

$$
\phi_D=1.
$$

否则：

$$
\phi_D=0.
$$

JSON：

```json
"initialProfile": {
  "type": "slottedDiskConeCosineHump",
  "radius": 0.15,
  "diskCenter": [0.5, 0.75],
  "slotHalfWidth": 0.025,
  "slotTopY": 0.85
}
```

JSON 到函数参数的对应关系：

```text
radius        -> radius
diskCenter    -> center
slotHalfWidth -> slot_half_width
slotTopY      -> slot_top_y
```

### 2.4 圆锥

圆心：

$$
(x_C,y_C)=(0.5,0.25).
$$

代码：

```text
advection_rotation.py::cone()
```

完整函数签名：

```python
def cone(
    x: float,
    y: float,
    radius: float = 0.15,
    center: tuple[float, float] = (0.5, 0.25),
) -> float:
```

实际调用：

```python
cone(
    x,
    y,
    params["radius"],
    params["cone_center"],
)
```

括号中的参数依次是：

```text
x                     当前点 x
y                     当前点 y
params["radius"]      半径
params["cone_center"] 圆锥中心
```

函数内部：

```python
xc, yc = center
distance = math.hypot(x - xc, y - yc)
if distance > radius:
    return 0.0
return 1.0 - distance / radius
```

含义是：

```text
距离超过半径 -> 返回 0
距离在半径内 -> 返回 1-distance/radius
```

公式：

$$
r_C=\sqrt{(x-0.5)^2+(y-0.25)^2},
$$

$$
\phi_C=
\begin{cases}
1-r_C/0.15,&r_C\leq0.15,\\
0,&r_C>0.15.
\end{cases}
$$

JSON：

```json
"coneCenter": [0.5, 0.25]
```

### 2.5 光滑凸峰

圆心：

$$
(x_H,y_H)=(0.25,0.5).
$$

代码：

```text
advection_rotation.py::cosine_hump()
```

完整函数签名：

```python
def cosine_hump(
    x: float,
    y: float,
    radius: float = 0.15,
    center: tuple[float, float] = (0.25, 0.5),
) -> float:
```

实际调用：

```python
cosine_hump(
    x,
    y,
    params["radius"],
    params["hump_center"],
)
```

括号中的参数依次是：

```text
x                     当前点 x
y                     当前点 y
params["radius"]      凸峰半径
params["hump_center"] 凸峰中心
```

函数内部：

```python
distance = math.hypot(x - xc, y - yc)
if distance > radius:
    return 0.0
return 0.25 * (
    1.0 + math.cos(math.pi * distance / radius)
)
```

其中：

```python
math.cos(math.pi * distance / radius)
```

对应：

$$
\cos\left(\pi\frac{r_H}{r_0}\right).
$$

公式：

$$
r_H=\sqrt{(x-0.25)^2+(y-0.5)^2},
$$

$$
\phi_H=
\begin{cases}
\dfrac14
\left[
1+\cos\left(\pi r_H/0.15\right)
\right],
&r_H\leq0.15,\\
0,&r_H>0.15.
\end{cases}
$$

JSON：

```json
"humpCenter": [0.25, 0.5]
```

### 2.6 第二案例的总初值

代码：

```text
advection_rotation.py::solid_rotation_profile()
```

完整函数签名：

```python
def solid_rotation_profile(
    x: float,
    y: float,
    config: dict[str, Any] | None = None,
) -> float:
```

参数含义：

| 参数 | 含义 |
|---|---|
| `x` | 当前 cell centre 的 x 坐标 |
| `y` | 当前 cell centre 的 y 坐标 |
| `config` | JSON 中的 `initialProfile` 字典 |

函数第一步：

```python
params = _profile_parameters(config or {})
```

`config or {}` 的含义是：如果没有传入 `config`，就用空字典。

配置解析函数的完整签名：

```python
def _profile_parameters(
    config: dict[str, Any],
) -> dict[str, float]:
```

这个函数的括号里只有一个参数：

```text
config -> JSON 中的 initialProfile 对象
```

它内部还有一个辅助函数：

```python
def point(
    name: str,
    default: tuple[float, float],
) -> tuple[float, float]:
```

`point(...)` 的两个参数分别是：

```text
name    -> JSON 字段名，例如 "diskCenter"
default -> 字段缺失时使用的默认二维坐标
```

例如：

```python
point("diskCenter", (0.5, 0.75))
```

执行过程是：

```python
value = config.get(name, default)
```

括号内：

```text
name    -> 要查找的键
default -> 找不到键时的备用值
```

然后：

```python
return (float(value[0]), float(value[1]))
```

把 JSON 列表：

```json
[0.5, 0.75]
```

转换成 Python 元组：

```python
(0.5, 0.75)
```

`_profile_parameters(...)` 读取 JSON，并把 JSON 的驼峰命名转换成 Python
内部名称：

```python
{
    "radius": float(config.get("radius", 0.15)),
    "disk_center": point("diskCenter", (0.5, 0.75)),
    "cone_center": point("coneCenter", (0.5, 0.25)),
    "hump_center": point("humpCenter", (0.25, 0.5)),
    "slot_half_width": float(config.get("slotHalfWidth", 0.025)),
    "slot_top_y": float(config.get("slotTopY", 0.85)),
}
```

例如：

```python
point("diskCenter", (0.5, 0.75))
```

括号中的两个参数是：

```text
"diskCenter"    JSON 中的字段名
(0.5,0.75)      字段不存在时的默认圆心
```

最终总函数调用：

```python
return (
    slotted_disk(
        x,
        y,
        params["radius"],
        params["disk_center"],
        params["slot_half_width"],
        params["slot_top_y"],
    )
    + cone(
        x,
        y,
        params["radius"],
        params["cone_center"],
    )
    + cosine_hump(
        x,
        y,
        params["radius"],
        params["hump_center"],
    )
)
```

这段代码和数学公式逐项对应：

```text
slotted_disk(...) -> phi_D(x,y)
cone(...)         -> phi_C(x,y)
cosine_hump(...)  -> phi_H(x,y)
三个返回值相加    -> phi_0(x,y)
```

因此：

$$
\phi_0(x,y)
=
\phi_D(x,y)+\phi_C(x,y)+\phi_H(x,y).
$$

四边形网格中，脚本按结构化 cell centre 写入：

```text
cases/01_advection_equation/04_solid_rotation_quad_upwind/N100/0.orig/T
```

完整的四边形初值写入函数：

```python
def write_case_initial_field(
    case: Path,
    nx: int,
    ny: int,
    domain: tuple[float, float, float, float] = (0.0, 1.0, 0.0, 1.0),
    profile_config: dict[str, Any] | None = None,
) -> Path:
```

它首先调用：

```python
_cell_centres(nx, ny, domain)
```

`_cell_centres(...)` 的完整签名：

```python
def _cell_centres(
    nx: int,
    ny: int,
    domain: tuple[float, float, float, float],
) -> list[tuple[float, float]]:
```

括号参数：

```text
nx, ny -> x/y 方向单元数量
domain -> (xmin, xmax, ymin, ymax)
```

内部计算：

```python
xmin, xmax, ymin, ymax = domain
dx = (xmax - xmin) / nx
dy = (ymax - ymin) / ny
```

第 `(i,j)` 个 cell centre：

```python
(
    xmin + (i + 0.5) * dx,
    ymin + (j + 0.5) * dy,
)
```

然后对每个中心调用：

```python
solid_rotation_profile(x, y, profile_config)
```

最后调用：

```python
write_scalar_field(
    case / "0.orig" / "T",
    values,
    SCALAR_BOUNDARY,
)
```

三个位置参数分别是：

```text
输出文件路径
所有 cell 的 phi_0 值
标量场边界条件
```

三角形网格若要支持第二案例，必须改成：

```text
真实三角形 cell centre
    -> solid_rotation_profile(xc,yc)
    -> 0.orig/T
```

不能继续假设规则的 $N\times N$ cell 排列。

### 2.7 第二案例边界

三个轮廓都位于区域内部，题目不是周期正弦波边界。

当前四边形旋转分支使用：

```foam
xMin/xMax/yMin/yMax  type patch;
zMin/zMax            type empty;
```

`T` 的外边界使用：

```foam
xMin { type fixedValue; value uniform 0; }
xMax { type fixedValue; value uniform 0; }
yMin { type fixedValue; value uniform 0; }
yMax { type fixedValue; value uniform 0; }
zMin { type empty; }
zMax { type empty; }
```

这是因为初始轮廓紧支撑在区域内部，外边界上的标量可以取 0。

### 2.8 第二案例最终图和指标

题面主要要求：

```text
t=2*pi 时的标量场等值线图
```

当前后处理：

```text
scripts/common/postprocess_case.py
    -> _postprocess_quad_solid_rotation_case()
    -> plot_rotation_final_contour()
```

输出：

```text
figures/01_advection_equation/cases/04_solid_rotation_quad_upwind/N100/contour_final.png
```

同时记录：

```text
cycleL1AgainstInitial
initialMass
finalMass
normalizedMassError
minCo
maxCo
```

注意：案例二的一圈后 `L1` 只是“与初始场的差异检查”，不是第一案例那种
有解析平移解的网格收敛阶主指标。

## 3. 两个案例哪些东西复用

| 内容 | 案例一 | 案例二 | 是否复用 |
|---|---|---|---|
| PDE | 线性对流 | 线性对流 | 是 |
| `fvc::flux(U)` | 使用 | 使用 | 是 |
| `fvc::div(phi,T)` | 使用 | 使用 | 是 |
| CFL | `maxCo=0.2` | `maxCo=0.2` | 是 |
| 时间推进 | 前向 Euler | 前向 Euler | 是 |
| `U` | 常量 | 空间变化 | 字段接口复用 |
| `T` 初值 | 正弦函数 | 三个轮廓之和 | 生成函数不同 |
| 边界 | 周期 | 外边界零值 | 配置分支不同 |
| 后处理 | 解析误差/收敛阶 | 等值线/一圈差异 | 后处理分支不同 |

所以第二案例不是重新写一个求解器，而是扩展：

```text
配置模型
速度场生成
初值场生成
边界生成
后处理
```

## 4. 当前代码边界

当前已经实现：

```text
案例二 + 四边形 + upwind
```

当前没有实现：

```text
案例二 + 三角形
案例二 + linearUpwind
案例二 + 多分辨率收敛分析
```

后面三类扩展都需要先检查：

1. 旋转速度是否按真实 cell centre 写入；
2. 三角形外边界是否不再错误地使用周期 patch；
3. `T` 的三角形后处理是否使用真实 `C` 和 `Vc`；
4. `linearUpwind` 的梯度格式是否能在该边界和网格上稳定运行；
5. 结果图是否按非结构网格几何绘制。
