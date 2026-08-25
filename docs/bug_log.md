# CFD 项目 Bug 与解决记录

本文档集中记录本项目开发、编译、运行和后处理过程中发现的问题。

记录每一个 Bug 的目的不是简单记住“某一行代码错了”，而是建立下面这条
可复用的排查链：

```text
现象
    -> 日志中的真实错误
    -> 所属层次
    -> 根本原因
    -> 修复位置
    -> 验证方式
    -> 后续预防措施
```

项目根目录：

```text
/home/a776/workdocuments/上交船舶/slover/student_project
```

## 0. 如何使用本日志

以后遇到新问题，按下面格式追加：

```text
## Bug N：简短标题

状态：已修复 / 待修复 / 已知限制
日期：YYYY-MM-DD
影响范围：quad / tri / 所有案例

### 1. 现象

### 2. 日志证据

### 3. 根因

### 4. 修复

### 5. 验证

### 6. 经验
```

注意区分三类情况：

```text
真正 Bug
    -> 程序或配置与预期设计不一致

环境警告
    -> 出现提示，但没有导致当前任务失败

数值方法特性
    -> 例如一阶迎风的数值耗散，不一定是程序错误
```

## Bug 1：OpenFOAM 14 中错误使用 `runTime.timeName()`

状态：已修复  
影响范围：求解器编译

### 1. 现象

初稿求解器在创建字段或输出时间目录时使用了：

```cpp
runTime.timeName()
```

在当前 OpenFOAM 14 环境中，这个调用形式与实际接口不匹配，导致编译失败
或接口解析错误。

### 2. 根因

OpenFOAM 的 `Time` 类中，当前运行对象取得时间目录名使用：

```cpp
runTime.name()
```

而 `Time::timeName(...)` 主要是把数值时间格式化成目录名的函数，不应按照
无参数实例方法随意调用。

### 3. 修复

将求解器中的相关调用改为：

```cpp
runTime.name()
```

并查阅本机源码：

```text
/opt/openfoam14/src/OpenFOAM/db/Time/Time.H
/opt/openfoam14/src/OpenFOAM/db/Time/Time.C
```

### 4. 验证

重新执行：

```bash
source /opt/openfoam14/etc/bashrc
sh scripts/build_student_solver.sh
```

求解器成功生成：

```text
build/bin/explicitAdvectionFoamStudent
```

### 5. 经验

OpenFOAM 不同版本的 API 不能只凭记忆使用。遇到成员函数不确定时，先查：

```bash
rg -n "timeName|name\\(" /opt/openfoam14/src/OpenFOAM/db/Time
```

## Bug 2：`linearUpwind` 生成的 `fvSchemes` 大括号层级错误

状态：已修复  
影响范围：三角形或四边形 `linearUpwind` 案例

### 1. 现象

第一次运行三角形 `linearUpwind` 案例时，Python 外层显示：

```text
subprocess.CalledProcessError
```

这个异常只是说明 `Allrun` 返回了非零状态，真正错误来自 OpenFOAM：

```text
FOAM FATAL IO ERROR:
keyword divSchemes is undefined in dictionary .../system/fvSchemes
```

### 2. 错误文件

生成的 `system/fvSchemes` 曾经类似：

```foam
gradSchemes
{
    default         Gauss linear;
}
    grad(T)         Gauss linear;
}
```

多出的关闭大括号使 `divSchemes` 不再处于顶层字典。

### 3. 根因

脚本使用简单文本替换向 `gradSchemes` 插入：

```foam
grad(T)         Gauss linear;
```

但没有正确识别 `gradSchemes` 的完整大括号范围，导致插入位置和关闭大括号
发生错位。

### 4. 修复

修复文件：

```text
scripts/common/foam_case.py
```

修复函数：

```python
_patch_fv_schemes()
```

新逻辑会：

```text
定位 gradSchemes 的开始
    -> 逐行计算大括号深度
    -> 找到唯一对应的结束大括号
    -> 将 grad(T) 插入结束大括号之前
```

### 5. 验证

重新准备案例：

```bash
python3 scripts/run_study.py \
    --config scripts/configs/03_sine_wave_tri_linearUpwind.json \
    --resolutions 10,20,40,80 \
    --prepare-only \
    --overwrite
```

检查：

```bash
sed -n '20,60p' cases/03_sine_wave_tri_linearUpwind/N10/system/fvSchemes
```

必须看到顶层结构：

```foam
gradSchemes
{
    ...
}

divSchemes
{
    ...
}
```

### 6. 经验

OpenFOAM 字典是有层级的结构化文本。修改字典时，不能只依赖脆弱的字符串
替换；至少要跟踪大括号层级。

## Bug 3：三角形案例的初值不能按规则 `N x N` 生成

状态：已修复  
影响范围：三角形网格案例

### 1. 现象

四边形案例可以用：

```text
x_i = (i+0.5)/N
y_j = (j+0.5)/N
```

直接生成 cell centre。

但三角形案例的 cell 编号、几何中心和实际单元数量不能简单等同于规则
`N x N` 数组。如果继续使用四边形逻辑，可能出现：

```text
T 数量与网格 cell 数不一致
初值写入错误 cell
后处理 reshape 失败
误差坐标错误
```

### 2. 根因

三角形网格必须先生成真实网格，再从 OpenFOAM 网格中读取：

```text
constant/C
    -> cell centre

constant/Vc
    -> cell volume
```

### 3. 修复

三角形流程改为：

```text
Gmsh 生成网格
    -> gmshToFoam
    -> foamPostProcess 写 C 和 Vc
    -> 读取真实 cell centre
    -> 在每个 (xc,yc) 上计算初值
    -> 写入 0.orig/T
```

正弦波使用：

```python
write_case_initial_field_from_centres(...)
```

旋转案例使用：

```python
write_rotation_initial_field_from_centres(...)
write_rotation_velocity_field_from_centres(...)
```

### 4. 验证

三角形后处理使用：

```text
mesh/mesh_geometry.json
constant/C
constant/Vc
```

而不是把场值强行 reshape 成 `N x N`。

### 5. 经验

`N` 在三角形案例中表示基准边划分数，不等于最终 OpenFOAM cell 总数。
所有初值、速度和体积加权误差都必须服从真实网格的 cell 编号。

## Bug 4：旋转三角形案例错误使用周期边界

状态：已修复  
影响范围：固体旋转三角形案例

### 1. 现象

正弦波案例需要：

```text
x/y 方向周期边界
```

但固体旋转案例的标量轮廓位于区域内部，外边界应使用零值边界。若直接
复用正弦波三角形的 `createPatchDict`，会把旋转案例外边界错误地生成为
周期 patch。

### 2. 根因

三角形网格生成本身可以复用，但边界条件属于物理问题配置，不能只由
`meshType=tri` 决定。

### 3. 修复

根据 JSON 中的：

```json
"boundaryCondition": "zeroScalarAtOuterBoundary"
```

生成：

```foam
xMin { type patch; }
xMax { type patch; }
yMin { type patch; }
yMax { type patch; }
zMin { type empty; }
zMax { type empty; }
```

旋转案例的 `T` 使用：

```foam
fixedValue uniform 0;
```

正弦波案例仍使用周期 patch，不受影响。

### 4. 验证

检查：

```bash
grep -A5 -E "xMin|xMax|yMin|yMax" \
    cases/04_solid_rotation_tri_upwind/N100/system/createPatchDict
```

并查看 `log.createPatch`，应出现：

```text
type patch
```

而不是 `type cyclic`。

## Bug 5：三角形旋转案例的后处理不能复用四边形绘图

状态：已修复  
影响范围：三角形固体旋转后处理

### 1. 现象

四边形旋转后处理可以把场值转换成二维数组：

```python
initial.reshape(ny, nx)
```

三角形网格没有规则的 `ny x nx` 排列，直接复用会导致数据布局错误，
或者无法生成图像。

### 2. 根因

三角形案例需要同时使用：

```text
三角形节点坐标
三角形连接关系
OpenFOAM cell centre
OpenFOAM cell field value
```

### 3. 修复

文件：

```text
scripts/common/postprocess_case.py
```

新增三角形旋转后处理分支：

```python
_postprocess_tri_solid_rotation_case(...)
```

新增绘图函数：

```python
plot_tri_rotation_field_comparison(...)
plot_tri_rotation_final_contour(...)
```

cell 值先通过三角形中心匹配到 Gmsh 三角形，再使用：

```python
axis.tripcolor(...)
axis.tricontour(...)
axis.triplot(...)
```

生成真实三角形网格上的图。

### 4. 输出

```text
field_comparison.png
contour_final.png
cfl_history.png
```

## Bug 6：旋转 `N=100` 运行到 `t≈4` 后没有图片

状态：已修复并验证  
日期：2026-08-25  
影响范围：长时间运行的固体旋转案例，四边形和三角形

### 1. 现象

运行：

```bash
python3 scripts/run_study.py \
    --config scripts/configs/04_solid_rotation_tri_upwind.json \
    --resolutions 100 \
    --overwrite
```

案例目录中存在大量时间目录和求解器日志，但没有：

```text
data/cases/04_solid_rotation_tri_upwind/N100/summary.json
figures/cases/04_solid_rotation_tri_upwind/N100/contour_final.png
```

### 2. 日志证据

`log.explicitAdvectionFoamStudent` 中出现：

```text
Time = 3.9989795918371969

FOAM Warning:
Increased the timePrecision from 17 to 18

FOAM FATAL ERROR:
Current time name 4.00000000000046185
The maximum time precision has been reached
```

同时日志显示：

```text
rate max    = 392
CFL deltaT  = 0.00102040816327
```

### 3. 根因

原时间循环使用：

```cpp
runTime++;
```

其内部等价于不断执行：

```text
t^(n+1) = t^n + deltaT
```

由于 `deltaT` 是浮点数，经过数千次累加后，时间值出现微小漂移。
OpenFOAM 为避免两个时间目录被格式化成相同名字，会自动提高
`timePrecision`；当达到最大精度仍无法区分时，程序主动终止。

这不是有限体积离散公式错误，也不是三角形网格质量错误。该案例的
`checkMesh` 已显示：

```text
Mesh OK.
```

### 4. 修复

文件：

```text
UDF/solver/explicitAdvectionFoamStudent/explicitAdvectionFoamStudent.C
```

时间循环改为根据整数步号计算目标时间：

```cpp
const scalar targetTime
(
    min
    (
        startTime + (step + 1)*deltaT,
        endTime
    )
);

const scalar stepDeltaT = targetTime - oldTime;

runTime.setDeltaT(stepDeltaT);
++step;
runTime.setTime(targetTime, step);
```

数学含义：

$$
t_{\mathrm{target}}^k
=
\min\left(
t_0+k\Delta t_{\mathrm{CFL}},
t_{\mathrm{end}}
\right).
$$

这样不再把上一步已经发生误差的时间继续累加，并且最后一步仍会准确截断
到 `endTime=2*pi`。

### 5. 验证

重新编译已经成功：

```bash
source /opt/openfoam14/etc/bashrc
sh scripts/build_student_solver.sh
```

随后使用修复后的求解器完成了四边形和三角形的多网格测试：

```bash
python3 scripts/run_study.py \
    --config scripts/configs/04_solid_rotation_quad_upwind.json \
    --resolutions 50,100,200 \
    --overwrite

python3 scripts/run_study.py \
    --config scripts/configs/04_solid_rotation_tri_upwind.json \
    --resolutions 50,100,200 \
    --overwrite
```

六个算例均满足：

```text
finalTime = 6.283185307179586
meshOK = true
solverEnded = true
solverFatal = false
maxCo = 0.2
```

并且每个网格均生成：

```text
summary.json
time_history.csv
field_data.csv
field_comparison.png
contour_final.png
cfl_history.png
```

结果摘要如下：

| 网格 | N | cell 数 | 一圈后 `cycleL1AgainstInitial` | 归一化质量误差 |
|---|---:|---:|---:|---:|
| quad | 50 | 2,500 | 1.3055035920 | 7.40e-14 |
| quad | 100 | 10,000 | 1.1065381807 | 1.21e-14 |
| quad | 200 | 40,000 | 0.8851521923 | 1.49e-14 |
| tri | 50 | 5,000 | 1.0755771057 | 2.20e-13 |
| tri | 100 | 20,000 | 0.8869767627 | 3.87e-14 |
| tri | 200 | 80,000 | 0.6856004461 | 1.27e-13 |

这里的 `cycleL1AgainstInitial` 是固体旋转一周后与初始场的体积加权
差异，用来观察数值耗散和形状变形；它不是解析解收敛阶。随着网格加密，
该差异下降，说明数值解的旋转轮廓得到改善。质量误差维持在约
`10^-13` 到 `10^-14`，说明周期推进中的守恒检查正常。

### 6. 经验

“程序已经运行了几千步”不等于“案例成功”。长时间案例必须检查：

```text
求解器是否正常结束
最终时间是否正确
最终场是否写出
后处理是否执行
```

只有求解器成功返回，`run_study.py` 才会进入后处理阶段。

## Bug 7：`opal_ifinit` 和 ParaView 配置提示

状态：已知环境警告，不影响当前编译  
影响范围：本机 OpenFOAM 环境初始化

### 1. 现象

编译时曾出现：

```text
opal_ifinit: socket() failed with errno=1
/opt/openfoam14/etc/config.sh/paraview:
[: !=: 需要一元运算符
```

### 2. 判断

随后 `wmake` 仍然完成了目标文件编译和链接，生成：

```text
build/bin/explicitAdvectionFoamStudent
```

因此这两个提示没有导致当前求解器构建失败。

### 3. 处理原则

先看最终退出状态和实际产物：

```bash
test -x build/bin/explicitAdvectionFoamStudent && echo solver-present
```

如果只是环境脚本提示，但编译和运行成功，不要把它误判成求解器算法
错误。若以后它导致 OpenFOAM 工具无法运行，再单独记录为工具链问题。

## 8. 当前排查顺序

遇到“没有图”时，按以下顺序检查：

### 第一步：看案例是否有时间目录

```bash
find cases/<caseName>/N<N> -maxdepth 1 -type d | sort -V | tail
```

没有时间目录，说明求解器可能没有启动。

### 第二步：看求解器日志结尾

```bash
tail -80 cases/<caseName>/N<N>/log.explicitAdvectionFoamStudent
```

重点找：

```text
FOAM FATAL ERROR
Stage 5 time loop completed.
End
```

### 第三步：看后处理目录

```bash
find data/cases/<caseName>/N<N> -maxdepth 1 -type f
find figures/cases/<caseName>/N<N> -maxdepth 1 -type f
```

### 第四步：判断外层 Python 是否提前退出

```text
subprocess.CalledProcessError
    -> Python 只是报告 Allrun 失败
    -> 必须继续向上查看 OpenFOAM 日志
```

### 第五步：按层次分类

```text
Python 报错
    -> 检查脚本参数、路径和 subprocess 调用

FOAM FATAL IO ERROR
    -> 检查 system/ 字典

checkMesh 失败
    -> 检查网格和 patch 拓扑

solver 运行中断
    -> 检查 C++ 时间循环、字段和数值稳定性

没有图但 solver 成功
    -> 检查 postprocess_case.py 和输出路径
```

## 9. 每次修复后的最低验证清单

```bash
# 1. Python 语法
python3 - <<'PY'
import ast
from pathlib import Path

for path in Path("scripts").rglob("*.py"):
    ast.parse(path.read_text(encoding="utf-8"), filename=str(path))
print("Python syntax OK")
PY

# 2. OpenFOAM 求解器编译
source /opt/openfoam14/etc/bashrc
sh scripts/build_student_solver.sh

# 3. 案例网格检查
grep -n "Mesh OK" cases/<caseName>/N<N>/log.checkMesh

# 4. 求解器结束检查
grep -E "Stage 5 time loop completed|final time|End|FOAM FATAL ERROR" \
    cases/<caseName>/N<N>/log.explicitAdvectionFoamStudent

# 5. 后处理产物检查
find data/cases/<caseName>/N<N> figures/cases/<caseName>/N<N> \
    -maxdepth 1 -type f
```

## 10. 目前最重要的经验总结

1. 外层 Python 异常通常不是根因，必须查看 `log.*`。
2. OpenFOAM 字典错误优先检查括号层级和关键字名称。
3. 三角形案例必须以真实 cell centre 和体积为准。
4. 物理边界条件不能只由网格类型决定。
5. 长时间积分要避免反复累加浮点时间。
6. 求解器成功结束是后处理生成图片的前提。
7. 编译成功、网格正确、数值合理、后处理成功是四个不同验收层次。
