# 调试指南

## 1. 编译错误先看哪一层

### 找不到头文件

检查：

- 是否执行了 `source /opt/openfoam14/etc/bashrc`；
- `Make/options` 是否包含 `finiteVolume/lnInclude`；
- 新增 API 后是否加入了对应头文件。

### 找不到符号或链接失败

检查 `Make/options` 中是否链接：

```make
-lfiniteVolume
-lmeshTools
```

### 移动或重命名源文件后出现奇怪依赖

在求解器目录执行：

```bash
wclean
```

然后重新执行构建脚本。

## 2. 案例错误先区分类型

### `cannot find file .../0/U`

说明 `0.orig/U` 没有复制到 `0/U`，或者 `Allclean` 后没有重新执行：

```bash
cd /home/a776/workdocuments/上交船舶/slover/student_project
python3 scripts/01_sine_wave_quad/create_initial_fields.py \
    --case-dir cases/01_sine_wave_quad/N20
cp -R cases/01_sine_wave_quad/N20/0.orig cases/01_sine_wave_quad/N20/0
```

### `cannot find T`

说明初始场脚本没有运行，或者 `0.orig/T` 没有生成。

### 周期边界错误

检查 `blockMeshDict` 的：

- `neighbourPatch`；
- 两个周期 patch 的面数量；
- `U` 和 `T` 的 patch 类型；
- `zMin/zMax` 是否为 `empty`。

### `div(phi,T) not found`

检查 C++ 中散度名称和 `fvSchemes` 是否完全一致：

```text
div(phi,T)
```

大小写、括号和逗号都不能随意改变。

## 3. 数值结果异常

### `Co` 超过目标值

检查：

- 通量率是否除以 `mesh.V()`；
- CFL 公式中的 `1/2` 是否处理一致；
- 是否在更新之前设置了正确的 `deltaT`；
- 最后一步是否截断到剩余时间。

### 场值迅速爆炸

检查：

- 更新符号是否是减号；
- 是否误用了 `fvm::div`；
- 是否把旧时间层提前覆盖；
- `phi` 是否是 `fvc::flux(U)`；
- `maxCo` 是否为正数。

### 场值变得很平但没有爆炸

这通常是迎风格式的数值耗散，不一定是程序错误。
要用网格收敛和误差脚本判断，而不是只看一组网格。

### 质量不守恒

检查：

- 对流项是否使用守恒形式；
- 周期 patch 是否正确配对；
- 是否只更新了内部场但没有正确处理边界；
- 积分时是否乘了 `mesh.V()`。

## 4. 调试原则

每次只改变一个因素，并保留日志：

```bash
sh scripts/build_student_solver.sh
sh cases/01_sine_wave_quad/N20/Allrun
```

把“代码改动、编译结果、运行结果、你的解释”一起记录。
这样才能判断问题来自 C++、字典、网格，还是数学公式。
