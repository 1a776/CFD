# Stage 6：结果可视化、误差和网格收敛

本文说明当前教学工程如何把 OpenFOAM 求解结果转换成：

- 直观的二维场图；
- 数值解和精确解的剖面对比；
- 波幅随时间变化图；
- Courant 数随时间变化图；
- 多个网格分辨率下的误差收敛图；
- 可以写进研究计划或汇报材料的 `summary.json` 和 CSV 数据。

本阶段的核心问题是：

> 求解器虽然能运行，但它算出的结果是否合理？误差是否随着网格加密而下降？图像能不能解释数值方法产生了什么误差？

## 1. 本次修改总览

本次修改涉及四类文件：

| 文件 | 修改内容 | 针对的问题 |
|---|---|---|
| `UDF/solver/explicitAdvectionFoamStudent/explicitAdvectionFoamStudent.C` | 改进时间循环结束判据 | 避免浮点误差导致额外的伪时间步 |
| `scripts/01_sine_wave_quad/plot_results.py` | 新增单算例后处理 | 只有日志数值，不容易直观看场和耗散 |
| `scripts/01_sine_wave_quad/plot_convergence.py` | 新增多网格收敛分析 | 单个网格不能证明数值方法收敛 |
| 本文档及索引文档 | 解释公式、代码、图像和验收标准 | 让结果可以被理解、复核和汇报 |

本次没有修改：

- `/opt/openfoam14/src`；
- OpenFOAM 的 `fvc::div`、`fvc::flux` 或迎风格式源码；
- 第 1 题的数学模型；
- 求解器的有限体积离散公式；
- 已有的网格和边界条件定义。

## 2. Stage 5 时间循环的修正

### 2.1 原来的问题

原来的循环条件是：

```cpp
while (runTime.value() < endTime - SMALL)
```

在当前算例中，理论上每一步为：

$$\Delta t=0.005.$$

因此从 $0$ 到 $1$ 应该正好有 $200$ 步。但是浮点数不能精确表示所有十进制小数。经过 $200$ 次加法后，内存中的时间可能不是数学意义上精确的 $1$，而是略小于或略大于 $1$。

如果内存中的时间仍被判断为“小于终止时间”，求解器就可能再进入一次循环。此时剩余时间类似：

```text
5.55111512313e-16
```

于是日志会出现一个没有实际物理意义的尾步：

```text
Time = 1  step = 201  deltaT = 5.55111512313e-16
```

这不是 PDE 离散公式错误，而是时间比较中的浮点舍入问题。

### 2.2 现在的修改

现在在时间循环前定义相对终止容差：

```cpp
const scalar timeTolerance
(
    1.0e-12*max(1.0, mag(endTime))
);
```

循环条件改为：

```cpp
while (endTime - runTime.value() > timeTolerance)
```

它对应的数学判断是：

$$t_{\mathrm{end}}-t^n>\varepsilon_t.$$

其中：

$$\varepsilon_t=10^{-12}\max(1,|t_{\mathrm{end}}|).$$

这意味着只有当“真正剩余的时间”大于一个与终止时间尺度相关的容差时，才继续推进。

### 2.3 为什么采用相对容差

直接使用固定的 `SMALL` 不够清晰，因为它是 OpenFOAM 内部的通用小量，不一定和当前时间尺度匹配。

现在的写法有三个优点：

1. 对 `endTime = 1`，容差约为 $10^{-12}$；
2. 对更大的终止时间，容差会随时间尺度适当放大；
3. 对接近终点的浮点舍入误差不会误开一个新时间步。

这次修改不改变正常时间步，只处理循环终点附近的数值判断。

### 2.4 你重新运行后应该看到什么

源码已经修改，但当前目录里已有的旧日志仍然记录的是修改前的运行结果。因此你需要重新编译并运行后再判断。

重新运行后，当前案例预期为：

```text
time steps = 200
final time = 1
```

而不是旧日志中的：

```text
time steps = 201
```

## 3. 单算例后处理脚本

文件：

```text
scripts/01_sine_wave_quad/plot_results.py
```

这个脚本只读取结果，不修改求解器，也不重新计算 OpenFOAM 场。

它读取两类输入：

### 3.1 OpenFOAM 时间目录中的 `T`

脚本读取：

```text
0/T
1/T
```

更一般地，它会寻找数值最大的时间目录作为最终时间：

```text
0.05/T
0.10/T
...
1/T
```

`0/T` 提供初始数值场，最终时间目录中的 `T` 提供数值解。

### 3.2 求解器日志

脚本读取：

```text
log.explicitAdvectionFoamStudent
```

并从每一步日志中提取：

```text
Time
step
deltaT
maxCo
residual integral
T min
T max
```

因此图像不是凭空生成的，而是由求解器已经输出的时间推进数据生成的。

## 4. 图像是怎样生成的

完整数据流是：

```text
0/T 和最终时间/T
    -> Python 解析 internalField
    -> 按 20 x 20 重排成二维数组
    -> 计算精确解和误差
    -> matplotlib 绘图
    -> figures/cases/01_sine_wave_quad/N20/*.png
```

同时：

```text
log.explicitAdvectionFoamStudent
    -> 正则表达式提取每一步的日志行
    -> 生成 time_history.csv
    -> 绘制 amplitude_history.png 和 cfl_history.png
```

当前案例是结构化的 $20\times20$ 网格。初值脚本按“$j$ 外层、$i$ 内层”写入单元值：

$$x_i=\frac{i+1/2}{N_x},\qquad y_j=\frac{j+1/2}{N_y}.$$

后处理脚本用同样的顺序把一维 `internalField` 重新整理成二维数组：

```python
values.reshape(NY, NX)
```

这一步很重要。如果数组顺序和网格单元顺序不一致，图像会出现条纹错位，即使场值本身没有错。

## 5. 四联场图

生成文件：

```text
figures/cases/01_sine_wave_quad/N20/field_comparison.png
```

四个子图分别是：

1. 初始场；
2. 最终数值场；
3. 同一时刻的精确解；
4. 数值解减精确解的误差场。

### 5.1 精确解

本算例速度为：

$$\boldsymbol U=(1,1,0).$$

初值为：

$$T(x,y,0)=\sin(2\pi(x+y)).$$

因此精确解为：

$$T(x,y,t)=\sin(2\pi(x+y-2t)).$$

脚本中的 `exact_field(final_time)` 就是在计算这个公式。

### 5.2 颜色范围

数值场、初始场和精确解使用相同的颜色范围。这样颜色深浅可以直接比较，而不会因为每张图自动缩放色标而产生错觉。

误差图使用以零为中心的对称色标：

$$e_c=T_c^{\mathrm{num}}-T_c^{\mathrm{exact}}.$$

红色和蓝色分别表示正误差和负误差，颜色接近白色表示误差接近零。

### 5.3 应该怎么看

重点观察以下现象：

- 数值波形的位置是否和精确解一致；
- 数值波峰和波谷是否变浅；
- 误差是否主要集中在波峰、波谷或相位变化快的区域；
- 是否出现超出初始范围的异常颜色；
- 是否出现棋盘格、孤立尖峰或明显条纹。

一阶迎风格式通常会产生数值耗散。它的典型表现是：

```text
波形位置大体正确，但波峰和波谷逐渐变平。
```

如果出现明显的波形整体错位，则除了耗散，还可能存在相位误差。

如果出现越来越大的孤立峰值、颜色范围不断扩大，才更像是稳定性或符号错误。

## 6. 对角线剖面图

生成文件：

```text
figures/cases/01_sine_wave_quad/N20/diagonal_profile.png
```

脚本取 $x\approx y$ 的单元：

```python
numerical[index, index]
```

然后画出数值值和精确值的曲线。

这张图把二维场压缩成一条曲线，适合回答：

- 波峰位置有没有移动；
- 波幅衰减了多少；
- 数值曲线是否平滑；
- 数值误差是整体偏差还是局部偏差。

二维色彩图适合看空间分布，剖面图适合做定量比较。两者是互补的。

## 7. 振幅历史图

生成文件：

```text
figures/cases/01_sine_wave_quad/N20/amplitude_history.png
```

脚本用日志中的最大值和最小值定义当前波幅：

$$A^n=\frac{\max_c(T_c^n)-\min_c(T_c^n)}{2}.$$

对于理想的周期平移正弦波，连续精确解的波幅始终为 $1$：

$$A_{\mathrm{exact}}(t)=1.$$

因此图中如果数值波幅从 $1$ 逐渐下降，说明存在数值耗散。

当前旧的 `20\times20` 结果中，最终场范围约为：

```text
min(T) = -0.20524287087
max(T) =  0.20524287087
```

所以最终波幅约为 $0.205$。这说明当前粗网格和一阶迎风格式耗散很强，但并不等于求解器发散。

## 8. CFL 历史图

生成文件：

```text
figures/cases/01_sine_wave_quad/N20/cfl_history.png
```

脚本使用日志中的：

```text
maxCo
```

并和目标值 `0.2` 画在一起。

本工程采用的 CFL 定义是：

$$\mathrm{Co}_{\max}=\frac{\Delta t}{2}\max_c\left(\frac{\sum_f|F_{cf}|}{V_c}\right).$$

验收时应该看到：

- 正常时间步的 `maxCo` 接近 `0.2`；
- 最后一步如果为了精确抵达 `endTime` 而缩短，`maxCo` 可以小于 `0.2`；
- 不应长期超过目标值；
- 不应出现 `nan` 或 `inf`。

CFL 图回答的是稳定性和时间步控制问题，不回答空间精度问题。

## 9. 数值误差的定义

在最终时间，把数值解和精确解之差定义为：

$$e_c=T_c^{\mathrm{num}}-T_c^{\mathrm{exact}}.$$

本工程使用体积加权归一化 L1 误差：

$$E_{L1}=\frac{\sum_cV_c|e_c|}{\sum_cV_c|T_c^{\mathrm{exact}}|}.$$

由于当前案例所有单元体积相同，体积因子在分子和分母中会抵消；但代码仍然保留体积加权写法，因为这才是一般有限体积网格上正确的定义。

同时脚本还计算：

$$E_{L2}=\left(\frac{\sum_cV_c e_c^2}{\sum_cV_c(T_c^{\mathrm{exact}})^2}\right)^{1/2}.$$

以及归一化最大误差：

$$E_{L\infty}=\frac{\max_c|e_c|}{\max_c|T_c^{\mathrm{exact}}|}.$$

三个误差的侧重点不同：

| 指标 | 主要反映 |
|---|---|
| L1 | 全域平均误差，最适合做总体比较 |
| L2 | 对较大误差更敏感 |
| Linf | 最坏单元误差 |

## 10. 为什么单个网格不能证明收敛

在 $20\times20$ 网格上得到一个误差，只能说明：

```text
这个网格、这个时间步、这个格式，在这个案例上的表现。
```

它不能证明：

- 网格加密后误差一定下降；
- 方法具有一阶精度；
- 误差主要来自空间离散还是时间离散；
- 结果不是偶然的。

因此需要至少使用：

```text
N = 10, 20, 40, 80
```

并比较：

$$E_N,\qquad E_{2N}.$$

## 11. 观察收敛阶

当网格分辨率从 $N$ 加密到 $2N$ 时，观察收敛阶定义为：

$$p_N=\frac{\log(E_N/E_{2N})}{\log 2}.$$

更一般地，如果分辨率比值为 $r$：

$$p_N=\frac{\log(E_N/E_{rN})}{\log r}.$$

如果误差满足：

$$E_N\approx C h^p,$$

并且 $h\propto1/N$，则：

$$E_N/E_{2N}\approx2^p.$$

所以：

- $p\approx1$ 表示一阶收敛；
- $p\approx2$ 表示二阶收敛；
- $p<0$ 表示加密后误差反而变大，需要检查；
- 前几层网格的阶数偏低并不罕见，应该观察趋势。

本题使用一阶迎风空间格式和前向 Euler 时间推进，理论上不能期待二阶总体精度。

## 12. 多网格收敛脚本

文件：

```text
scripts/01_sine_wave_quad/plot_convergence.py
```

输入文件：

```text
data/analysis/01_sine_wave_quad/raw_results.csv
```

该文件由 `collect_results.py` 根据各 N 目录中的 `summary.json` 生成。当前脚本入口不接收 CSV 路径参数，而是根据案例名称固定读取上述项目路径。

脚本执行以下步骤：

1. 读取 `raw_results.csv`；
2. 按 `resolution` 从小到大排序；
3. 根据相邻 N 的误差计算 `observedOrder`；
4. 写出 `convergence_summary.csv`；
5. 生成对数坐标误差图；
6. 生成观察收敛阶图和所有 N 的对比图。

输出：

```text
data/analysis/01_sine_wave_quad/convergence_summary.csv
figures/analysis/01_sine_wave_quad/convergence_errors.png
figures/analysis/01_sine_wave_quad/convergence_order.png
figures/analysis/01_sine_wave_quad/all_N_comparison.png
```

误差图使用 log-log 坐标。若一阶参考线与数据曲线大致平行，说明误差随网格加密呈现接近一阶的下降趋势。

## 13. 当前已有收敛结果的解释

下面的数值是历史教学示例，不代表当前目录中仍然保留这些多网格结果。当前多网格结果已按用户要求清理；重新运行四个 N 后，应以新生成的 `data/analysis/01_sine_wave_quad/convergence_summary.csv` 为准。

历史示例为：

| $N$ | 归一化 L1 误差 | 观察阶 |
|---:|---:|---:|
| 10 | 0.9661795 | - |
| 20 | 0.7980079 | 0.2759 |
| 40 | 0.5469567 | 0.5450 |
| 80 | 0.3263462 | 0.7450 |

这组结果可以作如下判断：

1. 误差随着网格加密持续下降；
2. 观察阶从约 `0.28` 增长到约 `0.75`；
3. 收敛趋势逐渐接近一阶；
4. 目前还不能把 `0.75` 直接写成严格理论阶，因为网格仍可能不够细，且时间步随网格变化；
5. 结果支持“方法具有趋近一阶的收敛趋势”，但不支持“已经达到理想的一阶渐近区间”这种更强的表述。

这是一种合格的数值分析表达：结论和证据强度相匹配。

## 14. 守恒和精度不是同一件事

周期边界下，内部面通量应成对抵消。因此离散总量：

$$M^n=\sum_cV_cT_c^n.$$

应该近似保持不变：

$$M^{n+1}-M^n\approx0.$$

但守恒不代表误差很小。

当前 $20\times20$ 结果可以同时满足：

- `residual integral` 接近 $0$；
- 总量误差很小；
- 波幅从 $1$ 衰减到约 $0.205$；
- L1 误差仍然较大。

原因是：

```text
守恒约束控制总量的收支；
精度指标衡量局部场和精确解的差别。
```

一个数值方法可以严格守恒，但由于数值耗散或相位误差，局部解仍然不够准确。

## 15. 运行命令

先重新编译并运行当前学生版求解器：

```bash
source /opt/openfoam14/etc/bashrc
cd /home/a776/workdocuments/上交船舶/slover/student_project
sh scripts/build_student_solver.sh
build/bin/explicitAdvectionFoamStudent -help
sh cases/01_sine_wave_quad/N20/Allrun
```

然后生成单算例图片：

```bash
python3 scripts/01_sine_wave_quad/plot_results.py
```

生成多网格收敛图：

```bash
cd /home/a776/workdocuments/上交船舶/slover/student_project
python3 scripts/01_sine_wave_quad/collect_results.py
python3 scripts/01_sine_wave_quad/plot_convergence.py
```

## 16. 单算例验收标准

执行：

```bash
cd /home/a776/workdocuments/上交船舶/slover/student_project
python3 scripts/01_sine_wave_quad/plot_results.py
```

应生成：

```text
data/cases/01_sine_wave_quad/N20/summary.json
data/cases/01_sine_wave_quad/N20/time_history.csv
figures/cases/01_sine_wave_quad/N20/field_comparison.png
figures/cases/01_sine_wave_quad/N20/diagonal_profile.png
figures/cases/01_sine_wave_quad/N20/amplitude_history.png
figures/cases/01_sine_wave_quad/N20/cfl_history.png
```

检查：

- `summary.json` 中 `finalTime` 接近 `1`；
- `normalizedL1` 是有限数；
- `maxCo` 不超过目标值的明显范围；
- `maxAbsResidualIntegral` 很小；
- 四联图能够显示初值、数值解、精确解和误差；
- 振幅图不应出现发散；
- CFL 图不应出现异常尖峰；
- 对角线图中数值曲线和精确曲线能够直接比较。

## 17. 多网格验收标准

执行：

```bash
python3 /home/a776/workdocuments/上交船舶/slover/student_project/scripts/01_sine_wave_quad/collect_results.py
python3 /home/a776/workdocuments/上交船舶/slover/student_project/scripts/01_sine_wave_quad/plot_convergence.py
```

应生成：

```text
data/analysis/01_sine_wave_quad/convergence_summary.csv
figures/analysis/01_sine_wave_quad/convergence_errors.png
figures/analysis/01_sine_wave_quad/convergence_order.png
figures/analysis/01_sine_wave_quad/all_N_comparison.png
```

检查：

- `N=10,20,40,80` 都有数据；
- L1 误差总体随 $N$ 增大而下降；
- 观察阶为有限数；
- 收敛曲线整体接近一阶参考线；
- `meshOK=True`；
- `solverEnded=True`；
- `solverFatal=False`；
- 最大 CFL 数保持在目标值附近。

## 18. 当前阶段的边界

本阶段完成的是：

- 单算例结果解释；
- 误差定义；
- 多网格收敛分析；
- 结果图和数据表的自动生成。

本阶段尚未完成：

- 三角形网格的学生版求解器案例；
- 固体旋转案例；
- 时间步和空间网格分别独立加密的误差分离实验；
- 二阶空间格式对比；
- ParaView 交互式场景。

因此当前最稳妥的汇报表述是：

> 已完成守恒型线性对流方程的一阶迎风显式有限体积求解器，并在周期正弦波算例上验证了时间推进、周期守恒、场值有界性和网格加密下的误差下降趋势。当前粗网格结果存在明显数值耗散，观察收敛阶逐渐接近一阶。
