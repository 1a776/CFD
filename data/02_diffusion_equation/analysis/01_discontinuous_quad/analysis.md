# 01_discontinuous_quad 网格收敛性分析

本研究使用第二题第一个具有间断初值的二维扩散算例，比较不同 均匀结构化四边形 网格分辨率。
所有网格使用相同的区域 [-5,5]^2、扩散系数 μ=1、齐次 Neumann 边界、显式扩散稳定时间步和终止时间 t=0.2。

$$p = \frac{\log(E_N/E_{2N})}{\log(2)}$$

其中 $E_N$ 是分辨率为 $N$ 时的误差，$p$ 是观察收敛阶。

## 汇总表

| N | cells | L1 | L2 | Linf | L1 order | diffusionCo | final range |
|---:|---:|---:|---:|---:|---:|---:|---:|
| 10 | 100 | 2.10899769e-01 | 1.72401247e-01 | 1.53398257e-01 | - | 0.450000 | 6.95532413e-01 |

## 结论

L1 误差从 `2.10899769e-01` 变为 `2.10899769e-01`。
当前只有一个有效网格，暂时无法计算观察收敛阶。
每个 N 的详细场数据、时间历史、日志和单案例图保存在 `cases/02_diffusion_equation/<case>/Nxx`、`data/02_diffusion_equation/cases/<case>/Nxx` 和 `figures/02_diffusion_equation/cases/<case>/Nxx` 目录。
