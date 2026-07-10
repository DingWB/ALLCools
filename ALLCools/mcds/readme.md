# Beta-Binomial 矩估计（估计甲基化 Beta 先验的 α, β）

## 1. 模型与参数化

每个位点 $i$ 有覆盖 $n_i$、甲基化计数 $k_i$。假设位点真实甲基化率 $p_i \sim \text{Beta}(\alpha, \beta)$，观测 $k_i \mid p_i \sim \text{Binomial}(n_i, p_i)$。

用两个更好解释的参数：

$$
\mu = \frac{\alpha}{\alpha + \beta}\ (\text{均值}), \qquad
\kappa = \alpha + \beta\ (\text{浓度}), \qquad
\rho = \frac{1}{\kappa + 1}\ (\text{位点内相关 / 过离散})
$$

反过来：

$$
\alpha = \mu\kappa, \qquad \beta = (1-\mu)\kappa, \qquad \kappa = \frac{1-\rho}{\rho}
$$

## 2. Beta-Binomial 的边际矩（关键）

$$
E[k_i] = n_i \mu
$$

$$
\operatorname{Var}(k_i) = n_i\, \mu(1-\mu)\, \bigl[\, 1 + (n_i - 1)\rho \,\bigr]
$$

比纯 Binomial 多了过离散因子 $[1 + (n_i - 1)\rho]$——这正是"用 cov 校正"的来源。$\rho = 0$ 退化为 Binomial。

## 3. 矩估计（可处理不等覆盖，推荐）

**第一矩 → μ**（池化，即"全局 mc/cov"）：

$$
\hat\mu = \frac{\sum_i k_i}{\sum_i n_i}
$$

**第二矩 → ρ**：构造一个 Pearson 型过离散统计量

$$
X = \sum_i \frac{(k_i - n_i \hat\mu)^2}{n_i\, \hat\mu(1-\hat\mu)}
$$

因为

$$
E\!\left[\frac{(k_i - n_i \mu)^2}{n_i \mu(1-\mu)}\right] = 1 + (n_i - 1)\rho
$$

所以 $E[X] = N + \rho \sum_i (n_i - 1)$（$N$ = 位点数）。令 $X = E[X]$ 解出

$$
\boxed{\ \hat\rho = \frac{X - N}{\sum_i (n_i - 1)}\ }
$$

（更严的自由度修正把分子 $N$ 换成 $N-1$，即 Tarone 1979 的形式；数据量大时差别可忽略。）

**转成 α, β**：

$$
\hat\kappa = \frac{1 - \hat\rho}{\hat\rho}, \qquad
\hat\alpha = \hat\mu\, \hat\kappa, \qquad
\hat\beta = (1 - \hat\mu)\, \hat\kappa
$$

## 4. 等价的"频率方差"直觉版（覆盖近似相等时）

设 $f_i = k_i / n_i$，则

$$
\operatorname{Var}(f_i) = \mu(1-\mu)\left[\frac{1}{n_i} + \left(1 - \frac{1}{n_i}\right)\rho\right]
$$

若各位点覆盖 $\approx \bar n$，用观测方差 $s^2 = \operatorname{Var}(f_i)$ 解：

$$
\hat\rho = \frac{\dfrac{s^2}{\mu(1-\mu)} - \dfrac{1}{\bar n}}{1 - \dfrac{1}{\bar n}}
$$

这就明确显示了"从率方差里**扣掉 binomial 抽样项 $1/\bar n$**"再得到真实的 Beta 离散——这是它比 ALLCools 纯 Beta-MoM（不扣抽样噪声）更准的地方。覆盖差异大时应该用第 3 节的 $X$ 版本（对每个 $n_i$ 加权），不要用这个近似。

## 5. 边界与实现要点

- **裁剪** $\hat\rho \in (\epsilon,\ 1-\epsilon)$：
  - $\hat\rho \le 0$（欠离散 / 被抽样噪声压过）→ 令 $\rho \to \epsilon$，即 $\kappa$ 很大（强 prior）；
  - $\hat\rho \to 1$ → $\kappa \to 0$（几乎无收缩）。
- $\hat\mu$ 也裁剪到 $(\epsilon, 1-\epsilon)$，避免 $\log$ 发散。
- **过滤 $n_i = 0$**，最好 $n_i \ge 2$（$n_i = 1$ 对 $\sum(n_i - 1)$ 无贡献且噪声大）。
- **CpG / CpH 分开各估一套**（背景完全不同）。
- 在**深测序 pseudobulk 参考**上估，而非浅测序 query cell。

## 附：与 ALLCools 的区别

ALLCools（`calculate_posterior_mc_frac`）用的是**纯 Beta 矩估计**：直接对观测率 $f_i = mc/cov$ 求均值 $\mu$ 和方差 $\sigma^2$，套用

$$
\alpha = \mu\left[\frac{\mu(1-\mu)}{\sigma^2} - 1\right], \qquad
\beta = \alpha\left(\frac{1}{\mu} - 1\right)
$$

它**没有扣除 binomial 抽样噪声**（$\sigma^2$ 里混了抽样方差），所以在低覆盖数据上会高估方差、低估浓度 $\kappa$。上面的 Beta-Binomial 版本通过第 2 节的过离散因子对覆盖做了校正。

> 说明：当前 `ALLCools/mcds/utilities.py` 里的 `calculate_posterior_mc_frac` 已改为本文第 3 节的 Beta-Binomial 矩估计实现；下方 `calculate_posterior_mc_frac_deprecated` 保留的是旧的纯 Beta-MoM。

## 6. 新方法 vs 旧方法：好在哪里

**结论：新方法（Beta-Binomial 矩估计）在单细胞甲基化这种低/不等覆盖的场景下更准，且不会带来额外的显著计算开销。** 具体优势如下。

### 6.1 正确扣除了 binomial 抽样噪声

观测率 $f_i = k_i/n_i$ 的方差包含两部分：

$$
\operatorname{Var}(f_i) = \underbrace{\mu(1-\mu)\rho}_{\text{生物学真实离散}} + \underbrace{\frac{\mu(1-\mu)(1-\rho)}{n_i}}_{\text{binomial 抽样噪声}}
$$

- **旧方法**把整个 $\operatorname{Var}(f_i)$ 都当作 Beta 离散，于是**高估方差、低估浓度 $\kappa$**——先验被当成"比真实更宽"，收缩太弱。
- **新方法**用第 2 节的过离散因子 $[1+(n_i-1)\rho]$ 显式扣掉了 $1/n_i$ 那一项，得到的是**真实的位点间离散 $\rho$**。

覆盖越低（$n_i$ 越小），抽样噪声项越大，旧方法的偏差越严重——而单细胞甲基化恰恰是低覆盖场景。

### 6.2 正确处理不等覆盖（按 $n_i$ 加权）

- 旧方法对每个位点的 $f_i$ **等权**求均值/方差：一个 $cov=1$ 的位点（$f_i$ 只能取 0 或 1，噪声极大）和一个 $cov=100$ 的位点被同等对待。
- 新方法的 $\hat\mu=\sum k_i/\sum n_i$ 是**按覆盖加权**的池化估计，$X$ 统计量对每个 $n_i$ 单独加权，天然地让高覆盖位点贡献更多、低覆盖位点贡献更少。

### 6.3 均值估计更稳、无 $0/0$ 病态

- 旧方法先对每个位点算 $f_i=mc/cov$ 再平均，$cov=0$ 会产生 `NaN`、$cov$ 很小时 $f_i$ 抖动剧烈。
- 新方法先求和再相除（$\sum k_i / \sum n_i$），对 $cov=0$ 的位点自然免疫，估计更平滑。

### 6.4 参数有明确、可解释的物理含义

新方法把先验拆成 $\mu$（背景甲基化水平）和 $\rho$（位点内相关 / 过离散），二者都可解释、可分别裁剪；CpG / CpH 背景差异也能自然体现在各自的 $\rho$ 上。

### 6.5 下游后验一致

先验 $\alpha,\beta$ 更准 $\Rightarrow$ 共轭后验 $\text{Beta}(\alpha+k_i,\ \beta+n_i-k_i)$ 的**收缩强度更合理**：低覆盖位点向背景 $\mu$ 收缩得当，高覆盖位点保留自身信息。后验均值 `post_frac` 和后验标准差 `post_sigma` 的公式不变，但因为 $\alpha,\beta$ 更准，数值也更可靠。

### 6.6 什么时候两者差别不大

当所有位点覆盖都很高且大致相等（$n_i \gg 1$、$\bar n$ 接近）时，抽样噪声项 $1/n_i \to 0$，两种方法趋于一致。**差别主要出现在低覆盖、覆盖高度不均的数据上**——也就是单细胞甲基化的典型情形，因此推荐新方法。

### 局限与注意

- 新方法只估**一套全局 $(\alpha,\beta)$/每个 mc_type**，属于矩估计，不是完整贝叶斯/最大似然；若需要更精细可上 Beta-Binomial MLE，但计算更贵。
- 仍需注意在**深测序 pseudobulk 参考**上估计，而非浅测序 query cell（见第 5 节）。
