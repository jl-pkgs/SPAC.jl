#import "@local/modern-cug-report:0.1.3": *
#show: doc => template(doc, size: 12pt, footer: "CUG水文气象学2025", header: "")

#let int = $integral$
#let VCmax = $"V"_{c max}$
#let VCmax0 = $"V"_{c max, 0}$


= 1 *复杂冠层结构*

== 1.1 羧化能力$bold(VCmax)$

N含量（$g thin m^{-2}$）指数下降，顶部N含量最高。$chi_n$定义了$VCmax$对N含量的敏感性，可认为是回归系数。


$ N(L) = N_0 e^{-k_n L} $

$ VCmax = VCmax0 chi_n N(L) $

$ f_"sun" (L)= e^{-k L} $

$ f_"sha" (L) = 1 - e^{-k L} $

辐射衰减系数，Assuming a spherical leaf angle distribution：

$ k = 0.5 Ω / "CosZs" $

// 求每一层的$"V"_{c max, "sun"}$和$"V"_{c max, "sha"}$

*阳叶*

$
  "V"_{c max, "sun"} (L)
  = & (int_0^L VCmax0 chi_n N(L) f_"sun" (L) d L ) / (int_0^L f_"sun" (L) d L ) \
  = & (int_0^L VCmax0 chi_n N_0 e^{-k_n L} e^{-k L} d L ) / (int_0^L e^{-k L} d L ) \
  = & (VCmax0 chi_n N_0 (e^{-(k_n + k) L} - 1) / (-(k_n + k))) / ((e^{-k L} - 1) / (-k)) \
  = & #h(0.4em) boxed(VCmax0 chi_n N_0 k (1 - e^{-(k_n + k) L}) / ( (k_n + k) (1 - e^{-k L}) ))
$

*阴叶*

$
  "V"_{c max, "sha"} (L)
  &= (int_0^L VCmax0 chi_n N(L) f_"sha" (L) d L ) / (int_0^L f_"sha" (L) d L ) \
  &= (int_0^L VCmax0 chi_n N_0 e^{-k_n L} (1 - e^{-k L}) d L ) / (int_0^L (1 - e^{-k L}) d L ) \
  &= #h(0.4em) boxed(VCmax0 chi_n N_0 k ( ( (1 - e^{-k_n L}) / k_n - (1 - e^{-(k_n + k) L}) / (k_n + k) )) / (k L - (1 - e^{-k L}) ))
$

每一层的$"V"_{c max, "sha"} (L)$与$"V"_{c max, "sun"} (L)$为：

$ "V"_{c max, "sha", i} = "V"_{c max, "sha"} (L_i^{+}) - "V"_{c max, "sha"} (L_i^{-}) $

$ "V"_{c max, "sun", i} = "V"_{c max, "sun"} (L_i^{+}) - "V"_{c max, "sun"} (L_i^{-}) $

其中$L_i^{+}$与$L_i^{-}$分别为第i层的上界与下界。


=== 1.1.1 状态变量state

- 胞内CO2浓度：$c_{i, "sha"}, c_{i, "sun"}$


// == 参数params

== 1.3 辐射传输

本阶段先不用 Norman 全参数方案，先采用简化 Beer 方法，重点解决阳叶/阴叶辐射比例。

给定总短波辐射 $R_s$ 与总叶面积指数 $L = L_"sun" + L_"sha"$。

设散射比例为 $f_"dif"$（固定参数，默认 0.2）：

$
  R_"dir" = (1 - f_"dif") R_s,
  quad R_"dif" = f_"dif" R_s
$

简化分配规则：

- 直射部分优先进入阳叶，但允许一小部分分给阴叶（表示多次散射）；
- 散射吸收按叶面积比例在阳叶与阴叶之间分配。

写成单位叶面积辐射：

$
  R_"sun" = (R_"dir,sun" + R_"dif" L_"sun" / L) / L_"sun"
$

$
  R_"sha" = (R_"dir,sha" + R_"dif" L_"sha" / L) / L_"sha"
$

该方案满足守恒：

$
  L_"sun" R_"sun" + L_"sha" R_"sha" = R_s
$

优势是参数少、可辨识性高，可直接用于 TwoLeaf / TwoBigLeaf 的光合与导度计算。

实现上将 `enable_partition` 作为开关：默认 `false` 保持历史兼容（阳叶/阴叶接收相同辐射）；
当 `enable_partition = true` 时启用阳叶/阴叶辐射差异。
