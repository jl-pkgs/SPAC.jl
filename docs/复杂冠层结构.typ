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

- 需要知道地表温度。

- 如何构建一个地表温度，替代模型？
  
  基于5公分土壤温度数据？
