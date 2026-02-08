# VCmax-Nitrogen calculations
# 
# VCmax25 = 62.5          # maximum capacity of Rubisco at 25C-VCmax	
# N_leaf = 3.10 + 1.35    # leaf Nitrogen content	mean value + 1 SD g/m2 
# slope = 20.72 / 62.5       # slope of VCmax-N curve
function VCmax(lai::FT, Ω::FT, CosZs::FT, VCmax25::FT, N_leaf::FT, slope::FT) where {FT<:Real}
  CosZs <= 0 && return 0.0, 0.0 # 光合仅发生在白天
  K = 0.5 * Ω / CosZs # assuming a spherical leaf angle distribution
  Kn = 0.3            # 0.713/2.4

  expr1 = 1 - exp(-K * lai)
  expr2 = 1 - exp(-lai * (Kn + K))
  expr3 = 1 - exp(-Kn * lai)

  # Formulas based on Chen et al., 2012, GBC
  if (expr1 > 0)
    VCmax_sunlit = VCmax25 * slope * N_leaf * K * expr2 / (Kn + K) / expr1
  else
    VCmax_sunlit = VCmax25
  end

  if (K > 0 && lai > expr1 / K)
    VCmax_shaded = VCmax25 * slope * N_leaf *
                   (expr3 / Kn - expr2 / (Kn + K)) / (lai - expr1 / K)
  else
    VCmax_shaded = VCmax25
  end
  VCmax_sunlit, VCmax_shaded
end


"""
    VCmax_profile(nlyr::Int, lai::FT, Ω::FT, CosZs::FT, VCmax25::FT, N_leaf::FT, slope::FT) where {FT<:Real}

计算冠层Vcmax的垂向分布（多层模型）

# 参数
- `nlyr::Int`: 冠层垂向层数
- `lai::FT`: 总叶面积指数 [m² m⁻²]
- `Ω::FT`: 遮挡指数（clumping index）[-]
- `CosZs::FT`: 太阳天顶角余弦 [-]
- `VCmax25::FT`: 25°C下的最大羧化速率 [μmol m⁻² s⁻¹]
- `N_leaf::FT`: 冠层平均叶片氮含量 [g m⁻²]
- `slope::FT`: Vcmax-N曲线的斜率 [μmol g⁻¹ s⁻¹]

# 返回值
- `VCmax_sunlit::Vector{FT}`: 每层向阳叶的Vcmax [μmol m⁻² s⁻¹]
- `VCmax_shaded::Vector{FT}`: 每层背阴叶的Vcmax [μmol m⁻² s⁻¹]
- `LAI_sunlit::Vector{FT}`: 每层向阳叶的LAI [m² m⁻²]
- `LAI_shaded::Vector{FT}`: 每层背阴叶的LAI [m² m⁻²]

# 原理
基于Chen et al., 2012, GBC的氮素垂向分布模型：
- 叶片氮含量随累积LAI呈指数衰减：N(L) = N₀ * exp(-Kn * L)
- Vcmax与氮含量成正比：Vcmax(L) = Vcmax25 * slope * N(L)
- 考虑光的消光效应，区分向阳和背阴叶

# 示例
```julia
nlyr = 10
lai = 5.0
Ω = 1.0
CosZs = 0.866  # 30° 太阳天顶角
VCmax25 = 62.5
N_leaf = 3.10
slope = 20.72 / 62.5

VCmax_sunlit, VCmax_shaded, LAI_sunlit, LAI_shaded =
    VCmax_profile(nlyr, lai, Ω, CosZs, VCmax25, N_leaf, slope)
```
"""
function VCmax_profile(
  nlyr::Int,
  lai::FT,
  Ω::FT,
  CosZs::FT,
  VCmax25::FT,
  N_leaf::FT,
  slope::FT
) where {FT<:Real}

  # 初始化输出数组
  VCmax_sunlit = zeros(FT, nlyr)
  VCmax_shaded = zeros(FT, nlyr)
  LAI_sunlit = zeros(FT, nlyr)
  LAI_shaded = zeros(FT, nlyr)

  # 每层的LAI
  dLAI = lai / nlyr

  # 夜间处理：没有光合作用，所有LAI都是背阴叶
  if CosZs <= 0
    fill!(LAI_shaded, dLAI)
    return VCmax_sunlit, VCmax_shaded, LAI_sunlit, LAI_shaded
  end

  # 消光系数
  K = 0.5 * Ω / CosZs  # 光的消光系数（假设球形叶倾角分布）
  Kn = 0.3              # 氮素的消光系数 (0.713/2.4)

  # 对每一层进行计算（从冠层顶部到底部）
  for i in 1:nlyr
    # 当前层顶部的累积LAI（从冠层顶部开始）
    L_top = (i - 1) * dLAI
    # 当前层底部的累积LAI
    L_bottom = i * dLAI
    # 当前层中心的累积LAI
    L_mid = (L_top + L_bottom) / 2

    # 计算当前层的向阳和背阴叶面积
    # 基于 Beer 定律微分：每层向阳叶面积 = dLAI * K * exp(-K * L_mid)
    # 这表示在深度 L_mid 处，叶片被光照射的比例

    if K > 1e-6
      # 微分形式：每层中点处的向阳叶比例
      LAI_sunlit[i] = dLAI * K * exp(-K * L_mid)
    else
      # 当 K 很小时，所有叶片都是向阳的
      LAI_sunlit[i] = dLAI
    end

    # 背阴叶面积 = 总叶面积 - 向阳叶面积
    LAI_shaded[i] = dLAI - LAI_sunlit[i]

    # 确保非负
    LAI_sunlit[i] = max(0.0, LAI_sunlit[i])
    LAI_shaded[i] = max(0.0, LAI_shaded[i])

    # 计算当前层中心的氮含量（基于氮素垂向分布）
    # N(L) = N₀ * exp(-Kn * L)
    # 这里N_leaf是冠层平均氮含量，需要转换为冠层顶部氮含量N₀
    # ∫N₀*exp(-Kn*L)dL from 0 to LAI = N₀/Kn * (1 - exp(-Kn*LAI))
    # N_leaf * LAI = N₀/Kn * (1 - exp(-Kn*LAI))
    # N₀ = N_leaf * LAI * Kn / (1 - exp(-Kn*LAI))

    if Kn > 1e-6 && lai > 1e-6
      N0 = N_leaf * lai * Kn / (1 - exp(-Kn * lai))
    else
      N0 = N_leaf
    end

    # 当前层的氮含量
    N_layer = N0 * exp(-Kn * L_mid)

    # 基于氮含量计算Vcmax
    # Vcmax = Vcmax25 * slope * N
    # 这里的slope是归一化的斜率（相对于参考氮含量）
    VCmax_layer = VCmax25 * slope * N_layer

    # 向阳和背阴叶的Vcmax可以认为是相同的（同一层的叶片氮含量相同）
    # 但也可以考虑光适应效应，这里先用相同的值
    VCmax_sunlit[i] = VCmax_layer
    VCmax_shaded[i] = VCmax_layer

    # 另一种方法：考虑光适应，向阳叶的Vcmax略高
    # 这里使用简化假设：同层叶片Vcmax相同
  end

  return VCmax_sunlit, VCmax_shaded, LAI_sunlit, LAI_shaded
end


"""
    VCmax_profile_mean(nlyr::Int, lai::FT, Ω::FT, CosZs::FT, VCmax25::FT, N_leaf::FT, slope::FT) where {FT<:Real}

计算考虑垂向分布后的冠层平均Vcmax（用于大叶模型）

# 返回值
- `VCmax_mean::FT`: 冠层平均Vcmax，按LAI加权 [μmol m⁻² s⁻¹]
- `VCmax_sunlit_mean::FT`: 向阳叶平均Vcmax [μmol m⁻² s⁻¹]
- `VCmax_shaded_mean::FT`: 背阴叶平均Vcmax [μmol m⁻² s⁻¹]

# 示例
```julia
VCmax_mean, VCmax_sunlit_mean, VCmax_shaded_mean =
    VCmax_profile_mean(10, 5.0, 1.0, 0.866, 62.5, 3.10, 20.72/62.5)
```
"""
function VCmax_profile_mean(
  nlyr::Int,
  lai::FT,
  Ω::FT,
  CosZs::FT,
  VCmax25::FT,
  N_leaf::FT,
  slope::FT
) where {FT<:Real}

  # 获取垂向分布
  VCmax_sunlit, VCmax_shaded, LAI_sunlit, LAI_shaded =
    VCmax_profile(nlyr, lai, Ω, CosZs, VCmax25, N_leaf, slope)

  # 计算加权平均
  total_sunlit = sum(LAI_sunlit)
  total_shaded = sum(LAI_shaded)

  if total_sunlit > 1e-6
    VCmax_sunlit_mean = sum(VCmax_sunlit .* LAI_sunlit) / total_sunlit
  else
    VCmax_sunlit_mean = 0.0
  end

  if total_shaded > 1e-6
    VCmax_shaded_mean = sum(VCmax_shaded .* LAI_shaded) / total_shaded
  else
    VCmax_shaded_mean = 0.0
  end

  # 总体平均
  if lai > 1e-6
    VCmax_mean = (sum(VCmax_sunlit .* LAI_sunlit) + sum(VCmax_shaded .* LAI_shaded)) / lai
  else
    VCmax_mean = VCmax25
  end

  return VCmax_mean, VCmax_sunlit_mean, VCmax_shaded_mean
end
