export Photosynthesis_Rong2018
export photosynthesis


@bounds @units @with_kw mutable struct Photosynthesis_Rong2018{FT} <: AbstractPhotosynthesisModel{FT}
  "initial slope of the light response curve to assimilation rate, (i.e., quantum efficiency)"
  α::FT = 0.06 | (0.01, 0.10) | "μmol CO2 [μmol PAR]⁻¹"

  "initial slope of the CO2 response curve to assimilation rate, (i.e., carboxylation efficiency)"
  η::FT = 0.04 | (0.01, 0.07) | "μmol m⁻² s⁻¹ [μmol m⁻² s⁻¹]⁻¹"

  "carbon saturated rate of photosynthesis at 25 °C"
  VCmax25::FT = 50.00 | (5.00, 120.00) | "μmol m⁻² s⁻¹"

  "photoperiod constraint"
  d_PC::FT = 2.0 | (0.0, 5.0) | "-"

  "extinction coefficients for visible radiation" # 植被光合参数
  kQ::FT = 0.45 | (0.10, 1.0) | "-"

  watercons::AbstractWaterConsGPPModel{FT} = β_GPP_Zhang2019{FT}()

  PC_photo::Bool = false | (NaN, NaN) | "-" # 是否开启光周期, 不参与参数优化
end


"""
    photosynthesis(Tavg::T, Rs::T, VPD::T, LAI::T, Ca=380.0; par)

# Example
```julia
# Ag, Rd = photosynthesis(photo, Tavg, Rs, VPD, LAI, Ca, PC)
```
"""
function photosynthesis(
  photo::Photosynthesis_Rong2018{T},
  Tavg::T, Rs::T, VPD::T, LAI::T, Ca::T=380.0, PC::T=1.0) where {T<:Real}
  (; α, η, VCmax25, d_PC, kQ) = photo

  PAR = 0.45 * Rs # W m-2, taken as 0.45 time of solar radiation
  PAR_mol = PAR * 4.57 # 1 W m-2 = 4.57 umol m-2 s-1

  Vm = VCmax25 * T_adjust_Vm25(Tavg) * PC^d_PC # * data$dhour_norm^2 
  Am = Vm # 认为最大光合速率 = 最大羧化能力

  P1 = Am * α * η * PAR_mol
  P2 = Am * α * PAR_mol
  P3 = Am * η * Ca
  P4 = α * η * PAR_mol * Ca

  ## canopy conductance in (mol m-2 s-1)
  Ags = Ca * P1 / (P2 * kQ + P4 * kQ) * (
    kQ * LAI + log((P2 + P3 + P4) / (P2 + P3 * exp(kQ * LAI) + P4))) # [umol m-2 s-1]
  Ag = Ags * β_GPP(VPD, photo.watercons) # gross assimilation rate in [umol m-2 s-1]

  Rd = 0.15Vm # Collatz et al. (1991); CLM5 Tech Note
  Ag, Rd # [umol m-2 s-1]
end


"""
    photosynthesis_single_layer(photo, Tavg, Rs, VPD, LAI, Ca, PC; is_sunlit)

单层光合作用计算（用于向阳叶或背阴叶）

# 参数
- `photo`: 光合作用模型
- `Tavg`: 平均温度 [℃]
- `Rs`: 太阳辐射 [W m⁻²]
- `VPD`: 饱和水汽压差 [kPa]
- `LAI`: 叶面积指数 [m² m⁻²]
- `Ca`: 大气 CO2 浓度 [μmol mol⁻¹]
- `PC`: 光周期约束 [-]
- `is_sunlit`: 是否为向阳叶（Bool）

# 返回值
- `Ag`: 光合速率 [μmol m⁻² s⁻¹]
- `Rd`: 呼吸速率 [μmol m⁻² s⁻¹]
"""
function photosynthesis_single_layer(
  photo::Photosynthesis_Rong2018{T},
  Tavg::T, Rs::T, VPD::T, LAI::T, Ca::T, PC::T;
  is_sunlit::Bool=true) where {T<:Real}

  (; α, η, VCmax25, d_PC, kQ) = photo

  # 边界条件检查
  if Rs <= 0 || LAI <= 0
    # 无辐射或无叶片，返回零
    Vm = VCmax25 * T_adjust_Vm25(Tavg) * PC^d_PC
    Rd = T(0.15) * Vm
    return T(0.0), Rd
  end

  # PAR 计算
  PAR = T(0.45) * Rs  # W m-2
  PAR_mol = PAR * T(4.57)  # μmol m-2 s-1

  # 调整 PAR 强度：向阳叶=直射辐射，背阴叶=散射辐射
  # 散射辐射约占总 PAR 的 15%
  PAR_factor = is_sunlit ? one(T) : T(0.15)
  PAR_mol_adjusted = PAR_mol * PAR_factor

  # 温度调整后的 VCmax
  Vm = VCmax25 * T_adjust_Vm25(Tavg) * PC^d_PC
  Am = Vm  # 最大光合速率 = 最大羧化能力

  # P 参数
  P1 = Am * α * η * PAR_mol_adjusted
  P2 = Am * α * PAR_mol_adjusted
  P3 = Am * η * Ca
  P4 = α * η * PAR_mol_adjusted * Ca

  # 冠层导度 [μmol m⁻² s⁻¹]
  Ags = Ca * P1 / (P2 * kQ + P4 * kQ) * (
    kQ * LAI + log((P2 + P3 + P4) / (P2 + P3 * exp(kQ * LAI) + P4)))

  # 水分胁迫下的实际光合速率
  Ag = Ags * β_GPP(VPD, photo.watercons)

  # 呼吸速率
  Rd = T(0.15) * Vm

  return Ag, Rd
end


"""
    photosynthesis(photo, Tavg, Rs, VPD, canopy::TwoLeaf, Ca, PC)

TwoLeaf 双叶模型光合作用

分别计算向阳叶和背阴叶的光合速率，然后求和。

# 参数
- `photo`: 光合作用模型
- `Tavg`: 平均温度 [℃]
- `Rs`: 太阳辐射 [W m⁻²]
- `VPD`: 饱和水汽压差 [kPa]
- `canopy::TwoLeaf`: 双叶冠层结构（需要先调用 allocate_LAI! 分配 LAI）
- `Ca`: 大气 CO2 浓度 [μmol mol⁻¹]
- `PC`: 光周期约束 [-]

# 返回值
- `Ag_total`: 总光合速率 [μmol m⁻² s⁻¹]
- `Rd_total`: 总呼吸速率 [μmol m⁻² s⁻¹]

# 示例
```julia
using SPAC

photo = Photosynthesis_Rong2018{Float64}()
canopy = TwoLeaf{Float64}(Lai=5.0, Ω=1.0)

# 先分配 LAI
CosZs = 0.866
allocate_LAI!(canopy, CosZs)

# 计算光合作用
Ag, Rd = photosynthesis(photo, 25.0, 200.0, 2.0, canopy, 380.0, 1.0)
```
"""
function photosynthesis(
  photo::Photosynthesis_Rong2018{T},
  Tavg::T, Rs::T, VPD::T,
  canopy::TwoLeaf{T},
  Ca::T=T(380.0), PC::T=T(1.0)) where {T<:Real}

  (; Lai_sunlit, Lai_shaded) = canopy

  # 边界条件：零 LAI
  total_LAI = Lai_sunlit + Lai_shaded
  if total_LAI <= T(1e-6)
    return T(0.0), T(0.0)
  end

  Rs_sunlit, Rs_shaded = partition_sunshade_radiation_beer(Lai_sunlit, Lai_shaded, Rs)

  # 向阳叶光合作用
  if Lai_sunlit > T(1e-6)
    Ag_sunlit, Rd_sunlit = photosynthesis_single_layer(
      photo, Tavg, Rs_sunlit, VPD, Lai_sunlit, Ca, PC; is_sunlit=true)
  else
    Ag_sunlit = T(0.0)
    Rd_sunlit = T(0.0)
  end

  # 背阴叶光合作用
  if Lai_shaded > T(1e-6)
    Ag_shaded, Rd_shaded = photosynthesis_single_layer(
      photo, Tavg, Rs_shaded, VPD, Lai_shaded, Ca, PC; is_sunlit=false)
  else
    Ag_shaded = T(0.0)
    Rd_shaded = T(0.0)
  end

  # 总光合 = 向阳 + 背阴
  Ag_total = Ag_sunlit + Ag_shaded
  Rd_total = Rd_sunlit + Rd_shaded

  return Ag_total, Rd_total
end
