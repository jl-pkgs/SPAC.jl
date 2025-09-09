export Photosynthesis_Rong2018


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
