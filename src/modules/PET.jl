# const Cp = 4.2 * 0.242 / 1000 # MJ kg-1 K-1
const Cp = 1.013 * 1e-3 # [MJ kg-1 K-1]
const ϵ = 0.622 # ratio of molecular weight of water vapor to dry air
const atm = 101.325 # kPa
const K0 = 273.15


"""
- `Rn`  : [W m-2]
- `Tair`: [℃]
- `Pa`  : [kPa]
"""
function ET0_eq(Rn::T, Tair::T, Pa::T=atm, ignored...) where {T<:Real}
  λ::T = cal_lambda(Tair)     # [MJ kg-1]
  Δ::T = cal_slope(Tair)      # [kPa degC-1]
  γ::T = Cp * Pa / (ϵ * λ)    # [kPa degC-1], Psychrometric constant

  Eeq::T = Δ / (Δ + γ) * Rn |> x -> W2mm(x, λ)
  Eeq = max(Eeq, 0)
  Eeq, λ, Δ, γ
end

# - Rn: W m-2
# - α: PT coefficient for water saturated surface
function ET0_PT72(Rn::T, Tair::T, Pa::T=atm; α=1.26) where {T<:Real}
  Eeq = ET0_eq(Rn, Tair, Pa)[1]
  α * Eeq # [mm d-1]
end


"""
ET0_Penman48(Rn, Tair, VPD, Uz)
"""
function ET0_Penman48(Rn::T, Tair::T, VPD::T, Uz::T, Pa::T=atm; z_wind=2) where {T<:Real}
  Eeq, λ, Δ, γ = ET0_eq(Rn, Tair, Pa)
  U2::T = cal_U2(Uz, z_wind)
  Evp::T = γ / (Δ + γ) * 6.43 * (1 + 0.536 * U2) * VPD / λ
  ET0::T = Eeq + Evp
  ET0
end


"""
    ET0_FAO98(Rn::T, Tair::T, VPD::T, wind::T, Pa::T=atm; z_wind=2, tall_crop=false)
# Examples
```julia
ET0_Penman48(200., 20., 2., 2.)
PET = ET0_FAO98(200.0, 20.0, 2.0, 2.0) # mm
PET = ET0_FAO98(Rn, Tair, VPD, U2, Pa; tall_crop = false) # mm
```
"""
function ET0_FAO98(Rn::T, Tair::T, VPD::T, wind::T, Pa::T=atm; z_wind=2, tall_crop=false) where {T<:Real}
  # T = eltype(Rn)
  U2 = cal_U2(wind, z_wind)

  if tall_crop
    p1 = T(1600.0)
    p2 = T(0.38)
  else
    p1 = T(900.0)
    p2 = T(0.34)
  end
  Eeq, λ, Δ, γ = ET0_eq(Rn, Tair, Pa)
  Eeq::T = Δ / (Δ + (γ * (1.0 + p2 * U2))) * Rn |> x -> W2mm(x, λ)
  Evp::T = γ * p1 / (Tair + 273.15) * U2 * VPD / (Δ + (γ * (1.0 + p2 * U2)))
  ET = Eeq + Evp
  (; ET, Eeq, Evp)
end



# β: 土壤水分限制因子, [0~1], 1充分供水。
"""
ET, Eeq, Evp, ra = ET0_Monteith65(Rn, Tair, VPD, Uz, Pa; rs, hc, z_wind=2.0)
"""
function ET0_Monteith65(Rn::T, Tair::T, VPD::T, Uz::T, Pa::T=atm;
  β::T=1.0, z_wind::T=2.0, rs::T=70.0, hc::T=0.12, z_obs::T=T(15.0)) where {T<:Real}
  Φ = T(0)
  rs <= 0.01 && return (; Ec=Φ, Eeq=Φ, Evp=Φ, ra=Φ) # 无冠层 无冠层蒸发

  Eeq, λ, Δ, γ = ET0_eq(Rn, Tair, Pa)
  U2 = cal_U2(Uz, z_wind)
  ra = aerodynamic_resistance(U2, hc; z_obs)
  rs = rs / β

  ρₐ = cal_rho_a(Tair, Pa) # FAO56, Eq. 3-5, [kg m-3]
  # Cp = 1.013 * 1e-3 # MJ kg-1 degC-1
  # ρₐ * Cp * dT * gH (in MJ m-2 s-1)
  # = kg m-3 * MJ kg-1 degC-1 * degC * m s-1
  # = MJ m-2 s-1
  Eeq = Δ / (Δ + (γ * (1 + rs / ra))) * Rn |> x -> W2mm(x, λ)
  Evp = (ρₐ * Cp * VPD / ra) / (Δ + (γ * (1 + rs / ra))) * 86400 / λ
  Ec = Eeq + Evp
  (; Ec, Eeq, Evp, ra)
end


# 空气动力学阻力
function aerodynamic_resistance(U2::T, hc::T=0.12; z_obs::T=T(2.0)) where {T<:Real}
  k = 0.41           # von Karman's constant, Allen 1998, Eq. 4
  d = 2 / 3 * hc     # zero plane displacement height
  z_om = 0.123 * hc  # roughness length governing momentum transfer 
  z_oh = 0.1 * z_om  # roughness length governing transfer of heat and vapour
  # z_m = 2.0        # height of wind measurements 
  # z_h = 2.0        # height of humidity measurements
  Uz = cal_Uz(U2, z_obs)
  z_m = z_h = z_obs

  ra = safe_log((z_m - d) / z_om) * safe_log((z_h - d) / z_oh) / (k^2 * Uz)
  isnan(ra) && return 208.0 / U2 # 70s/m, default 
  ra
end

safe_log(x::T) where {T<:Real} = x < 0 ? T(NaN) : log(x)


# λ: latent heat of vaporization (MJ kg-1)
W2mm(w::T, λ::T) where {T<:Real} = w * 0.086400 / λ

MJ2W(x::T) where {T<:Real} = x / 0.086400  # x * 1e6 / 86400

# λ: latent heat of vaporization (MJ kg-1)
cal_lambda(Ta::T) where {T<:Real} = 2.501 - 2.361 / 1000 * Ta # FAO65, Eq. 3-1


# γ: "kPa / K"
function cal_gamma(Tair::T, Pa::T=atm) where {T<:Real}
  λ = cal_lambda(Tair)
  Cp * Pa / (ϵ * λ)
end

function cal_bowen(Tair::T, Pa::T=atm) where {T<:Real}
  Δ = cal_slope(Tair)
  γ = cal_gamma(Tair, Pa)
  γ / Δ
end

# cal_gamma(Pa, λ) = (Cp * Pa) / (ϵ * λ)
cal_es(Ta::T) where {T<:Real} = 0.6108 * exp((17.27 * Ta) / (Ta + 237.3))

# "kPa K-1"
cal_slope(Ta::T) where {T<:Real} = 4098 * cal_es(Ta) / ((Ta + 237.3)^2)

function cal_Uz(U2::T, z)::T where {T<:Real}
  z == 2 && (return U2)
  log(67.8 * z - 5.42) / 4.87 * U2
end

function cal_U2(Uz::T, z=10.0) where {T<:Real}
  z == 2 && (return Uz)
  Uz * 4.87 / log(67.8 * z - 5.42)
end

# kg m-3
cal_rho_a(Tair, Pa) = 3.486 * Pa / cal_TvK(Tair)

# # rho_a: kg m-3
# function cal_rho_a(Tair, q, Pa)
#   3.486 * Pa / cal_TvK(Tair, q) # FAO56, Eq. 3-5
# end

## 几种计算虚温的方法

# FAO56, Eq. 3-7
cal_TvK(Tair) = 1.01 * (Tair + 273)

# 这个是最精确的版本
# FAO56, Eq. 3-6
cal_TvK(Tair, ea, Pa) = (Tair + K0) * (1 + (1 - epsilon) * ea / Pa)

# https://github.com/CUG-hydro/class2022_CUG_HydroMet/blob/master/ch02_大气的基本特征.md
# q ≈ ϵ*ea/Pa
# q = ϵ*ea/(Pa - (1 - ϵ)*ea)
function cal_TvK(Tair, q)
  # ea / Pa = q/(ϵ + (1 - ϵ) * q)
  (Tair + K0) * (1 + (1 - epsilon) * q / (ϵ + (1 - ϵ) * q))
end


export aerodynamic_resistance
export cal_rho_a
export ET0_eq, ET0_Penman48, ET0_Monteith65, ET0_PT72, ET0_FAO98
export cal_es
