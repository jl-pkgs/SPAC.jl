@bounds @units @with_kw mutable struct Stomatal_Yu2004{FT} <: AbstractStomatalModel{FT}
  "stomatal conductance coefficient, [unitless]" # 气孔导度斜率参数
  g1::FT = 10.00 | (2.00, 100.00) | "-" # (2-20) in Rong2018

  "水汽压参数"  # leuning 2008, Yu 2004
  D0::FT = 0.7 | (0.50, 2.0) | "kPa"

  # PC_g0::Bool = false | (NaN, NaN) | "-" # ignored

  PC_g1::Bool = false | (NaN, NaN) | "-" # PC施加到斜率项
end


"""
- `g0`: 气孔导度截距项
  + `CLM5`      : 0.0001 [mol m-2 s-1] for C3 and C4
  + `Medlyn2011`: [-0.04, 0.04]
"""
@bounds @units @with_kw mutable struct Stomatal_Medlyn2011{FT} <: AbstractStomatalModel{FT}
  "stomatal conductance coefficient, `mol m⁻² s⁻¹`" # 气孔导度斜率参数
  g0::FT = 0.0001 | (0.0, 0.04) | "mol m⁻² s⁻¹" # 100 μmol m⁻² s⁻¹

  "stomatal conductance coefficient, `sqrt(kPa)`" # 气孔导度斜率参数
  g1::FT = 2.00 | (1.0, 6.0) | "sqrt(kPa)"

  PC_g0::Bool = false | (NaN, NaN) | "-" # PC施加到截距项

  PC_g1::Bool = false | (NaN, NaN) | "-" # PC施加到斜率项
end


# Ball-Berry 气孔导度模型
@bounds @units @with_kw mutable struct Stomatal_BallBerry1987{FT} <: AbstractStomatalModel{FT}
  "气孔导度截距项 [mol m⁻² s⁻¹]"
  g0::FT = 0.0001 | (0.0, 0.04) | "mol m⁻² s⁻¹"

  "气孔导度斜率参数 [unitless]"
  g1::FT = 9.0 | (0.0, 20.0) | "-"

  PC_g0::Bool = false | (NaN, NaN) | "-"
  PC_g1::Bool = false | (NaN, NaN) | "-"
end


"""
    stomatal_conductance(stomatal::Stomatal_Medlyn2011, Ag::T, Rd::T, D::T, Ca::T)
    stomatal_conductance(stomatal::Stomatal_Yu2004, Ag::T, Rd::T, D::T, Ca::T)

# Arguments
- `An`: net assimilation rate, [umol m-2 s-1], `An = Ag - Rd`
- `An`: net assimilation rate, [umol m-2 s-1]
- `Rd`: 0.15Vcmax, [umol m-2 s-1]
- `D` : air vapor pressure, [kPa]
- `Ca`: ~380 ppm, 380 [10^-6 mol mol-1]

> 1.6是`Gc` to `Gc_w`
Gs在0.1 ~ 0.4 [mol m-2 s-1]是合理的值域，对应的阻力在 100~400 [s m-1]

```julia
stomatal_conductance(stomatal, Ag, Rd, VPD, Ca)
```
"""
function stomatal_conductance(stomatal::Stomatal_Medlyn2011, Ag::T, Rd::T, D::T, Ca::T, PC::T, ignored...) where {T<:Real}
  (; g0, g1, PC_g0, PC_g1) = stomatal
  _g0::T = PC_g0 ? g0 * PC : g0
  _g1::T = PC_g1 ? g1 * PC : g1

  An::T = Ag - Rd
  gs::T = _g0 + 1.6(1 + _g1 / sqrt(D)) * (An / Ca) # second term in [mol m-2 s-1]
  return gs # [mol m-2 s-1]
end


function stomatal_conductance(stomatal::Stomatal_Yu2004, Ag::T, Rd::T, D::T, Ca::T, PC::T, ignored...) where {T<:Real}
  (; g1, D0, PC_g1) = stomatal
  _g1::T = PC_g1 ? g1 * PC : g1

  gs::T = 1.6_g1 * Ag / (Ca * (1 + D / D0)) # Yu 2004; Rong 2018 Eq. A3
  return gs # [mol m-2 s-1]
end

function stomatal_conductance(stomatal::Stomatal_BallBerry1987, Ag::T, Rd::T, D::T, Ca::T, PC::T, Ta::T=T(25.0)) where {T<:Real}
  (; g0, g1, PC_g0, PC_g1) = stomatal
  _g0::T = PC_g0 ? g0 * PC : g0
  _g1::T = PC_g1 ? g1 * PC : g1

  es = cal_es(Ta)
  ea = es - D
  RH = ea / es
  
  An::T = Ag - Rd # [umol m-2 s-1] 
  # Ca: umol mol-1
  # Gan Rong, 2018, A1; He 2017, JGR-b, Eq. 1
  gs::T = _g0 + _g1 * (An * RH / Ca) # [mol m-2 s-1]
  return gs  # [mol m-2 s-1]
end


export Stomatal_Yu2004, Stomatal_Medlyn2011, Stomatal_BallBerry1987
export stomatal_conductance
