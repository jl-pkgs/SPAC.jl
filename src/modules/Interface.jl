export LandModel, OverUnderCanopy
export AbstractWaterConsGPPModel, AbstractStomatalModel, AbstractPhotosynthesisModel, AbstractEvapotranspirationModel
export AirLayer, update!


using Parameters
import FieldMetadata: @bounds, bounds, @units, units

# abstract type AbstractModel{T} end

# 水分限制-GPP
abstract type AbstractWaterConsGPPModel{T} <: AbstractModel{T} end

# 气孔导度
abstract type AbstractStomatalModel{T} <: AbstractModel{T} end

# 光合
abstract type AbstractPhotosynthesisModel{T} <: AbstractModel{T} end

# 蒸发
abstract type AbstractEvapotranspirationModel{T} <: AbstractModel{T} end



@units @with_kw mutable struct AirLayer{T<:Real}
  Prcp::T = 0.0 | "mm d⁻¹"
  Tavg::T = 20.0 | "℃"
  Rs::T = 200.0 | "W m⁻²"
  Rn::T = 150 | "W m⁻²"
  VPD::T = 2.0 | "kPa"
  U2::T = 1.0 | "m s⁻¹"
  Pa::T = atm | "kPa"
  Ca::T = 380.0 | "μmol mol⁻¹" # air CO2 content
  PC::T = 1.0 | "-" # photoperiod constaint

  # auxiliary vairables
  λ::T = cal_lambda(Tavg) | "MJ kg⁻¹"
  γ::T = cal_gamma(Tavg, Pa) | "kPa K⁻¹"
  Δ::T = cal_slope(Tavg) | "kPa K⁻¹"
  β::T = γ / Δ | "-" # Bowen ratio, H / LE, γ / Δ
  ρₐ::T = cal_rho_a(Tavg, Pa) | "kg m⁻³"
end


function update!(air::AirLayer{T}; Tavg::T, Pa::T, kw...) where {T<:Real}
  λ::T = cal_lambda(Tavg) # "MJ kg⁻¹"
  γ::T = cal_gamma(Tavg, Pa) # "kPa K⁻¹"
  Δ::T = cal_slope(Tavg) # "kPa K⁻¹"
  β::T = γ / Δ # "-" # Bowen ratio, H / LE, γ / Δ
  ρₐ::T = cal_rho_a(Tavg, Pa) # "kg m⁻³"
  @pack! air = Tavg, Pa, λ, γ, Δ, β, ρₐ, kw...
end


function update!(air::AirLayer{T}, Prcp::T, Tavg::T, Rs::T, Rn::T, VPD::T, U2::T, Pa::T;
  Ca::T=T(380.0), PC=T(1.0)) where {T<:Real}

  λ::T = cal_lambda(Tavg) # "MJ kg⁻¹"
  γ::T = cal_gamma(Tavg, Pa) # "kPa K⁻¹"
  Δ::T = cal_slope(Tavg) # "kPa K⁻¹"
  β::T = γ / Δ # "-" # Bowen ratio, H / LE, γ / Δ
  ρₐ::T = cal_rho_a(Tavg, Pa) # "kg m⁻³"
  @pack! air = Tavg, Pa, λ, γ, Δ, β, ρₐ, Prcp, Rs, Rn, VPD, U2, Ca, PC
end


@with_kw_noshow mutable struct LandModel{T} <: AbstractModel{T}
  # atmospheric_forcing::AtmosphericForcing
  evap::AbstractEvapotranspirationModel{T}
  photo::AbstractPhotosynthesisModel{T}
  stomatal::AbstractStomatalModel{T}
end


function Base.show(io::IO, x::LandModel{T}) where {T<:Real}
  (; evap, photo, stomatal) = x
  
  printstyled(io, "LandModel{$T}: \n", color=:green, bold=true)

  printstyled(io, "photo{$T}: ", color=:blue, bold=true)
  print(io, photo)
  
  printstyled(io, "stomatal{$T}: ", color=:blue, bold=true)
  print(io, evap)

  printstyled(io, "stomatal{$T}: ", color=:blue, bold=true)
  print(io, stomatal)
  return nothing
end
