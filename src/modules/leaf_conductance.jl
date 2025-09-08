function leaf_conductance(
  air::AirLayer{T},
  canopy::AbstractLeaf{T},
  photo::AbstractPhotosynthesisModel{T},
  stomatal::AbstractStomatalModel{T}) where {T}

  (; Tavg, Rs, VPD, Pa, Ca, PC) = air
  (; Lai) = canopy
  Lai <= 0.01 && return T(0.0), T(0.0) # GPP, rs # 无冠层无rs

  Ag, Rd = photosynthesis(photo, Tavg, Rs, VPD, Lai, Ca, PC)
  GPP = umol2gC(Ag)

  gs = stomatal_conductance(stomatal, Ag, Rd, VPD, Ca) # [mol m-2 s-1], 0.1~0.4
  rs = 1 / mol2m(gs, Tavg, Pa) #[s m-1]
  GPP, rs
end
