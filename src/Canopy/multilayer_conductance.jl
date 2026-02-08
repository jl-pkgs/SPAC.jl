"""
    stomatal_conductance_multilayer!(leaves::Leaves{T}, photo::AbstractPhotosynthesisModel{T},
                                     stomatal::AbstractStomatalModel{T}, air::AirLayer{T},
                                     VCmax_sunlit::Vector{T}, VCmax_shaded::Vector{T}) where {T<:Real}

计算多层冠层的气孔导度
"""
function stomatal_conductance_multilayer!(
  leaves::Leaves{T},
  photo::AbstractPhotosynthesisModel{T},
  stomatal::AbstractStomatalModel{T},
  air::AirLayer{T},
  VCmax_sunlit::Vector{T},
  VCmax_shaded::Vector{T}
) where {T<:Real}

  (; Tavg, Rs, VPD, Pa, Ca, PC) = air
  (; Lai_sunlit, Lai_shaded, nlyr) = leaves

  # 计算太阳天顶角余弦
  CosZs = Rs > 0 ? T(0.866) : T(0.0)

  # 计算消光系数
  Ω = T(1.0)
  K = CosZs > 0 ? T(0.5) * Ω / CosZs : T(0.5)

  # PAR 转换系数
  PAR_factor = T(0.45)
  PAR_mol_factor = T(4.57)

  # 初始化
  gs_weighted_sum = T(0.0)
  total_LAI = T(0.0)
  GPP_total = T(0.0)

  # 对每一层计算气孔导度
  for i in 1:nlyr
    # 当前层顶部的累积 LAI
    dLAI_total = sum(Lai_sunlit) + sum(Lai_shaded)
    if dLAI_total > 0
      LAI_per_layer = dLAI_total / nlyr
      L_cumulative_top = (i - 1) * LAI_per_layer
    else
      L_cumulative_top = T(0.0)
    end

    # PAR 在当前层的衰减
    PAR_transmission = exp(-K * L_cumulative_top)
    PAR_layer = Rs * PAR_factor * PAR_transmission
    PAR_mol_layer = PAR_layer * PAR_mol_factor

    # 向阳叶气孔导度
    if Lai_sunlit[i] > T(1e-6) && VCmax_sunlit[i] > T(0.0)
      Ag_sunlit, Rd_sunlit = photosynthesis_single_layer_vcmax(
        photo, Tavg, PAR_mol_layer, VPD, Lai_sunlit[i], Ca, PC,
        VCmax_sunlit[i]; is_sunlit=true)
      gs_sunlit = stomatal_conductance(stomatal, Ag_sunlit, Rd_sunlit, VPD, Ca, PC, Tavg)
      GPP_sunlit = umol2gC(Ag_sunlit)
    else
      gs_sunlit = T(0.0)
      GPP_sunlit = T(0.0)
    end

    # 背阴叶气孔导度
    if Lai_shaded[i] > T(1e-6) && VCmax_shaded[i] > T(0.0)
      Ag_shaded, Rd_shaded = photosynthesis_single_layer_vcmax(
        photo, Tavg, PAR_mol_layer, VPD, Lai_shaded[i], Ca, PC,
        VCmax_shaded[i]; is_sunlit=false)
      gs_shaded = stomatal_conductance(stomatal, Ag_shaded, Rd_shaded, VPD, Ca, PC, Tavg)
      GPP_shaded = umol2gC(Ag_shaded)
    else
      gs_shaded = T(0.0)
      GPP_shaded = T(0.0)
    end

    # 更新 leaves 结构
    leaves.gs_sunlit[i] = gs_sunlit
    leaves.gs_shaded[i] = gs_shaded
    leaves.GPP_sunlit[i] = GPP_sunlit
    leaves.GPP_shaded[i] = GPP_shaded

    # LAI 加权平均
    layer_LAI = Lai_sunlit[i] + Lai_shaded[i]
    gs_weighted_sum += (gs_sunlit * Lai_sunlit[i] + gs_shaded * Lai_shaded[i])
    total_LAI += layer_LAI
    GPP_total += GPP_sunlit + GPP_shaded
  end

  # 计算冠层平均导度
  if total_LAI > T(1e-6)
    gs_canopy = gs_weighted_sum / total_LAI
  else
    gs_canopy = T(0.0)
  end

  # 转换为阻力
  rs = gs_canopy > T(0.0) ? T(1.0) / mol2m(gs_canopy, Tavg, Pa) : T(0.0)

  return GPP_total, rs
end
