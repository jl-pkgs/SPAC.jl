using Statistics

"""
    photosynthesis_multilayer!(leaves::Leaves{T}, photo::AbstractPhotosynthesisModel{T},
                               air::AirLayer{T}, VCmax_sunlit::Vector{T}, VCmax_shaded::Vector{T}) where {T<:Real}

计算多层冠层的光合作用，每层使用不同的 VCmax 值（垂向分布）
"""
function photosynthesis_multilayer!(
  leaves::Leaves{T},
  photo::AbstractPhotosynthesisModel{T},
  air::AirLayer{T},
  VCmax_sunlit::Vector{T},
  VCmax_shaded::Vector{T}
) where {T<:Real}

  (; Tavg, Rs, VPD, Ca, PC) = air
  (; Lai_sunlit, Lai_shaded, nlyr) = leaves

  # 计算太阳天顶角余弦（简化假设）
  CosZs = Rs > 0 ? T(0.866) : T(0.0)

  # 初始化总 GPP 和 Rd
  GPP_total = T(0.0)
  Rd_total = T(0.0)

  # 计算消光系数
  Ω = T(1.0)
  K = CosZs > 0 ? T(0.5) * Ω / CosZs : T(0.5)

  # PAR 转换系数
  PAR_factor = T(0.45)
  PAR_mol_factor = T(4.57)

  # 对每一层计算光合作用
  for i in 1:nlyr
    # 当前层顶部的累积 LAI
    dLAI = sum(Lai_sunlit) + sum(Lai_shaded)
    if dLAI > 0
      LAI_per_layer = dLAI / nlyr
      L_cumulative_top = (i - 1) * LAI_per_layer
      L_cumulative_mid = (i - 0.5) * LAI_per_layer
    else
      L_cumulative_top = T(0.0)
      L_cumulative_mid = T(0.0)
    end

    # PAR 在当前层的衰减
    PAR_transmission = exp(-K * L_cumulative_top)
    PAR_layer = Rs * PAR_factor * PAR_transmission
    PAR_mol_layer = PAR_layer * PAR_mol_factor

    # 向阳叶光合作用
    if Lai_sunlit[i] > T(1e-6) && VCmax_sunlit[i] > T(0.0)
      Ag_sunlit, Rd_sunlit = photosynthesis_single_layer_vcmax(
        photo, Tavg, PAR_mol_layer, VPD, Lai_sunlit[i], Ca, PC,
        VCmax_sunlit[i]; is_sunlit=true)
      GPP_sunlit = umol2gC(Ag_sunlit)
      Rd_sunlit_gC = umol2gC(Rd_sunlit)
    else
      GPP_sunlit = T(0.0)
      Rd_sunlit_gC = T(0.0)
    end

    # 背阴叶光合作用
    if Lai_shaded[i] > T(1e-6) && VCmax_shaded[i] > T(0.0)
      Ag_shaded, Rd_shaded = photosynthesis_single_layer_vcmax(
        photo, Tavg, PAR_mol_layer, VPD, Lai_shaded[i], Ca, PC,
        VCmax_shaded[i]; is_sunlit=false)
      GPP_shaded = umol2gC(Ag_shaded)
      Rd_shaded_gC = umol2gC(Rd_shaded)
    else
      GPP_shaded = T(0.0)
      Rd_shaded_gC = T(0.0)
    end

    # 更新 leaves 结构
    leaves.GPP_sunlit[i] = GPP_sunlit
    leaves.GPP_shaded[i] = GPP_shaded

    # 累加
    GPP_total += GPP_sunlit + GPP_shaded
    Rd_total += Rd_sunlit_gC + Rd_shaded_gC
  end

  return GPP_total, Rd_total
end


"""
    photosynthesis_single_layer_vcmax(photo, Tavg, PAR_mol, VPD, LAI, Ca, PC, VCmax25_layer; is_sunlit)

单层光合作用计算，使用指定的 VCmax 值
"""
function photosynthesis_single_layer_vcmax(
  photo::Photosynthesis_Rong2018{T},
  Tavg::T,
  PAR_mol::T,
  VPD::T,
  LAI::T,
  Ca::T,
  PC::T,
  VCmax25_layer::T;
  is_sunlit::Bool=true
) where {T<:Real}

  (; α, η, d_PC, kQ) = photo

  # 调整 PAR 强度
  PAR_factor = is_sunlit ? one(T) : T(0.15)
  PAR_mol_adjusted = PAR_mol * PAR_factor

  # 温度调整后的 VCmax
  Vm = VCmax25_layer * T_adjust_Vm25(Tavg) * PC^d_PC
  Am = Vm

  # P 参数
  P1 = Am * α * η * PAR_mol_adjusted
  P2 = Am * α * PAR_mol_adjusted
  P3 = Am * η * Ca
  P4 = α * η * PAR_mol_adjusted * Ca

  # 冠层导度
  if LAI > T(1e-6)
    Ags = Ca * P1 / (P2 * kQ + P4 * kQ) * (
      kQ * LAI + log((P2 + P3 + P4) / (P2 + P3 * exp(kQ * LAI) + P4)))
  else
    Ags = T(0.0)
  end

  # 水分胁迫下的实际光合速率
  Ag = Ags * β_GPP(VPD, photo.watercons)

  # 呼吸速率
  Rd = T(0.15) * Vm

  return Ag, Rd
end
