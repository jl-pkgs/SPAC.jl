"""
    Evapotranspiration_PML{FT<:Real} <: AbstractEvapotranspirationModel{FT}

# Fields
$(TYPEDFIELDS)
"""
@bounds @with_kw mutable struct Evapotranspiration_PML{FT<:Real} <: AbstractEvapotranspirationModel{FT}
  "extinction coefficients for available energy"
  kA::FT = 0.70 | (0.50, 0.9)

  "Specific leaf storage, van Dijk, A.I.J.M, 2001, Eq2"
  S_sls::FT = 0.1 | (0.01, 1.0)

  "Canopy cover fraction related parameter"
  fER0::FT = 0.1 | (0.01, 0.5)

  "canopy height, `[m]`"
  hc::FT = 1.0 | (0.01, 20.0)
  # LAIref::FT = 4.0   | (1.0, 6.0)      # 
  # frame::Integer = 10.0  | (6.0, 14.0) # 8-day moving window
end



"""
    PMLV2 (Penman–Monteith–Leuning Version 2) Evapotranspiration model

# Arguments

# References
1. Gan Rong, 2018, Ecohydrology
2. Zhang Yongqiang, 2019, RSE
3. Kong Dongdong, 2019, ISPRS
"""
function evapotranspiration!(
  output::SpacOutput{T},
  evap::AbstractEvapotranspirationModel{T},
  photo::AbstractPhotosynthesisModel{T},
  stomatal::AbstractStomatalModel{T},
  air::AirLayer{T},
  canopy::AbstractLeaf{T}) where {T<:Real}

  (; Prcp, Tavg, Rn, VPD, U2, Pa, Rs) = air
  (; fER0, S_sls, kA, hc) = evap

  ## 水面蒸发
  λ::T = cal_lambda(Tavg) # [MJ kg-1]
  γ::T = cal_gamma(Tavg)
  Δ::T = cal_slope(Tavg)

  Eeq_water = Δ / (Δ + γ) * Rn |> x -> W2mm(x, λ)
  Evp_water = γ / (Δ + γ) * 6.43 * (1 + 0.536 * U2) * VPD / λ # [mm]
  ET_water = Eeq_water + Evp_water

  ## Intercepted Evaporation (Ei)
  # Get LAI from different canopy types
  Lai_total = if isa(canopy, BigLeaf)
    canopy.Lai
  elseif isa(canopy, TwoLeaf)
    canopy.Lai
  elseif isa(canopy, TwoBigLeaf)
    canopy.Lai_sunlit + canopy.Lai_shaded
  elseif isa(canopy, Leaves)
    sum(canopy.Lai_sunlit) + sum(canopy.Lai_shaded)
  elseif isa(canopy, OverUnderCanopy)
    canopy.LAI_over + canopy.LAI_under
  else
    T(0.0)
  end

  Ei = cal_Ei_Dijk2021(Prcp, Lai_total, fER0, S_sls)
  Pi = Prcp - Ei # 应扣除这一部分消耗的能量
  # Rn = Rn - MJ2W(λ * Pi)

  # 根据冠层类型进行不同的处理
  if isa(canopy, BigLeaf)
    # 大叶模型：直接进行辐射传输
    radiative_transfer!(canopy, Rn; kA)
    (; Rn_c, Rn_s) = canopy
    GPP, rs = leaf_conductance(air, canopy, photo, stomatal)
    Ec, Ecr, Eca, ra = ET0_Monteith65(Rn_c, Tavg, VPD, U2, Pa; rs, hc, z_wind=2.0)

  elseif isa(canopy, TwoLeaf) || isa(canopy, TwoBigLeaf)
    # 双叶模型：先分配 LAI，再进行辐射传输
    CosZs = calculate_CosZs(Rs, Tavg, Pa)
    allocate_LAI!(canopy, CosZs)
    radiative_transfer!(canopy, Rn; kA)
    (; Rn_c, Rn_s) = canopy
    GPP, rs = leaf_conductance(air, canopy, photo, stomatal)
    Ec, Ecr, Eca, ra = ET0_Monteith65(Rn_c, Tavg, VPD, U2, Pa; rs, hc, z_wind=2.0)

  elseif isa(canopy, Leaves)
    # 多层模型：初始化多层结构
    # Note: Leaves 类型需要专门的辐射传输处理
    # 这里暂时使用简化方法
    radiative_transfer!(canopy, Rn; kA)
    GPP, rs = leaf_conductance(air, canopy, photo, stomatal)
    # 对于多层模型，Rn_c 需要从冠层顶层获取
    Rn_c = canopy.Rn[1]  # 第一层的净辐射
    Rn_s = canopy.Rn[end]  # 土壤层的净辐射
    Ec, Ecr, Eca, ra = ET0_Monteith65(Rn_c, Tavg, VPD, U2, Pa; rs, hc, z_wind=2.0)

  elseif isa(canopy, OverUnderCanopy)
    # 双层模型：使用专门的辐射传输函数
    radiative_transfer_2layer!(canopy, Rn; kA)
    GPP, rs = leaf_conductance(air, canopy, photo, stomatal)
    # 对于双层模型，使用上层林冠的净辐射
    Ec, Ecr, Eca, ra = ET0_Monteith65(canopy.Rn_over, Tavg, VPD, U2, Pa; rs, hc, z_wind=2.0)
  end

  ## 需要根据前期 Pi/Es_eq 去计算土壤水限制β，因此无法立即给出Es
  # 根据冠层类型获取土壤净辐射
  Rn_s = if isa(canopy, BigLeaf) || isa(canopy, TwoLeaf) || isa(canopy, TwoBigLeaf)
    canopy.Rn_s
  elseif isa(canopy, Leaves)
    canopy.Rn[end]
  elseif isa(canopy, OverUnderCanopy)
    canopy.Rn_soil
  else
    T(0.0)
  end

  Es_eq = Δ / (Δ + γ) * Rn_s |> x -> W2mm(x, λ) # 土壤均衡蒸发
  @pack! output = GPP, Ec, Ecr, Eca, Ei, Pi, Es_eq, ET_water, rs, ra
  return output
end


"""
    PMLV2(Prcp, Tavg, Rs, Rn, VPD, U2, LAI, Pa, Ca; par=param0, frame=3)

# Notes
一个站点的计算。注意，不同植被类型，参数不同。

# Arguments
- `frame`: in 8-days
"""
function evapotranspiration(
  evap::AbstractEvapotranspirationModel{T},
  photo::AbstractPhotosynthesisModel{T},
  stomatal::AbstractStomatalModel{T},
  Lai::V,
  Prcp::V, Tavg::V, Rs::V, Rn::V,
  VPD::V, U2::V,
  Pa::V, Ca::V, PC::Union{T,V}=T(1.0);
  frame::Int=3,
  res::Union{Nothing,SpacOutputs}=nothing) where {T<:Real,V<:AbstractVector{T}}

  # 是否开启光周期
  (; PC_photo) = photo
  !isa(PC, Vector) && (PC_photo = false)

  canopy = BigLeaf{T}() # 先测试大叶模型
  air = AirLayer{T}()
  r = SpacOutput{T}()


  ntime = length(Prcp)
  res === nothing && (res = SpacOutputs{T}(; ntime))

  @inbounds for i = 1:ntime
    _PC = PC_photo ? PC[i] : T(1.0)
    update!(air, Prcp[i], Tavg[i], Rs[i], Rn[i], VPD[i], U2[i], Pa[i]; Ca=Ca[i], PC=_PC)
    canopy.Lai = Lai[i]

    evapotranspiration!(r, evap, photo, stomatal, air, canopy)
    res[i] = r
  end

  res.β_Es = movmean2(res.Pi, frame, 0) ./ movmean2(res.Es_eq, frame, 0)
  clamp!(res.β_Es, 0, 1)

  res.Es .= res.β_Es .* res.Es_eq
  res.ET .= res.Ec .+ res.Ei .+ res.Es
  res
end


function get_forcing(d::AbstractDataFrame, evap::Evapotranspiration_PML{T}) where {T}
  (; Prcp, Tavg, Rs, Rn, VPD, U2, LAI, Pa, Ca) = d
  PC = "PC" ∈ names(d) ? d.PC : T(1.0)
  (; LAI, Prcp, Tavg, Rs, Rn, VPD, U2, Pa, Ca, PC)
end
