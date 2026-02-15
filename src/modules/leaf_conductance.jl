export leaf_conductance

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

  gs = stomatal_conductance(stomatal, Ag, Rd, VPD, Ca, PC, Tavg) # [mol m-2 s-1], 0.1~0.4
  rs = 1 / mol2m(gs, Tavg, Pa) #[s m-1]
  GPP, rs
end


"""
    leaf_conductance(air, canopy::TwoLeaf, photo, stomatal)

TwoLeaf 双叶模型叶片导度计算

分别计算向阳叶和背阴叶的导度，然后按 LAI 加权平均得到冠层导度。

# 算法
1. 分别计算向阳叶和背阴叶的光合作用（调用 photosynthesis_single_layer）
2. 分别计算向阳叶和背阴叶的气孔导度
3. 冠层平均导度 = (gs_sunlit*LAI_sunlit + gs_shaded*LAI_shaded) / total_LAI
4. 更新 canopy.GPP_sunlit 和 canopy.GPP_shaded 字段

# 返回值
- `GPP_total`: 总 GPP [gC m⁻² d⁻¹]
- `rs`: 冠层平均阻力 [s m⁻¹]
"""
function leaf_conductance(
  air::AirLayer{T},
  canopy::TwoLeaf{T},
  photo::AbstractPhotosynthesisModel{T},
  stomatal::AbstractStomatalModel{T}) where {T}

  (; Tavg, Rs, VPD, Pa, Ca, PC) = air
  (; Lai_sunlit, Lai_shaded) = canopy

  total_LAI = Lai_sunlit + Lai_shaded
  total_LAI <= T(0.01) && return T(0.0), T(0.0)

  Rs_sunlit, Rs_shaded = partition_sunshade_radiation_beer(Lai_sunlit, Lai_shaded, Rs)

  # 向阳叶导度
  if Lai_sunlit > T(0.01)
    Ag_sunlit, Rd_sunlit = photosynthesis_single_layer(
      photo, Tavg, Rs_sunlit, VPD, Lai_sunlit, Ca, PC; is_sunlit=true)
    gs_sunlit = stomatal_conductance(stomatal, Ag_sunlit, Rd_sunlit, VPD, Ca, PC, Tavg)
    GPP_sunlit = umol2gC(Ag_sunlit)
  else
    gs_sunlit = T(0.0)
    GPP_sunlit = T(0.0)
  end

  # 背阴叶导度
  if Lai_shaded > T(0.01)
    Ag_shaded, Rd_shaded = photosynthesis_single_layer(
      photo, Tavg, Rs_shaded, VPD, Lai_shaded, Ca, PC; is_sunlit=false)
    gs_shaded = stomatal_conductance(stomatal, Ag_shaded, Rd_shaded, VPD, Ca, PC, Tavg)
    GPP_shaded = umol2gC(Ag_shaded)
  else
    gs_shaded = T(0.0)
    GPP_shaded = T(0.0)
  end

  # 更新 canopy 字段（用于输出）
  canopy.GPP_sunlit = GPP_sunlit
  canopy.GPP_shaded = GPP_shaded

  # TwoLeaf 策略：导度不积分，按 LAI 加权平均（陈镜明 BEPS 方法）
  # 冠层平均导度 = (gs_sunlit*LAI_sunlit + gs_shaded*LAI_shaded) / total_LAI
  gs_canopy = (gs_sunlit * Lai_sunlit + gs_shaded * Lai_shaded) / total_LAI
  rs = T(1.0) / mol2m(gs_canopy, Tavg, Pa)

  GPP_total = GPP_sunlit + GPP_shaded
  return GPP_total, rs
end


"""
    leaf_conductance(air, canopy::TwoBigLeaf, photo, stomatal)

TwoBigLeaf 两大叶模型叶片导度计算（CLM 风格导度积分）

分别计算向阳叶和背阴叶的导度，然后积分到冠层尺度。

与 TwoLeaf 的关键区别：
- TwoLeaf: 导度按 LAI 加权平均，rs = 1 / mol2m(gs_canopy)
- TwoBigLeaf: 导度积分到冠层，蒸腾 T = T(Gs_sunlit) + T(Gs_shaded)

# 算法
1. 分别计算向阳叶和背阴叶的光合作用（调用 photosynthesis_single_layer）
2. 分别计算向阳叶和背阴叶的气孔导度 [mol m⁻² s⁻¹]
3. 将导度积分到冠层尺度：
   - Gs_sunlit_canopy = gs_sunlit * LAI_sunlit [mol m⁻² ground s⁻¹]
   - Gs_shaded_canopy = gs_shaded * LAI_shaded [mol m⁻² ground s⁻¹]
   - Gs_total = Gs_sunlit_canopy + Gs_shaded_canopy
4. 冠层阻力 rs = 1 / mol2m(Gs_total)
5. 更新 canopy.GPP_sunlit 和 canopy.GPP_shaded 字段

# 返回值
- `GPP_total`: 总 GPP [gC m⁻² d⁻¹]
- `rs`: 冠层阻力 [s m⁻¹]

# 参考文献
- CLM5 Tech Note: 导度积分方法
- Dai et al. (2004): 双大叶模型
"""
function leaf_conductance(
  air::AirLayer{T},
  canopy::TwoBigLeaf{T},
  photo::AbstractPhotosynthesisModel{T},
  stomatal::AbstractStomatalModel{T}) where {T}

  (; Tavg, Rs, VPD, Pa, Ca, PC) = air
  (; Lai_sunlit, Lai_shaded) = canopy

  total_LAI = Lai_sunlit + Lai_shaded
  total_LAI <= T(0.01) && return T(0.0), T(0.0)

  Rs_sunlit, Rs_shaded = partition_sunshade_radiation_beer(Lai_sunlit, Lai_shaded, Rs)

  # 向阳叶导度
  if Lai_sunlit > T(0.01)
    Ag_sunlit, Rd_sunlit = photosynthesis_single_layer(
      photo, Tavg, Rs_sunlit, VPD, Lai_sunlit, Ca, PC; is_sunlit=true)
    gs_sunlit = stomatal_conductance(stomatal, Ag_sunlit, Rd_sunlit, VPD, Ca, PC, Tavg)
    GPP_sunlit = umol2gC(Ag_sunlit)
  else
    gs_sunlit = T(0.0)
    GPP_sunlit = T(0.0)
  end

  # 背阴叶导度
  if Lai_shaded > T(0.01)
    Ag_shaded, Rd_shaded = photosynthesis_single_layer(
      photo, Tavg, Rs_shaded, VPD, Lai_shaded, Ca, PC; is_sunlit=false)
    gs_shaded = stomatal_conductance(stomatal, Ag_shaded, Rd_shaded, VPD, Ca, PC, Tavg)
    GPP_shaded = umol2gC(Ag_shaded)
  else
    gs_shaded = T(0.0)
    GPP_shaded = T(0.0)
  end

  # 更新 canopy 字段（用于输出）
  canopy.GPP_sunlit = GPP_sunlit
  canopy.GPP_shaded = GPP_shaded
  canopy.gs_sunlit = gs_sunlit
  canopy.gs_shaded = gs_shaded

  # TwoBigLeaf 策略：导度积分到冠层（CLM 风格）
  # Gs = gs_sunlit * LAI_sunlit + gs_shaded * LAI_shaded [mol m-2 ground s-1]
  # 蒸腾 T = T(Gs_sunlit) + T(Gs_shaded)
  Gs_canopy = gs_sunlit * Lai_sunlit + gs_shaded * Lai_shaded
  rs = T(1.0) / mol2m(Gs_canopy, Tavg, Pa)

  GPP_total = GPP_sunlit + GPP_shaded
  return GPP_total, rs
end


"""
    leaf_conductance(air, canopy::Leaves, photo, stomatal)

Leaves 多层模型叶片导度计算

对每一层分别计算向阳叶和背阴叶的导度，然后积分到冠层尺度。

# 算法
1. 对每一层 i：
   - 计算向阳叶和背阴叶的光合作用
   - 计算向阳叶和背阴叶的气孔导度
   - 累加 GPP：GPP_total += GPP_sunlit[i] + GPP_shaded[i]
2. 导度积分到冠层：
   - Gs_canopy = Σ(gs_sunlit[i] * LAI_sunlit[i] + gs_shaded[i] * LAI_shaded[i])
3. 冠层阻力 rs = 1 / mol2m(Gs_canopy)

# 返回值
- `GPP_total`: 总 GPP [gC m⁻² d⁻¹]
- `rs`: 冠层平均阻力 [s m⁻¹]
"""
function leaf_conductance(
  air::AirLayer{T},
  canopy::Leaves{T},
  photo::AbstractPhotosynthesisModel{T},
  stomatal::AbstractStomatalModel{T}) where {T}

  (; Tavg, Rs, VPD, Pa, Ca, PC) = air
  (; Lai_sunlit, Lai_shaded, nlyr) = canopy

  total_LAI = sum(Lai_sunlit) + sum(Lai_shaded)
  total_LAI <= T(0.01) && return T(0.0), T(0.0)

  GPP_total = T(0.0)
  Gs_canopy = T(0.0)

  # 对每一层计算导度
  for i in 1:nlyr
    # 向阳叶导度
    if Lai_sunlit[i] > T(0.01)
      Ag_sunlit, Rd_sunlit = photosynthesis_single_layer(
        photo, Tavg, Rs, VPD, Lai_sunlit[i], Ca, PC; is_sunlit=true)
      gs_sunlit = stomatal_conductance(stomatal, Ag_sunlit, Rd_sunlit, VPD, Ca, PC, Tavg)
      GPP_sunlit = umol2gC(Ag_sunlit)
      
      # 累加到冠层
      GPP_total += GPP_sunlit
      Gs_canopy += gs_sunlit * Lai_sunlit[i]
      
      # 更新 canopy 字段
      canopy.GPP_sunlit[i] = GPP_sunlit
      canopy.gs_sunlit[i] = gs_sunlit
    else
      canopy.GPP_sunlit[i] = T(0.0)
      canopy.gs_sunlit[i] = T(0.0)
    end

    # 背阴叶导度
    if Lai_shaded[i] > T(0.01)
      Ag_shaded, Rd_shaded = photosynthesis_single_layer(
        photo, Tavg, Rs, VPD, Lai_shaded[i], Ca, PC; is_sunlit=false)
      gs_shaded = stomatal_conductance(stomatal, Ag_shaded, Rd_shaded, VPD, Ca, PC, Tavg)
      GPP_shaded = umol2gC(Ag_shaded)
      
      # 累加到冠层
      GPP_total += GPP_shaded
      Gs_canopy += gs_shaded * Lai_shaded[i]
      
      # 更新 canopy 字段
      canopy.GPP_shaded[i] = GPP_shaded
      canopy.gs_shaded[i] = gs_shaded
    else
      canopy.GPP_shaded[i] = T(0.0)
      canopy.gs_shaded[i] = T(0.0)
    end
  end

  # 计算冠层阻力
  rs = T(1.0) / mol2m(Gs_canopy, Tavg, Pa)

  return GPP_total, rs
end
