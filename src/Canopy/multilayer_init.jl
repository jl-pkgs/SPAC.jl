"""
    initialize_multilayer!(leaves::Leaves{T}, LAI::T, Ω::T, CosZs::T, VCmax25::T, N_leaf::T, slope::T) where {T<:Real}

初始化多层冠层模型，分配 LAI 和 VCmax 到每一层
"""
function initialize_multilayer!(
  leaves::Leaves{T},
  LAI::T,
  Ω::T,
  CosZs::T,
  VCmax25::T,
  N_leaf::T,
  slope::T
) where {T<:Real}

  nlyr = leaves.nlyr

  # 创建 VCmax 输出数组
  VCmax_sunlit = zeros(T, nlyr)
  VCmax_shaded = zeros(T, nlyr)

  # 边界条件：零LAI
  if LAI <= 0
    fill!(leaves.Lai_sunlit, T(0.0))
    fill!(leaves.Lai_shaded, T(0.0))
    return leaves, VCmax_sunlit, VCmax_shaded
  end

  # 调用 VCmax_profile 计算垂向分布
  VCmax_sunlit_layer, VCmax_shaded_layer, LAI_sunlit_layer, LAI_shaded_layer =
    VCmax_profile(nlyr, LAI, Ω, CosZs, VCmax25, N_leaf, slope)

  # 存储到 leaves 结构中
  copyto!(leaves.Lai_sunlit, LAI_sunlit_layer)
  copyto!(leaves.Lai_shaded, LAI_shaded_layer)

  # 复制 VCmax 到输出数组
  copyto!(VCmax_sunlit, VCmax_sunlit_layer)
  copyto!(VCmax_shaded, VCmax_shaded_layer)

  return leaves, VCmax_sunlit, VCmax_shaded
end


"""
    initialize_multilayer_with_vcmax!(leaves::Leaves{T}, LAI::T, Ω::T, CosZs::T, VCmax25::T, N_leaf::T, slope::T,
                                     VCmax_sunlit::Vector{T}, VCmax_shaded::Vector{T}) where {T<:Real}

初始化多层冠层模型，并返回 VCmax 数组（带 VCmax 输出的版本）
"""
function initialize_multilayer_with_vcmax!(
  leaves::Leaves{T},
  LAI::T,
  Ω::T,
  CosZs::T,
  VCmax25::T,
  N_leaf::T,
  slope::T,
  VCmax_sunlit::Vector{T},
  VCmax_shaded::Vector{T}
) where {T<:Real}

  nlyr = leaves.nlyr

  # 边界条件：零LAI
  if LAI <= 0
    fill!(leaves.Lai_sunlit, T(0.0))
    fill!(leaves.Lai_shaded, T(0.0))
    fill!(VCmax_sunlit, T(0.0))
    fill!(VCmax_shaded, T(0.0))
    return leaves, VCmax_sunlit, VCmax_shaded
  end

  # 调用 VCmax_profile 计算垂向分布
  VCmax_sunlit_layer, VCmax_shaded_layer, LAI_sunlit_layer, LAI_shaded_layer =
    VCmax_profile(nlyr, LAI, Ω, CosZs, VCmax25, N_leaf, slope)

  # 存储到 leaves 结构中
  copyto!(leaves.Lai_sunlit, LAI_sunlit_layer)
  copyto!(leaves.Lai_shaded, LAI_shaded_layer)

  # 复制 VCmax 到输出数组
  copyto!(VCmax_sunlit, VCmax_sunlit_layer)
  copyto!(VCmax_shaded, VCmax_shaded_layer)

  return leaves, VCmax_sunlit, VCmax_shaded
end
