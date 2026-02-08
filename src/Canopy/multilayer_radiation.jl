"""
    radiative_transfer_multilayer!(leaves::Leaves{T}, Rn::T; kA::T=T(0.7)) where {T<:Real}

多层冠层的辐射传输计算（简化方案，使用 Beer 定律）
"""
function radiative_transfer_multilayer!(leaves::Leaves{T}, Rn::T; kA::T=T(0.7)) where {T<:Real}
  (; Lai_sunlit, Lai_shaded, nlyr) = leaves

  # 计算总 LAI 和每层 LAI
  total_LAI = sum(Lai_sunlit) + sum(Lai_shaded)

  # 边界条件：零 LAI
  if total_LAI <= T(1e-6)
    fill!(leaves.Rn, T(0.0))
    leaves.Rn[nlyr + 1] = Rn
    return leaves
  end

  # 计算每层的 LAI
  LAI_per_layer = total_LAI / nlyr

  # 计算每层边界的净辐射
  for i in 1:(nlyr + 1)
    L_cumulative = (i - 1) * LAI_per_layer
    leaves.Rn[i] = Rn * exp(-kA * L_cumulative)
  end

  return leaves
end


"""
    radiative_transfer_multilayer_absorbed!(leaves::Leaves{T}, Rn::T; kA::T=T(0.7)) where {T<:Real}

多层冠层的辐射传输计算（返回每层吸收的辐射）
"""
function radiative_transfer_multilayer_absorbed!(leaves::Leaves{T}, Rn::T; kA::T=T(0.7)) where {T<:Real}
  (; Lai_sunlit, Lai_shaded, nlyr) = leaves

  # 计算总 LAI 和每层 LAI
  total_LAI = sum(Lai_sunlit) + sum(Lai_shaded)

  # 边界条件：零 LAI
  if total_LAI <= T(1e-6)
    fill!(leaves.Rn, T(0.0))
    leaves.Rn[nlyr + 1] = Rn
    return leaves
  end

  # 计算每层的 LAI
  LAI_per_layer = total_LAI / nlyr

  # 先计算每层边界的净辐射
  Rn_boundary = zeros(T, nlyr + 1)
  for i in 1:(nlyr + 1)
    L_cumulative = (i - 1) * LAI_per_layer
    Rn_boundary[i] = Rn * exp(-kA * L_cumulative)
  end

  # 计算每层吸收的辐射
  for i in 1:nlyr
    leaves.Rn[i] = Rn_boundary[i] - Rn_boundary[i + 1]
  end

  # 最后一项是到达地面的辐射
  leaves.Rn[nlyr + 1] = Rn_boundary[nlyr + 1]

  return leaves
end
