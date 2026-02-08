export allocate_LAI!

"""
    allocate_LAI!(leaf::TwoLeaf{T}, CosZs::T) where {T<:Real}

分配总LAI到向阳叶和背阴叶（单层双叶模型）

基于 Norman (1982) 和 Chen et al. (2012) 的方法，使用 Beer 定律计算向阳叶面积。

# 参数
- `leaf::TwoLeaf{T}`: 双叶冠层结构，会更新 Lai_sunlit 和 Lai_shaded 字段
- `CosZs::T`: 太阳天顶角余弦值

# 算法
消光系数 K = 0.5 * Ω / CosZs （球形叶倾角分布）

向阳叶 LAI = (1 - exp(-K * LAI_total)) / K
背阴叶 LAI = LAI_total - 向阳叶 LAI

# 边界条件
- CosZs ≤ 0（夜间）：所有叶片为背阴叶
- LAI_total ≤ 0：返回零值

# 返回值
更新后的 `leaf` 对象（支持链式调用）

# 示例
```julia
using SPAC

canopy = TwoLeaf{Float64}(Lai=5.0, Ω=1.0)
CosZs = 0.866  # 太阳天顶角 30°
allocate_LAI!(canopy, CosZs)

println("向阳叶 LAI: ", canopy.Lai_sunlit)
println("背阴叶 LAI: ", canopy.Lai_shaded)
println("LAI 守恒: ", canopy.Lai_sunlit + canopy.Lai_shaded ≈ canopy.Lai)
```
"""
function allocate_LAI!(leaf::TwoLeaf{T}, CosZs::T) where {T<:Real}
  (; Lai, Ω) = leaf

  # 边界条件：零LAI
  if Lai <= 0
    leaf.Lai_sunlit = T(0.0)
    leaf.Lai_shaded = T(0.0)
    return leaf
  end

  # 边界条件：夜间（CosZs ≤ 0）
  if CosZs <= 0
    leaf.Lai_sunlit = T(0.0)
    leaf.Lai_shaded = Lai
    return leaf
  end

  # 消光系数（球形叶倾角分布）
  K = T(0.5) * Ω / CosZs

  # 向阳叶 LAI（解析解）
  # LAI_sunlit = ∫₀^LAI K * exp(-K * L) dL = 1 - exp(-K * LAI)
  # 但这里需要除以K得到实际的叶面积
  sunlit_LAI = (one(T) - exp(-K * Lai)) / K

  # 背阴叶 LAI
  shaded_LAI = Lai - sunlit_LAI

  # 确保非负
  leaf.Lai_sunlit = max(T(0.0), sunlit_LAI)
  leaf.Lai_shaded = max(T(0.0), shaded_LAI)

  return leaf
end


# TwoBigLeaf allocation is in twobigleaf_allocation.jl
