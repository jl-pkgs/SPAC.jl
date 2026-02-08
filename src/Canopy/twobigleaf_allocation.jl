export allocate_LAI!

"""
    allocate_LAI!(leaf::TwoBigLeaf{T}, CosZs::T) where {T<:Real}

分配总LAI到向阳叶和背阴叶（两大叶模型）

TwoBigLeaf 使用与 TwoLeaf 相同的 Beer 定律算法，但存储到 TwoBigLeaf 结构中。
与 TwoLeaf 的区别在于后续的导度计算策略：
- TwoLeaf: 导度按 LAI 加权平均（陈镜明 BEPS 方法）
- TwoBigLeaf: 导度积分到冠层（戴永久 CLM 方法）

# 参数
- `leaf::TwoBigLeaf{T}`: 两大叶冠层结构（包含 Lai 和 Ω 字段）
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

canopy = TwoBigLeaf{Float64}(Lai=5.0, Ω=1.0)
CosZs = 0.866  # 太阳天顶角 30°
allocate_LAI!(canopy, CosZs)

println("向阳叶 LAI: ", canopy.Lai_sunlit)
println("背阴叶 LAI: ", canopy.Lai_shaded)
println("LAI 守恒: ", canopy.Lai_sunlit + canopy.Lai_shaded ≈ canopy.Lai)
```
"""
function allocate_LAI!(leaf::TwoBigLeaf{T}, CosZs::T) where {T<:Real}
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
  sunlit_LAI = (one(T) - exp(-K * Lai)) / K

  # 背阴叶 LAI
  shaded_LAI = Lai - sunlit_LAI

  # 确保非负
  leaf.Lai_sunlit = max(T(0.0), sunlit_LAI)
  leaf.Lai_shaded = max(T(0.0), shaded_LAI)

  return leaf
end
