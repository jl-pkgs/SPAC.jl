"""
    CanopyInterface.jl

统一冠层模型接口，为不同冠层方案提供一致的方法调用。

这个模块定义了一套统一的接口函数，使得 `evapotranspiration!` 等核心函数
可以以统一的方式处理不同的冠层类型（BigLeaf, TwoLeaf, TwoBigLeaf, Leaves, OverUnderCanopy）。

## 接口函数

每个冠层类型需要实现以下函数：
- `get_lai(canopy)` - 获取总叶面积指数
- `allocate_lai!(canopy, CosZs)` - 分配向阳叶和背阴叶LAI（可选）
- `get_rn_canopy(canopy)` - 获取冠层净辐射
- `get_rn_soil(canopy)` - 获取土壤净辐射
- `requires_lai_allocation(canopy)` - 是否需要LAI分配
"""

export get_lai, allocate_lai_if_needed!
export get_rn_canopy, get_rn_soil
export requires_lai_allocation, canopy_type_name
export radiative_transfer_for_canopy!


"""
    get_lai(canopy::AbstractLeaf{T}) where {T}

获取冠层的总叶面积指数 [m² m⁻²]
"""
function get_lai(canopy::BigLeaf{T}) where {T}
  canopy.Lai
end

function get_lai(canopy::TwoLeaf{T}) where {T}
  canopy.Lai
end

function get_lai(canopy::TwoBigLeaf{T}) where {T}
  canopy.Lai_sunlit + canopy.Lai_shaded
end

function get_lai(canopy::Leaves{T}) where {T}
  sum(canopy.Lai_sunlit) + sum(canopy.Lai_shaded)
end

function get_lai(canopy::OverUnderCanopy{T}) where {T}
  canopy.LAI_over + canopy.LAI_under
end


"""
    requires_lai_allocation(canopy::AbstractLeaf{T}) where {T}

判断冠层类型是否需要进行向阳叶/背阴叶LAI分配。

- BigLeaf: 不需要
- TwoLeaf: 需要
- TwoBigLeaf: 需要
- Leaves: 需要（在多层辐射传输中处理）
- OverUnderCanopy: 不需要（双层结构，不需要sunlit/shaded区分）
"""
requires_lai_allocation(::BigLeaf) = false
requires_lai_allocation(::TwoLeaf) = true
requires_lai_allocation(::TwoBigLeaf) = true
requires_lai_allocation(::Leaves) = true
requires_lai_allocation(::OverUnderCanopy) = false


"""
    allocate_lai_if_needed!(canopy::AbstractLeaf{T}, CosZs::T) where {T}

根据冠层类型，如果需要则进行LAI分配。

对于 TwoLeaf 和 TwoBigLeaf，调用 `allocate_LAI!` 进行向阳叶/背阴叶分配。
对于其他类型，直接返回 canopy（无操作）。
"""
function allocate_lai_if_needed!(canopy::BigLeaf{T}, ::T) where {T}
  canopy
end

function allocate_lai_if_needed!(canopy::TwoLeaf{T}, CosZs::T) where {T}
  allocate_LAI!(canopy, CosZs)
end

function allocate_lai_if_needed!(canopy::TwoBigLeaf{T}, CosZs::T) where {T}
  allocate_LAI!(canopy, CosZs)
end

function allocate_lai_if_needed!(canopy::Leaves{T}, ::T) where {T}
  # Leaves 类型的LAI分配在多层辐射传输中处理
  canopy
end

function allocate_lai_if_needed!(canopy::OverUnderCanopy{T}, ::T) where {T}
  # OverUnderCanopy 不需要 sunlit/shaded 分配
  canopy
end


"""
    get_rn_canopy(canopy::AbstractLeaf{T}) where {T}

获取冠层的净辐射 [W m⁻²]
"""
function get_rn_canopy(canopy::Union{BigLeaf, TwoLeaf, TwoBigLeaf})
  canopy.Rn_c
end

function get_rn_canopy(canopy::Leaves)
  canopy.Rn[1]  # 第一层的净辐射
end

function get_rn_canopy(canopy::OverUnderCanopy)
  canopy.Rn_over
end


"""
    get_rn_soil(canopy::AbstractLeaf{T}) where {T}

获取土壤的净辐射 [W m⁻²]
"""
function get_rn_soil(canopy::Union{BigLeaf, TwoLeaf, TwoBigLeaf})
  canopy.Rn_s
end

function get_rn_soil(canopy::Leaves)
  canopy.Rn[end]  # 最后一个元素是地面
end

function get_rn_soil(canopy::OverUnderCanopy)
  canopy.Rn_soil
end


"""
    canopy_type_name(canopy::AbstractLeaf)

获取冠层类型的字符串名称，用于日志和调试。
"""
canopy_type_name(::BigLeaf) = "BigLeaf"
canopy_type_name(::TwoLeaf) = "TwoLeaf"
canopy_type_name(::TwoBigLeaf) = "TwoBigLeaf"
canopy_type_name(::Leaves) = "Leaves"
canopy_type_name(::OverUnderCanopy) = "OverUnderCanopy"


"""
    radiative_transfer_for_canopy!(canopy::AbstractLeaf{T}, Rn::T; kA::T=T(0.7)) where {T}

根据冠层类型选择合适的辐射传输函数。
"""
function radiative_transfer_for_canopy!(canopy::Union{BigLeaf, TwoLeaf, TwoBigLeaf}, Rn::T; kA::T=T(0.7)) where {T}
  radiative_transfer!(canopy, Rn; kA)
end

function radiative_transfer_for_canopy!(canopy::Leaves, Rn::T; kA::T=T(0.7)) where {T}
  # Leaves 类型使用多层辐射传输函数（计算边界辐射）
  radiative_transfer_multilayer!(canopy, Rn; kA)
end

function radiative_transfer_for_canopy!(canopy::OverUnderCanopy, Rn::T; kA::T=T(0.7)) where {T}
  # 双层模型使用专门的辐射传输函数
  radiative_transfer_2layer!(canopy, Rn; kA)
end
