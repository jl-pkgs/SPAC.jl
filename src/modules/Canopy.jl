export AbstractLeaf
export BigLeaf, TwoLeaf, TwoBigLeaf, Leaves


abstract type AbstractLeaf{FT<:AbstractFloat} <: AbstractModel{FT} end

"""
    struct BigLeaf{FT<:AbstractFloat} <: AbstractLeaf{FT}  

G_s = gs_sunlit * LAI_sunlit + gs_shaded * L_shaded
T = T(G_s)
"""
@bounds @units @with_kw mutable struct BigLeaf{FT} <: AbstractLeaf{FT}
  "Leaf area index [m² m⁻²]"
  Lai::FT = 2.0 | () | "m² m⁻²"

  "Stomatal conductance for H2O [m s⁻¹]"
  gs::FT = 0.0 | () | "m s⁻¹"

  # Rs_c::FT = 0.0
  # Rs_s::FT = 0.0
  "Canopy net radiation [W m⁻²]"
  Rn_c::FT = 0.0 | () | "W m⁻²"

  "Soil net radiation [W m⁻²]"
  Rn_s::FT = 0.0 | () | "W m⁻²"

  # H_s::FT = 0.0   # _s: soil
  # LE_s::FT = 0.0
  # H_c::FT = 0.0   # _c: canopy
  "Transpiration [mm d⁻¹]"
  T::FT = 0.0 | () | "mm d⁻¹"

  "Gross Primary Productivity [gC m⁻² d⁻¹]"
  GPP::FT = 0.0 | () | "gC m⁻² d⁻¹"
end

function radiative_transfer!(leaf::BigLeaf{T}, Rn::T; kA::T=T(0.7)) where {T}
  (; Lai) = leaf
  τ = Lai <= 0.01 ? 1.0 : exp(-kA * Lai) # 穿过的比例
  leaf.Rn_c = (1 - τ) * Rn
  leaf.Rn_s = τ * Rn
end


# 导度不进行积分，直接乘到冠层，陈镜明BEPS
"""
T = T(gs_sunlit) * LAI_sunlit + T(gs_shaded) * Lai_shaded
"""
@bounds @units @with_kw mutable struct TwoLeaf{FT} <: AbstractLeaf{FT}
  "Leaf area index [m² m⁻²]"
  Lai::FT = 2.0 | () | "m² m⁻²"

  "Clumping index [unitless]"
  Ω::FT = 1.0 | (0.1, 1.0) | "-"

  "Sunlit leaf area index [unitless]"
  Lai_sunlit::FT = 0.0 | () | "-"

  "Shaded leaf area index [unitless]"
  Lai_shaded::FT = 0.0 | () | "-"

  "Sunlit transpiration [mm d⁻¹]"
  T_sunlit::FT = 0.0 | () | "mm d⁻¹"

  "Shaded transpiration [mm d⁻¹]"
  T_shaded::FT = 0.0 | () | "mm d⁻¹"

  "Sunlit stomatal conductance [m s⁻¹]"
  gs_sunlit::FT = 0.0 | () | "m s⁻¹"

  "Shaded stomatal conductance [m s⁻¹]"
  gs_shaded::FT = 0.0 | () | "m s⁻¹"

  "Sunlit GPP [gC m⁻² d⁻¹]"
  GPP_sunlit::FT = 0.0 | () | "gC m⁻² d⁻¹"

  "Shaded GPP [gC m⁻² d⁻¹]"
  GPP_shaded::FT = 0.0 | () | "gC m⁻² d⁻¹"

  "Canopy net radiation [W m⁻²]"
  Rn_c::FT = 0.0 | () | "W m⁻²"

  "Soil net radiation [W m⁻²]"
  Rn_s::FT = 0.0 | () | "W m⁻²"
end

# 导度积分到冠层，戴永久CLM
"""
T = T(Gs_sunlit) + T(Gs_shaded)
"""
@bounds @units @with_kw mutable struct TwoBigLeaf{FT} <: AbstractLeaf{FT}
  "Total leaf area index [m² m⁻²]"
  Lai::FT = 2.0 | () | "m² m⁻²"

  "Clumping index [unitless]"
  Ω::FT = 1.0 | (0.1, 1.0) | "-"

  "Sunlit leaf area index [unitless]"
  Lai_sunlit::FT = 0.0 | () | "-"

  "Shaded leaf area index [unitless]"
  Lai_shaded::FT = 0.0 | () | "-"

  "Sunlit stomatal conductance [m s⁻¹]"
  gs_sunlit::FT = 0.0 | () | "m s⁻¹"

  "Shaded stomatal conductance [m s⁻¹]"
  gs_shaded::FT = 0.0 | () | "m s⁻¹"

  "Sunlit transpiration [mm d⁻¹]"
  T_sunlit::FT = 0.0 | () | "mm d⁻¹"

  "Shaded transpiration [mm d⁻¹]"
  T_shaded::FT = 0.0 | () | "mm d⁻¹"

  "Sunlit GPP [gC m⁻² d⁻¹]"
  GPP_sunlit::FT = 0.0 | () | "gC m⁻² d⁻¹"

  "Shaded GPP [gC m⁻² d⁻¹]"
  GPP_shaded::FT = 0.0 | () | "gC m⁻² d⁻¹"

  "Canopy net radiation [W m⁻²]"
  Rn_c::FT = 0.0 | () | "W m⁻²"

  "Soil net radiation [W m⁻²]"
  Rn_s::FT = 0.0 | () | "W m⁻²"
end


"""
    radiative_transfer!(leaf::TwoLeaf{T}, Rn::T; kA::T=T(0.7)) where {T}

双叶模型辐射传输（使用总LAI进行简单的Beer定律消光）

与 BigLeaf 相同的算法，因为双叶模型主要用于光合作用的区分，
辐射传输仍然使用总LAI计算。
"""
function radiative_transfer!(leaf::TwoLeaf{T}, Rn::T; kA::T=T(0.7)) where {T}
  (; Lai) = leaf
  τ = Lai <= 0.01 ? 1.0 : exp(-kA * Lai) # 穿过的比例
  leaf.Rn_c = (1 - τ) * Rn
  leaf.Rn_s = τ * Rn
end


"""
    radiative_transfer!(leaf::TwoBigLeaf{T}, Rn::T; kA::T=T(0.7)) where {T}

两大叶模型辐射传输（使用总LAI进行简单的Beer定律消光）

与 BigLeaf 相同的算法，因为双叶模型主要用于光合作用的区分，
辐射传输仍然使用总LAI计算。
"""
function radiative_transfer!(leaf::TwoBigLeaf{T}, Rn::T; kA::T=T(0.7)) where {T}
  (; Lai) = leaf
  τ = Lai <= 0.01 ? 1.0 : exp(-kA * Lai) # 穿过的比例
  leaf.Rn_c = (1 - τ) * Rn
  leaf.Rn_s = τ * Rn
end


@with_kw mutable struct Leaves{FT} <: AbstractLeaf{FT}
  nlyr::Int = 10

  Lai_sunlit::Vector{FT} = zeros(FT, nlyr)
  Lai_shaded::Vector{FT} = zeros(FT, nlyr)

  Tleaf_sunlit::Vector{FT} = zeros(FT, nlyr)
  Tleaf_shaded::Vector{FT} = zeros(FT, nlyr)

  gs_sunlit::Vector{FT} = zeros(FT, nlyr) # stomatal conductance for h2o
  gs_shaded::Vector{FT} = zeros(FT, nlyr) # stomatal conductance for h2o

  Rs::Vector{FT} = zeros(FT, nlyr + 1) # the last is ground
  Rl::Vector{FT} = zeros(FT, nlyr + 1)
  Rn::Vector{FT} = zeros(FT, nlyr + 1)

  LE_sunlit::Vector{FT} = zeros(FT, nlyr)
  LE_shaded::Vector{FT} = zeros(FT, nlyr)

  GPP_sunlit::Vector{FT} = zeros(FT, nlyr)
  GPP_shaded::Vector{FT} = zeros(FT, nlyr)
end


# Ec_CanopyTrans
# Es_EvapSoil

"""
    OverUnderCanopy{FT}

A two-layer (overstory and understory) canopy model structure.
"""
@with_kw mutable struct OverUnderCanopy{FT} <: AbstractLeaf{FT}
  # Inputs
  LAI_over::FT = 1.5   # Overstory LAI [m2 m-2]
  LAI_under::FT = 0.5  # Understory LAI [m2 m-2]

  # Partitioned Radiation
  Rn_over::FT = 0.0    # Net radiation for overstory [W m-2]
  Rn_under::FT = 0.0   # Net radiation for understory [W m-2]
  Rn_soil::FT = 0.0    # Net radiation for soil [W m-2]

  # Outputs
  GPP_over::FT = 0.0   # GPP from overstory [gC m-2 d-1]
  GPP_under::FT = 0.0  # GPP from understory [gC m-2 d-1]
  Ec_over::FT = 0.0    # Transpiration from overstory [mm d-1]
  Ec_under::FT = 0.0   # Transpiration from understory [mm d-1]
  Es::FT = 0.0         # Soil evaporation [mm d-1]
end

"""
    radiative_transfer_2layer!(canopy::OverUnderCanopy{T}, Rn::T; kA::T=T(0.7)) where {T}

Partition net radiation `Rn` into overstory, understory, and soil based on Beer's Law.
"""
function radiative_transfer_2layer!(canopy::OverUnderCanopy{T}, Rn::T; kA::T=T(0.7)) where {T}
  (; LAI_over, LAI_under) = canopy
  
  # Radiation partitioning based on Beer's Law
  τ_over = exp(-kA * LAI_over)
  τ_under = exp(-kA * LAI_under)

  # Energy absorbed by overstory
  canopy.Rn_over = Rn * (1 - τ_over)
  
  # Energy that penetrates overstory and reaches understory
  Rn_penetrated = Rn * τ_over
  
  # Energy absorbed by understory
  canopy.Rn_under = Rn_penetrated * (1 - τ_under)

  # Energy that penetrates both layers and reaches the soil
  canopy.Rn_soil = Rn_penetrated * τ_under
end


# Include unified canopy interface
include("../Canopy/canopy_interface.jl")
# Include VCmax vertical distribution functions
include("../Canopy/VCmax.jl")
# Include TwoLeaf LAI allocation functions
include("../Canopy/twoleaf_allocation.jl")
# Include TwoBigLeaf LAI allocation functions
include("../Canopy/twobigleaf_allocation.jl")
# Note: Multilayer functions are included in modules/multilayer.jl after photosynthesis.jl
