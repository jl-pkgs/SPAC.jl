using SPAC

println("\n=== TwoBigLeaf LAI Allocation Test ===")
FT = Float64

# 测试正常条件
canopy = TwoBigLeaf{FT}(Lai=5.0, Ω=1.0)
CosZs = 0.866  # 太阳天顶角 30°

allocate_LAI!(canopy, CosZs)

println("Total LAI: ", canopy.Lai)
println("Sunlit LAI: ", canopy.Lai_sunlit)
println("Shaded LAI: ", canopy.Lai_shaded)
println("LAI 守恒检查: ", canopy.Lai_sunlit + canopy.Lai_shaded ≈ canopy.Lai)

@assert canopy.Lai_sunlit > 0 "向阳叶 LAI 应该大于零"
@assert canopy.Lai_shaded > 0 "背阴叶 LAI 应该大于零"
@assert canopy.Lai_sunlit + canopy.Lai_shaded ≈ canopy.Lai rtol=1e-6 "LAI 总和应该等于输入值"

# 测试夜间条件
println("\n=== TwoBigLeaf Nighttime Test ===")
canopy_night = TwoBigLeaf{FT}(Lai=5.0, Ω=1.0)
allocate_LAI!(canopy_night, -0.1)  # 负值表示夜间

println("Nighttime - Sunlit LAI: ", canopy_night.Lai_sunlit)
println("Nighttime - Shaded LAI: ", canopy_night.Lai_shaded)

@assert canopy_night.Lai_sunlit ≈ 0.0 "夜间向阳叶 LAI 应该为零"
@assert canopy_night.Lai_shaded ≈ 5.0 "夜间背阴叶 LAI 应该等于总 LAI"

# 测试零 LAI
println("\n=== TwoBigLeaf Zero LAI Test ===")
canopy_zero = TwoBigLeaf{FT}(Lai=0.0, Ω=1.0)
allocate_LAI!(canopy_zero, CosZs)

println("Zero LAI - Sunlit: ", canopy_zero.Lai_sunlit)
println("Zero LAI - Shaded: ", canopy_zero.Lai_shaded)

@assert canopy_zero.Lai_sunlit ≈ 0.0 "零 LAI 时向阳叶应该为零"
@assert canopy_zero.Lai_shaded ≈ 0.0 "零 LAI 时背阴叶应该为零"

# 测试不同 Ω 值
println("\n=== TwoBigLeaf Clumping Index Test ===")
canopy_clustered = TwoBigLeaf{FT}(Lai=5.0, Ω=0.7)  # 聚集冠层
allocate_LAI!(canopy_clustered, CosZs)

println("聚集冠层 (Ω=0.7) - Sunlit LAI: ", canopy_clustered.Lai_sunlit)
println("聚集冠层 (Ω=0.7) - Shaded LAI: ", canopy_clustered.Lai_shaded)

canopy_random = TwoBigLeaf{FT}(Lai=5.0, Ω=1.0)  # 随机分布
allocate_LAI!(canopy_random, CosZs)

println("随机分布 (Ω=1.0) - Sunlit LAI: ", canopy_random.Lai_sunlit)
println("随机分布 (Ω=1.0) - Shaded LAI: ", canopy_random.Lai_shaded)

# 聚集冠层应该有更多的向阳叶（消光系数更小）
@assert canopy_clustered.Lai_sunlit > canopy_random.Lai_sunlit "聚集冠层向阳叶应该更多"

# 测试叶片导度
println("\n=== TwoBigLeaf Leaf Conductance Test ===")
stomatal = Stomatal_Yu2004{FT}()
photo = Photosynthesis_Rong2018{FT}()
canopy2 = TwoBigLeaf{FT}(Lai=3.5, Ω=0.9)

air = AirLayer{FT}(
  Prcp=0.0, Tavg=25.0, Rs=200.0, Rn=150.0,
  VPD=1.5, U2=2.0, Pa=101.3, Ca=410.0)

allocate_LAI!(canopy2, CosZs)

GPP, rs = leaf_conductance(air, canopy2, photo, stomatal)

println("总 GPP: ", GPP, " gC m⁻² d⁻¹")
println("向阳叶 GPP: ", canopy2.GPP_sunlit, " gC m⁻² d⁻¹")
println("背阴叶 GPP: ", canopy2.GPP_shaded, " gC m⁻² d⁻¹")
println("向阳叶导度: ", canopy2.gs_sunlit, " mol m⁻² s⁻¹")
println("背阴叶导度: ", canopy2.gs_shaded, " mol m⁻² s⁻¹")
println("冠层阻力 rs: ", rs, " s m⁻¹")

@assert GPP > 0 "GPP 应该大于零"
@assert rs > 0 "阻力应该大于零"
@assert canopy2.GPP_sunlit > 0 "向阳叶应该有 GPP"
@assert canopy2.GPP_shaded >= 0 "背阴叶可能有少量 GPP"
@assert canopy2.GPP_sunlit + canopy2.GPP_shaded ≈ GPP rtol=1e-6 "GPP 守恒"
@assert canopy2.gs_sunlit > 0 "向阳叶导度应该大于零"
@assert canopy2.gs_shaded >= 0 "背阴叶导度应该非负"

# 测试夜间导度
println("\n=== TwoBigLeaf Nighttime Conductance Test ===")
canopy3 = TwoBigLeaf{FT}(Lai=3.5, Ω=0.9)
air_night = AirLayer{FT}(
  Prcp=0.0, Tavg=20.0, Rs=0.0, Rn=50.0,
  VPD=0.5, U2=1.0, Pa=101.3, Ca=410.0)

allocate_LAI!(canopy3, -0.1)  # 夜间
GPP_night, rs_night = leaf_conductance(air_night, canopy3, photo, stomatal)

println("夜间 GPP: ", GPP_night, " gC m⁻² d⁻¹")
println("夜间向阳叶 GPP: ", canopy3.GPP_sunlit, " gC m⁻² d⁻¹")
println("夜间背阴叶 GPP: ", canopy3.GPP_shaded, " gC m⁻² d⁻¹")

@assert GPP_night ≈ 0.0 atol=1e-6 "夜间 GPP 应该接近零"
@assert canopy3.GPP_sunlit ≈ 0.0 "夜间向阳叶 GPP 应该为零"

# 对比 TwoLeaf 和 TwoBigLeaf 的结果
println("\n=== TwoLeaf vs TwoBigLeaf Comparison ===")
twoleaf = TwoLeaf{FT}(Lai=3.5, Ω=0.9)
twobigleaf = TwoBigLeaf{FT}(Lai=3.5, Ω=0.9)

allocate_LAI!(twoleaf, CosZs)
allocate_LAI!(twobigleaf, CosZs)

GPP_twoleaf, rs_twoleaf = leaf_conductance(air, twoleaf, photo, stomatal)
GPP_twobigleaf, rs_twobigleaf = leaf_conductance(air, twobigleaf, photo, stomatal)

println("TwoLeaf GPP: ", GPP_twoleaf, " gC m⁻² d⁻¹")
println("TwoBigLeaf GPP: ", GPP_twobigleaf, " gC m⁻² d⁻¹")
println("TwoLeaf rs: ", rs_twoleaf, " s m⁻¹")
println("TwoBigLeaf rs: ", rs_twobigleaf, " s m⁻¹")

println("GPP 差异: ", abs(GPP_twoleaf - GPP_twobigleaf) / GPP_twoleaf * 100, "%")
println("rs 差异: ", abs(rs_twoleaf - rs_twobigleaf) / rs_twoleaf * 100, "%")

# GPP 应该相同（使用相同的光合作用模型）
@assert GPP_twoleaf ≈ GPP_twobigleaf rtol=1e-3 "TwoLeaf 和 TwoBigLeaf 的 GPP 应该接近"
# rs 可能不同（导度计算策略不同）
println("注: rs 差异是由于导度积分策略不同（TwoLeaf 加权平均 vs TwoBigLeaf 积分）")

println("\n✅ 所有 TwoBigLeaf 测试通过！")
