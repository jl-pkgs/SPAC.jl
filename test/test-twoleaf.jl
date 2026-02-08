using SPAC

println("\n=== TwoLeaf LAI Allocation Test ===")
FT = Float64

# 测试正常条件
leaf = TwoLeaf{FT}(Lai=5.0, Ω=1.0)
CosZs = 0.866  # 太阳天顶角 30°

allocate_LAI!(leaf, CosZs)

println("Total LAI: ", leaf.Lai)
println("Sunlit LAI: ", leaf.Lai_sunlit)
println("Shaded LAI: ", leaf.Lai_shaded)
println("LAI 守恒检查: ", leaf.Lai_sunlit + leaf.Lai_shaded ≈ 5.0)

@assert leaf.Lai_sunlit > 0 "向阳叶 LAI 应该大于零"
@assert leaf.Lai_shaded > 0 "背阴叶 LAI 应该大于零"
@assert leaf.Lai_sunlit + leaf.Lai_shaded ≈ 5.0 rtol=1e-6 "LAI 总和应该等于输入值"

# 测试夜间条件
println("\n=== TwoLeaf Nighttime Test ===")
leaf_night = TwoLeaf{FT}(Lai=5.0, Ω=1.0)
allocate_LAI!(leaf_night, -0.1)  # 负值表示夜间

println("Nighttime - Sunlit LAI: ", leaf_night.Lai_sunlit)
println("Nighttime - Shaded LAI: ", leaf_night.Lai_shaded)

@assert leaf_night.Lai_sunlit ≈ 0.0 "夜间向阳叶 LAI 应该为零"
@assert leaf_night.Lai_shaded ≈ 5.0 "夜间背阴叶 LAI 应该等于总 LAI"

# 测试零 LAI
println("\n=== TwoLeaf Zero LAI Test ===")
leaf_zero = TwoLeaf{FT}(Lai=0.0, Ω=1.0)
allocate_LAI!(leaf_zero, CosZs)

println("Zero LAI - Sunlit: ", leaf_zero.Lai_sunlit)
println("Zero LAI - Shaded: ", leaf_zero.Lai_shaded)

@assert leaf_zero.Lai_sunlit ≈ 0.0 "零 LAI 时向阳叶应该为零"
@assert leaf_zero.Lai_shaded ≈ 0.0 "零 LAI 时背阴叶应该为零"

# 测试光合作用
println("\n=== TwoLeaf Photosynthesis Test ===")
photo = Photosynthesis_Rong2018{FT}()
canopy = TwoLeaf{FT}(Lai=5.0, Ω=1.0)

allocate_LAI!(canopy, CosZs)

Ag, Rd = photosynthesis(photo, 25.0, 200.0, 2.0, canopy, 380.0, 1.0)

println("总光合速率 Ag: ", Ag, " μmol m⁻² s⁻¹")
println("呼吸速率 Rd: ", Rd, " μmol m⁻² s⁻¹")

@assert Ag > 0 "白天应该有光合作用"
@assert Rd > 0 "呼吸应该大于零"

# 测试夜间光合作用
canopy_night = TwoLeaf{FT}(Lai=5.0, Ω=1.0)
allocate_LAI!(canopy_night, -0.1)

Ag_night, Rd_night = photosynthesis(photo, 25.0, 0.0, 2.0, canopy_night, 380.0, 1.0)

println("夜间光合速率: ", Ag_night)
println("夜间呼吸速率: ", Rd_night)

@assert Ag_night ≈ 0.0 atol=1e-6 "夜间光合应该接近零"

# 测试叶片导度
println("\n=== TwoLeaf Leaf Conductance Test ===")
stomatal = Stomatal_Yu2004{FT}()
canopy2 = TwoLeaf{FT}(Lai=3.5, Ω=0.9)

air = AirLayer{FT}(
  Prcp=0.0, Tavg=25.0, Rs=200.0, Rn=150.0,
  VPD=1.5, U2=2.0, Pa=101.3, Ca=410.0)

allocate_LAI!(canopy2, CosZs)

GPP, rs = leaf_conductance(air, canopy2, photo, stomatal)

println("总 GPP: ", GPP, " gC m⁻² d⁻¹")
println("向阳叶 GPP: ", canopy2.GPP_sunlit, " gC m⁻² d⁻¹")
println("背阴叶 GPP: ", canopy2.GPP_shaded, " gC m⁻² d⁻¹")
println("冠层阻力 rs: ", rs, " s m⁻¹")

@assert GPP > 0 "GPP 应该大于零"
@assert rs > 0 "阻力应该大于零"
@assert canopy2.GPP_sunlit > 0 "向阳叶应该有 GPP"
@assert canopy2.GPP_shaded >= 0 "背阴叶可能有少量 GPP"
@assert canopy2.GPP_sunlit + canopy2.GPP_shaded ≈ GPP rtol=1e-6 "GPP 守恒"

println("\n✅ 所有 TwoLeaf 测试通过！")
