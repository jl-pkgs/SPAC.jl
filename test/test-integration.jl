using SPAC, Test

println("\n=== Integration Test: Different Canopy Types ===")
FT = Float64

# 公共参数
evap = Evapotranspiration_PML{FT}()
photo = Photosynthesis_Rong2018{FT}()
stomatal = Stomatal_Yu2004{FT}()

# 创建大气强迫数据
air = AirLayer{FT}(
  Prcp=0.0, Tavg=25.0, Rs=200.0, Rn=150.0,
  VPD=1.5, U2=2.0, Pa=101.3, Ca=410.0
)

# 创建输出对象
output = SpacOutput{FT}()

println("\n--- Test 1: BigLeaf (original) ---")
canopy_bigleaf = BigLeaf{FT}(Lai=3.5)
evapotranspiration!(output, evap, photo, stomatal, air, canopy_bigleaf)
println("BigLeaf GPP: ", output.GPP)
println("BigLeaf Ec: ", output.Ec)
println("BigLeaf rs: ", output.rs)
@test output.GPP > 0
@test output.Ec > 0
@test output.rs > 0

println("\n--- Test 2: TwoLeaf ---")
canopy_twoleaf = TwoLeaf{FT}(Lai=3.5, Ω=1.0)
evapotranspiration!(output, evap, photo, stomatal, air, canopy_twoleaf)
println("TwoLeaf GPP: ", output.GPP)
println("TwoLeaf GPP_sunlit: ", canopy_twoleaf.GPP_sunlit)
println("TwoLeaf GPP_shaded: ", canopy_twoleaf.GPP_shaded)
println("TwoLeaf Ec: ", output.Ec)
println("TwoLeaf rs: ", output.rs)
@test output.GPP > 0
@test output.Ec > 0
@test output.rs > 0
@test canopy_twoleaf.GPP_sunlit > 0
@test canopy_twoleaf.GPP_shaded >= 0
@test canopy_twoleaf.GPP_sunlit + canopy_twoleaf.GPP_shaded ≈ output.GPP rtol=1e-6

println("\n--- Test 3: TwoBigLeaf ---")
canopy_twobigleaf = TwoBigLeaf{FT}(Lai=3.5, Ω=1.0)
evapotranspiration!(output, evap, photo, stomatal, air, canopy_twobigleaf)
println("TwoBigLeaf GPP: ", output.GPP)
println("TwoBigLeaf GPP_sunlit: ", canopy_twobigleaf.GPP_sunlit)
println("TwoBigLeaf GPP_shaded: ", canopy_twobigleaf.GPP_shaded)
println("TwoBigLeaf Ec: ", output.Ec)
println("TwoBigLeaf rs: ", output.rs)
@test output.GPP > 0
@test output.Ec > 0
@test output.rs > 0
@test canopy_twobigleaf.GPP_sunlit > 0
@test canopy_twobigleaf.GPP_shaded >= 0

println("\n--- Test 4: OverUnderCanopy (SKIPPED - needs leaf_conductance implementation) ---")
# TODO: OverUnderCanopy needs a special leaf_conductance implementation
# canopy_overunder = OverUnderCanopy{FT}(LAI_over=2.0, LAI_under=0.5)
# evapotranspiration!(output, evap, photo, stomatal, air, canopy_overunder)
# println("OverUnderCanopy GPP: ", output.GPP)
# println("OverUnderCanopy Ec: ", output.Ec)
# println("OverUnderCanopy rs: ", output.rs)
# @test output.GPP > 0
# @test output.Ec > 0
# @test output.rs > 0
println("Skipping OverUnderCanopy test for now")

println("\n--- Test 5: Nighttime conditions ---")
air_night = AirLayer{FT}(
  Prcp=0.0, Tavg=20.0, Rs=0.0, Rn=50.0,
  VPD=0.5, U2=1.0, Pa=101.3, Ca=410.0
)
canopy_night = TwoLeaf{FT}(Lai=3.5, Ω=1.0)
evapotranspiration!(output, evap, photo, stomatal, air_night, canopy_night)
println("Nighttime GPP: ", output.GPP)
println("Nighttime Lai_sunlit: ", canopy_night.Lai_sunlit)
println("Nighttime Lai_shaded: ", canopy_night.Lai_shaded)
@test canopy_night.Lai_sunlit ≈ 0.0  # 夜间无向阳叶
@test canopy_night.Lai_shaded ≈ 3.5  # 所有叶片为背阴叶

println("\n--- Test 6: Zero LAI ---")
air_zero = AirLayer{FT}(
  Prcp=0.0, Tavg=25.0, Rs=200.0, Rn=150.0,
  VPD=1.5, U2=2.0, Pa=101.3, Ca=410.0
)
canopy_zero = TwoLeaf{FT}(Lai=0.0, Ω=1.0)
evapotranspiration!(output, evap, photo, stomatal, air_zero, canopy_zero)
println("Zero LAI GPP: ", output.GPP)
println("Zero LAI rs: ", output.rs)
@test output.GPP ≈ 0.0  # 无冠层无光合
@test output.rs ≈ 0.0  # 无冠层无阻力

println("\n--- Test 7: calculate_CosZs function ---")
# 测试白天
CosZs_day = calculate_CosZs(200.0, 25.0, 101.3)
println("Daytime CosZs: ", CosZs_day)
@test 0.0 < CosZs_day <= 1.0

# 测试夜间
CosZs_night = calculate_CosZs(0.0, 25.0, 101.3)
println("Nighttime CosZs: ", CosZs_night)
@test CosZs_night ≈ 0.0

# 测试低辐射
CosZs_low = calculate_CosZs(10.0, 25.0, 101.3)
println("Low radiation CosZs: ", CosZs_low)
@test 0.0 <= CosZs_low <= 1.0

println("\n--- Test 8: Backward compatibility ---")
# 确保旧代码仍能工作
canopy_old = BigLeaf{FT}(Lai=2.5)
evapotranspiration!(output, evap, photo, stomatal, air, canopy_old)
println("Old code GPP: ", output.GPP)
println("Old code Ec: ", output.Ec)
@test output.GPP > 0
@test output.Ec > 0

println("\n✅ All integration tests passed!")
println("\nSummary:")
println("- BigLeaf: ✓")
println("- TwoLeaf: ✓")
println("- TwoBigLeaf: ✓")
println("- OverUnderCanopy: ⚠ (skipped - needs leaf_conductance implementation)")
println("- Nighttime conditions: ✓")
println("- Zero LAI: ✓")
println("- calculate_CosZs: ✓")
println("- Backward compatibility: ✓")
