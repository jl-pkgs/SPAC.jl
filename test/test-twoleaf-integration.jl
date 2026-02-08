using SPAC
using Printf

println("\n" * "="^70)
println("TwoLeaf 和 TwoBigLeaf 集成测试")
println("="^70)

FT = Float64

# 创建模型组件
evap = Evapotranspiration_PML{FT}()
photo = Photosynthesis_Rong2018{FT}()
stomatal = Stomatal_Yu2004{FT}()

println("\n创建测试气象数据...")
air = AirLayer{FT}(
  Prcp=1.5, Tavg=25.0, Rs=250.0, Rn=180.0,
  VPD=1.5, U2=2.0, Pa=101.3, Ca=410.0)

println("  Prcp = ", air.Prcp, " mm d⁻¹")
println("  Tavg = ", air.Tavg, " °C")
println("  Rs   = ", air.Rs, " W m⁻²")
println("  Rn   = ", air.Rn, " W m⁻²")
println("  VPD  = ", air.VPD, " kPa")

# 计算 CosZs
CosZs = calculate_CosZs(air.Rs, air.Tavg, air.Pa)
println("\n计算太阳天顶角余弦值: CosZs = ", CosZs)

###############################################################################
# 测试 1: TwoLeaf 模型
###############################################################################
println("\n" * "-"^70)
println("测试 1: TwoLeaf 模型")
println("-"^70)

canopy_two = TwoLeaf{FT}(Lai=3.5, Ω=0.9)
println("\n初始冠层参数:")
println("  总 LAI = ", canopy_two.Lai)
println("  遮挡指数 Ω = ", canopy_two.Ω)

# 分配 LAI
println("\n执行 LAI 分配...")
allocate_LAI!(canopy_two, CosZs)
println("  向阳叶 LAI = ", canopy_two.Lai_sunlit, " m² m⁻²")
println("  背阴叶 LAI = ", canopy_two.Lai_shaded, " m² m⁻²")
println("  LAI 守恒检查: ", canopy_two.Lai_sunlit + canopy_two.Lai_shaded ≈ canopy_two.Lai)

# 创建输出对象
output_two = SpacOutput{FT}()

# 运行蒸散发模拟
println("\n运行 PMLV2 蒸散发模拟...")
evapotranspiration!(output_two, evap, photo, stomatal, air, canopy_two)

println("\n模拟结果:")
println("  GPP          = ", output_two.GPP, " gC m⁻² d⁻¹")
println("  向阳叶 GPP   = ", canopy_two.GPP_sunlit, " gC m⁻² d⁻¹")
println("  背阴叶 GPP   = ", canopy_two.GPP_shaded, " gC m⁻² d⁻¹")
println("  冠层蒸腾 Ec  = ", output_two.Ec, " mm d⁻¹")
println("  截留蒸发 Ei  = ", output_two.Ei, " mm d⁻¹")
println("  土壤蒸发 Es  = ", output_two.Es_eq, " mm d⁻¹ (均衡值)")
println("  总 ET        = ", output_two.ET, " mm d⁻¹")
println("  冠层阻力 rs  = ", output_two.rs, " s m⁻¹")
println("  空气动力学阻力 ra = ", output_two.ra, " s m⁻¹")

# 验证
println("\n验证 TwoLeaf 结果:")
println("  GPP > 0: ", output_two.GPP > 0)
println("  Ec > 0: ", output_two.Ec > 0)
println("  Ei ≥ 0: ", output_two.Ei ≥ 0)
println("  向阳叶贡献 > 背阴叶: ", canopy_two.GPP_sunlit > canopy_two.GPP_shaded)
println("  GPP 守恒: ", canopy_two.GPP_sunlit + canopy_two.GPP_shaded ≈ output_two.GPP)

###############################################################################
# 测试 2: TwoBigLeaf 模型
###############################################################################
println("\n" * "-"^70)
println("测试 2: TwoBigLeaf 模型")
println("-"^70)

canopy_twobig = TwoBigLeaf{FT}(Lai=3.5, Ω=0.9)
println("\n初始冠层参数:")
println("  总 LAI = ", canopy_twobig.Lai)
println("  遮挡指数 Ω = ", canopy_twobig.Ω)

# 分配 LAI
println("\n执行 LAI 分配...")
allocate_LAI!(canopy_twobig, CosZs)
println("  向阳叶 LAI = ", canopy_twobig.Lai_sunlit, " m² m⁻²")
println("  背阴叶 LAI = ", canopy_twobig.Lai_shaded, " m² m⁻²")
println("  LAI 守恒检查: ", canopy_twobig.Lai_sunlit + canopy_twobig.Lai_shaded ≈ canopy_twobig.Lai)

# 创建输出对象
output_twobig = SpacOutput{FT}()

# 运行蒸散发模拟
println("\n运行 PMLV2 蒸散发模拟...")
evapotranspiration!(output_twobig, evap, photo, stomatal, air, canopy_twobig)

println("\n模拟结果:")
println("  GPP          = ", output_twobig.GPP, " gC m⁻² d⁻¹")
println("  向阳叶 GPP   = ", canopy_twobig.GPP_sunlit, " gC m⁻² d⁻¹")
println("  背阴叶 GPP   = ", canopy_twobig.GPP_shaded, " gC m⁻² d⁻¹")
println("  向阳叶 gs    = ", canopy_twobig.gs_sunlit, " mol m⁻² s⁻¹")
println("  背阴叶 gs    = ", canopy_twobig.gs_shaded, " mol m⁻² s⁻¹")
println("  冠层蒸腾 Ec  = ", output_twobig.Ec, " mm d⁻¹")
println("  截留蒸发 Ei  = ", output_twobig.Ei, " mm d⁻¹")
println("  土壤蒸发 Es  = ", output_twobig.Es_eq, " mm d⁻¹ (均衡值)")
println("  总 ET        = ", output_twobig.ET, " mm d⁻¹")
println("  冠层阻力 rs  = ", output_twobig.rs, " s m⁻¹")
println("  空气动力学阻力 ra = ", output_twobig.ra, " s m⁻¹")

# 验证
println("\n验证 TwoBigLeaf 结果:")
println("  GPP > 0: ", output_twobig.GPP > 0)
println("  Ec > 0: ", output_twobig.Ec > 0)
println("  Ei ≥ 0: ", output_twobig.Ei ≥ 0)
println("  向阳叶 gs > 0: ", canopy_twobig.gs_sunlit > 0)
println("  背阴叶 gs > 0: ", canopy_twobig.gs_shaded > 0)
println("  GPP 守恒: ", canopy_twobig.GPP_sunlit + canopy_twobig.GPP_shaded ≈ output_twobig.GPP)

###############################################################################
# 测试 3: 对比 BigLeaf, TwoLeaf, TwoBigLeaf
###############################################################################
println("\n" * "-"^70)
println("测试 3: 对比三种冠层模型")
println("-"^70)

canopy_big = BigLeaf{FT}(Lai=3.5)
output_big = SpacOutput{FT}()

evapotranspiration!(output_big, evap, photo, stomatal, air, canopy_big)

println("\n对比结果:")
println("\n模型         GPP          Ec           Ei          ET")
println("-" ^ 80)
println(@sprintf("BigLeaf      %6.2f       %6.2f       %6.2f       %6.2f",
  output_big.GPP, output_big.Ec, output_big.Ei, output_big.ET))
println(@sprintf("TwoLeaf      %6.2f       %6.2f       %6.2f       %6.2f",
  output_two.GPP, output_two.Ec, output_two.Ei, output_two.ET))
println(@sprintf("TwoBigLeaf   %6.2f       %6.2f       %6.2f       %6.2f",
  output_twobig.GPP, output_twobig.Ec, output_twobig.Ei, output_twobig.ET))

println("\n向阳叶/背阴叶贡献:")
println("\n模型         向阳叶GPP    背阴叶GPP    向阳叶比例   背阴叶比例")
println("-" ^ 80)
if canopy_two.GPP_sunlit + canopy_two.GPP_shaded > 0
  sunlit_frac_two = canopy_two.GPP_sunlit / output_two.GPP * 100
  shaded_frac_two = canopy_two.GPP_shaded / output_two.GPP * 100
  println(@sprintf("TwoLeaf      %6.2f       %6.2f       %6.1f%%       %6.1f%%",
    canopy_two.GPP_sunlit, canopy_two.GPP_shaded, sunlit_frac_two, shaded_frac_two))
end

if canopy_twobig.GPP_sunlit + canopy_twobig.GPP_shaded > 0
  sunlit_frac_twobig = canopy_twobig.GPP_sunlit / output_twobig.GPP * 100
  shaded_frac_twobig = canopy_twobig.GPP_shaded / output_twobig.GPP * 100
  println(@sprintf("TwoBigLeaf   %6.2f       %6.2f       %6.1f%%       %6.1f%%",
    canopy_twobig.GPP_sunlit, canopy_twobig.GPP_shaded, sunlit_frac_twobig, shaded_frac_twobig))
end

###############################################################################
# 测试 4: 夜间条件
###############################################################################
println("\n" * "-"^70)
println("测试 4: 夜间条件验证")
println("-"^70)

air_night = AirLayer{FT}(
  Prcp=0.0, Tavg=20.0, Rs=0.0, Rn=50.0,
  VPD=1.0, U2=1.0, Pa=101.3, Ca=410.0)

CosZs_night = calculate_CosZs(air_night.Rs, air_night.Tavg, air_night.Pa)
println("\n夜间气象条件:")
println("  Rs = ", air_night.Rs, " W m⁻² (夜间)")
println("  CosZs = ", CosZs_night, " (应该为0)")

canopy_two_night = TwoLeaf{FT}(Lai=3.5, Ω=0.9)
output_two_night = SpacOutput{FT}()

allocate_LAI!(canopy_two_night, CosZs_night)
println("\n夜间 LAI 分配:")
println("  向阳叶 LAI = ", canopy_two_night.Lai_sunlit, " (应该为0)")
println("  背阴叶 LAI = ", canopy_two_night.Lai_shaded, " (应该等于总LAI)")

evapotranspiration!(output_two_night, evap, photo, stomatal, air_night, canopy_two_night)

println("\n夜间模拟结果:")
println("  GPP = ", output_two_night.GPP, " gC m⁻² d⁻¹ (应该接近0)")
println("  Ec  = ", output_two_night.Ec, " mm d⁻¹ (夜间蒸腾)")
println("  Ei  = ", output_two_night.Ei, " mm d⁻¹ (截留)")

println("\n夜间验证:")
println("  GPP ≈ 0: ", isapprox(output_two_night.GPP, 0.0, atol=0.1))
println("  向阳叶 LAI = 0: ", isapprox(canopy_two_night.Lai_sunlit, 0.0))

###############################################################################
# 测试 5: 不同 LAI 条件
###############################################################################
println("\n" * "-"^70)
println("测试 5: 不同 LAI 条件敏感性")
println("-"^70)

LAI_values = [0.5, 1.0, 2.0, 3.0, 5.0, 8.0]

println("\nLAI    总GPP    向阳GPP  背阴GPP  向阳比例  背阴比例  冠层阻力")
println("-"^90)

for LAI_test in LAI_values
  canopy_test = TwoLeaf{FT}(Lai=LAI_test, Ω=0.9)
  output_test = SpacOutput{FT}()

  allocate_LAI!(canopy_test, CosZs)
  evapotranspiration!(output_test, evap, photo, stomatal, air, canopy_test)

  if output_test.GPP > 0
    sunlit_frac = canopy_test.GPP_sunlit / output_test.GPP * 100
    shaded_frac = canopy_test.GPP_shaded / output_test.GPP * 100
  else
    sunlit_frac = 0.0
    shaded_frac = 0.0
  end

  println(@sprintf("%3.1f   %6.2f    %6.2f   %6.2f    %5.1f%%     %5.1f%%   %6.1f",
    LAI_test, output_test.GPP, canopy_test.GPP_sunlit,
    canopy_test.GPP_shaded, sunlit_frac, shaded_frac, output_test.rs))
end

###############################################################################
# 最终验证
###############################################################################
println("\n" * "="^70)
println("最终验证与断言")
println("="^70)

passed = 0
total = 0

# TwoLeaf 验证
total += 1
if output_two.GPP > 0 && output_two.Ec > 0
  println("✅ TwoLeaf: GPP和Ec为正数")
  passed += 1
else
  println("❌ TwoLeaf: GPP或Ec为负数")
end

total += 1
if canopy_two.GPP_sunlit + canopy_two.GPP_shaded ≈ output_two.GPP rtol=1e-6
  println("✅ TwoLeaf: GPP守恒")
  passed += 1
else
  println("❌ TwoLeaf: GPP不守恒")
end

total += 1
if canopy_two.GPP_sunlit > canopy_two.GPP_shaded
  println("✅ TwoLeaf: 向阳叶GPP > 背阴叶GPP")
  passed += 1
else
  println("❌ TwoLeaf: 向阳叶GPP应该 > 背阴叶GPP")
end

# TwoBigLeaf 验证
total += 1
if output_twobig.GPP > 0 && output_twobig.Ec > 0
  println("✅ TwoBigLeaf: GPP和Ec为正数")
  passed += 1
else
  println("❌ TwoBigLeaf: GPP或Ec为负数")
end

total += 1
if canopy_twobig.GPP_sunlit + canopy_twobig.GPP_shaded ≈ output_twobig.GPP rtol=1e-6
  println("✅ TwoBigLeaf: GPP守恒")
  passed += 1
else
  println("❌ TwoBigLeaf: GPP不守恒")
end

total += 1
if canopy_twobig.gs_sunlit > 0 && canopy_twobig.gs_shaded > 0
  println("✅ TwoBigLeaf: 气孔导度为正数")
  passed += 1
else
  println("❌ TwoBigLeaf: 气孔导度为负数或零")
end

# 夜间验证
total += 1
if output_two_night.GPP ≈ 0.0 atol=0.1
  println("✅ 夜间: GPP接近零")
  passed += 1
else
  println("❌ 夜间: GPP应该接近零")
end

total += 1
if canopy_two_night.Lai_sunlit ≈ 0.0
  println("✅ 夜间: 向阳叶LAI为零")
  passed += 1
else
  println("❌ 夜间: 向阳叶LAI应该为零")
end

# 模型对比验证
total += 1
if abs(output_two.GPP - output_twobig.GPP) / output_big.GPP < 0.1
  println("✅ 模型对比: TwoLeaf和TwoBigLeaf的GPP差异 < 10%")
  passed += 1
else
  println("⚠️  模型对比: TwoLeaf和TwoBigLeaf的GPP差异较大 (预期差异)")
end

println("\n" * "="^70)
println("测试总结: $passed / $total 通过")
println("="^70)

if passed == total
  println("\n🎉 所有测试通过！TwoLeaf 和 TwoBigLeaf 实现正确！")
else
  println("\n⚠️  部分测试未通过，需要检查实现")
end
