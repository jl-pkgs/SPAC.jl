@testset "VCmax vertical distribution" begin
  using SPAC

  FT = Float64

  # 测试参数（基于Chen et al., 2012）
  nlyr = 10
  lai = 5.0
  Ω = 1.0
  CosZs = 0.866  # cos(30°) - 太阳天顶角30度
  VCmax25 = 62.5  # μmol m⁻² s⁻¹
  N_leaf = 3.10   # g m⁻²
  slope = 20.72 / 62.5  # 归一化斜率

  @testset "VCmax_profile basic functionality" begin
    VCmax_sunlit, VCmax_shaded, LAI_sunlit, LAI_shaded =
      VCmax_profile(nlyr, lai, Ω, CosZs, VCmax25, N_leaf, slope)

    # 检查输出数组维度
    @test length(VCmax_sunlit) == nlyr
    @test length(VCmax_shaded) == nlyr
    @test length(LAI_sunlit) == nlyr
    @test length(LAI_shaded) == nlyr

    # 检查LAI总和
    @test sum(LAI_sunlit) + sum(LAI_shaded) ≈ lai rtol=1e-3

    # 检查Vcmax递减趋势（从冠层顶部到底部）
    # 由于氮素随深度递减，Vcmax应该也递减
    @test VCmax_sunlit[1] >= VCmax_sunlit[end]
    @test VCmax_shaded[1] >= VCmax_shaded[end]

    # 检查所有值为正
    @test all(VCmax_sunlit .>= 0)
    @test all(VCmax_shaded .>= 0)
    @test all(LAI_sunlit .>= 0)
    @test all(LAI_shaded .>= 0)

    # 打印结果用于验证
    println("\n=== VCmax垂向分布测试结果 ===")
    println("层数: $nlyr, 总LAI: $lai")
    println("\n各层结果:")
    for i in 1:nlyr
      println("层 $i: LAI_sunlit=$(round(LAI_sunlit[i], digits=3)), " *
              "LAI_shaded=$(round(LAI_shaded[i], digits=3)), " *
              "VCmax_sunlit=$(round(VCmax_sunlit[i], digits=2)), " *
              "VCmax_shaded=$(round(VCmax_shaded[i], digits=2))")
    end
  end

  @testset "VCmax_profile_mean" begin
    VCmax_mean, VCmax_sunlit_mean, VCmax_shaded_mean =
      VCmax_profile_mean(nlyr, lai, Ω, CosZs, VCmax25, N_leaf, slope)

    # 检查平均值在合理范围内
    @test VCmax_mean > 0
    @test VCmax_sunlit_mean > 0
    @test VCmax_shaded_mean > 0

    # 平均值应该小于或等于VCmax25 * slope * N_leaf（冠层顶部最大值）
    max_vcmax = VCmax25 * slope * N_leaf * lai * 0.3 / (1 - exp(-0.3 * lai))
    @test VCmax_mean <= max_vcmax * 1.1  # 允许10%误差

    println("\n=== 冠层平均Vcmax ===")
    println("VCmax_mean: $(round(VCmax_mean, digits=2)) μmol m⁻² s⁻¹")
    println("VCmax_sunlit_mean: $(round(VCmax_sunlit_mean, digits=2)) μmol m⁻² s⁻¹")
    println("VCmax_shaded_mean: $(round(VCmax_shaded_mean, digits=2)) μmol m⁻² s⁻¹")
  end

  @testset "VCmax_profile comparison with original VCmax" begin
    # 比较新函数与原始向阳/背阴函数的结果
    VCmax_sunlit_old, VCmax_shaded_old = VCmax(lai, Ω, CosZs, VCmax25, N_leaf, slope)

    VCmax_mean, VCmax_sunlit_mean, VCmax_shaded_mean =
      VCmax_profile_mean(nlyr, lai, Ω, CosZs, VCmax25, N_leaf, slope)

    println("\n=== 新旧方法比较 ===")
    println("原始方法 - VCmax_sunlit: $(round(VCmax_sunlit_old, digits=2))")
    println("原始方法 - VCmax_shaded: $(round(VCmax_shaded_old, digits=2))")
    println("垂向分布 - VCmax_sunlit_mean: $(round(VCmax_sunlit_mean, digits=2))")
    println("垂向分布 - VCmax_shaded_mean: $(round(VCmax_shaded_mean, digits=2))")
    println("垂向分布 - VCmax_mean: $(round(VCmax_mean, digits=2))")
  end

  @testset "Edge cases" begin
    # 测试夜间情况（CosZs <= 0）
    VCmax_sunlit_night, VCmax_shaded_night, LAI_sunlit_night, LAI_shaded_night =
      VCmax_profile(nlyr, lai, Ω, -0.1, VCmax25, N_leaf, slope)

    @test all(VCmax_sunlit_night .== 0)
    @test all(VCmax_shaded_night .== 0)
    @test all(LAI_sunlit_night .== 0)
    @test all(LAI_shaded_night .== 0)

    # 测试零LAI情况
    VCmax_mean_zero, _, _ = VCmax_profile_mean(nlyr, 0.0, Ω, CosZs, VCmax25, N_leaf, slope)
    @test VCmax_mean_zero == VCmax25

    println("\n=== 边界条件测试 ===")
    println("夜间情况: 所有Vcmax = 0 ✓")
    println("零LAI情况: VCmax_mean = VCmax25 ✓")
  end

  @testset "Different layer numbers" begin
    # 测试不同层数的结果一致性
    for n in [5, 10, 20]
      VCmax_mean_n, _, _ = VCmax_profile_mean(n, lai, Ω, CosZs, VCmax25, N_leaf, slope)
      println("层数 $n: VCmax_mean = $(round(VCmax_mean_n, digits=2))")
    end

    # 层数越多，结果应该越接近真实的连续分布
    VCmax_mean_5, _, _ = VCmax_profile_mean(5, lai, Ω, CosZs, VCmax25, N_leaf, slope)
    VCmax_mean_20, _, _ = VCmax_profile_mean(20, lai, Ω, CosZs, VCmax25, N_leaf, slope)

    # 结果应该相对接近（层数增加后收敛）
    @test abs(VCmax_mean_5 - VCmax_mean_20) / VCmax_mean_20 < 0.05  # 相差小于5%
  end

  @testset "LAI conservation" begin
    # 详细检查LAI守恒
    VCmax_sunlit, VCmax_shaded, LAI_sunlit, LAI_shaded =
      VCmax_profile(nlyr, lai, Ω, CosZs, VCmax25, N_leaf, slope)

    total_LAI_computed = sum(LAI_sunlit) + sum(LAI_shaded)

    @test total_LAI_computed ≈ lai rtol=1e-6

    println("\n=== LAI守恒检查 ===")
    println("输入总LAI: $lai")
    println("计算总LAI: $(round(total_LAI_computed, digits=6))")
    println("向阳LAI: $(round(sum(LAI_sunlit), digits=3))")
    println("背阴LAI: $(round(sum(LAI_shaded), digits=3))")
  end
end
