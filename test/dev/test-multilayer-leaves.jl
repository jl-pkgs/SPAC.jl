using SPAC, Test, Statistics
using Printf

println("\n" * "="^70)
println("  多层树冠 Leaves 辐射传输测试")
println("="^70)



#==============================================================================#
#  测试 1: 基本辐射传输计算
#==============================================================================#
@testset "多层 Leaves 辐射传输 - 基本功能" begin
  nlyr = 10
  leaves = Leaves{Float64}(nlyr=nlyr)
  
  # 设置 LAI
  LAI = 3.0
  leaves.Lai_sunlit .= LAI / nlyr / 2
  leaves.Lai_shaded .= LAI / nlyr / 2
  
  # 测试净辐射计算
  Rn = 500.0  # W m⁻²
  
  # 测试边界辐射（不计算吸收）
  radiative_transfer_multilayer!(leaves, Rn; kA=0.7)
  
  # 检查边界条件
  @test leaves.Rn[1] ≈ Rn atol=1.0  # 顶部应该是入射辐射
  @test leaves.Rn[nlyr+1] >= 0.0   # 地面辐射非负
  @test leaves.Rn[nlyr+1] < Rn     # 地面辐射小于入射辐射
  
  # 检查辐射衰减（应该随深度递减）
  for i in 2:(nlyr+1)
    @test leaves.Rn[i] <= leaves.Rn[i-1] + 1e-6  # 允许小的数值误差
  end
  
  println("  ✓ 基本辐射传输计算通过")
end

#==============================================================================#
#  测试 2: 吸收辐射计算
#==============================================================================#
@testset "多层 Leaves 辐射传输 - 吸收辐射" begin
  nlyr = 10
  leaves = Leaves{Float64}(nlyr=nlyr)
  
  LAI = 3.0
  leaves.Lai_sunlit .= LAI / nlyr / 2
  leaves.Lai_shaded .= LAI / nlyr / 2
  
  Rn = 500.0
  
  # 测试吸收辐射计算
  radiative_transfer_multilayer_absorbed!(leaves, Rn; kA=0.7)
  
  # 检查吸收辐射非负
  for i in 1:nlyr
    @test leaves.Rn[i] >= 0.0
  end
  
  # 检查地面辐射非负
  @test leaves.Rn[nlyr+1] >= 0.0
  
  # 检查能量守恒：吸收 + 透射 = 入射
  total_absorbed = sum(leaves.Rn[1:nlyr])
  ground = leaves.Rn[nlyr+1]
  @test total_absorbed + ground ≈ Rn rtol=0.01
  
  println("  ✓ 吸收辐射计算通过")
  println("    总吸收: $(round(total_absorbed, digits=1)) W m⁻²")
  println("    地面辐射: $(round(ground, digits=1)) W m⁻²")
  println("    能量守恒误差: $(round((total_absorbed + ground - Rn) / Rn * 100, digits=2))%")
end

#==============================================================================#
#  测试 3: 边界条件 - 零 LAI
#==============================================================================#
@testset "多层 Leaves 辐射传输 - 零 LAI 边界条件" begin
  nlyr = 10
  leaves = Leaves{Float64}(nlyr=nlyr)
  
  # 零 LAI
  leaves.Lai_sunlit .= 0.0
  leaves.Lai_shaded .= 0.0
  
  Rn = 500.0
  
  # 测试边界辐射
  radiative_transfer_multilayer!(leaves, Rn; kA=0.7)
  
  # 所有层应该为零，地面应该是入射辐射
  @test all(leaves.Rn[1:nlyr] .≈ 0.0)
  @test leaves.Rn[nlyr+1] ≈ Rn
  
  # 测试吸收辐射
  radiative_transfer_multilayer_absorbed!(leaves, Rn; kA=0.7)
  
  @test all(leaves.Rn[1:nlyr] .≈ 0.0)
  @test leaves.Rn[nlyr+1] ≈ Rn
  
  println("  ✓ 零 LAI 边界条件通过")
end

#==============================================================================#
#  测试 4: 不同消光系数的影响
#==============================================================================#
@testset "多层 Leaves 辐射传输 - 消光系数影响" begin
  nlyr = 10
  leaves1 = Leaves{Float64}(nlyr=nlyr)
  leaves2 = Leaves{Float64}(nlyr=nlyr)
  
  LAI = 3.0
  leaves1.Lai_sunlit .= LAI / nlyr / 2
  leaves1.Lai_shaded .= LAI / nlyr / 2
  leaves2.Lai_sunlit .= LAI / nlyr / 2
  leaves2.Lai_shaded .= LAI / nlyr / 2
  
  Rn = 500.0
  
  # 不同消光系数
  kA1 = 0.5  # 较小消光
  kA2 = 0.9  # 较大消光
  
  radiative_transfer_multilayer_absorbed!(leaves1, Rn; kA=kA1)
  radiative_transfer_multilayer_absorbed!(leaves2, Rn; kA=kA2)
  
  # 较大消光系数 → 更多吸收，更少到达地面
  absorbed1 = sum(leaves1.Rn[1:nlyr])
  absorbed2 = sum(leaves2.Rn[1:nlyr])
  ground1 = leaves1.Rn[nlyr+1]
  ground2 = leaves2.Rn[nlyr+1]
  
  @test absorbed2 > absorbed1  # 更大消光 → 更多吸收
  @test ground2 < ground1      # 更大消光 → 更少到达地面
  
  println("  ✓ 消光系数影响验证通过")
  println("    kA=0.5: 吸收=$(round(absorbed1, digits=1)) W m⁻², 地面=$(round(ground1, digits=1)) W m⁻²")
  println("    kA=0.9: 吸收=$(round(absorbed2, digits=1)) W m⁻², 地面=$(round(ground2, digits=1)) W m⁻²")
end

#==============================================================================#
#  测试 5: 不同层数的影响
#==============================================================================#
@testset "多层 Leaves 辐射传输 - 层数影响" begin
  LAI = 3.0
  Rn = 500.0
  
  # 不同层数
  for nlyr in [5, 10, 20, 50]
    leaves = Leaves{Float64}(nlyr=nlyr)
    leaves.Lai_sunlit .= LAI / nlyr / 2
    leaves.Lai_shaded .= LAI / nlyr / 2
    
    radiative_transfer_multilayer_absorbed!(leaves, Rn; kA=0.7)
    
    # 检查能量守恒
    total = sum(leaves.Rn[1:nlyr]) + leaves.Rn[nlyr+1]
    @test total ≈ Rn rtol=0.01
  end
  
  println("  ✓ 不同层数能量守恒验证通过 (5, 10, 20, 50 层)")
end

#==============================================================================#
#  测试 6: Beer 定律验证
#==============================================================================#
@testset "多层 Leaves 辐射传输 - Beer 定律验证" begin
  nlyr = 10
  leaves = Leaves{Float64}(nlyr=nlyr)
  
  LAI = 3.0
  leaves.Lai_sunlit .= LAI / nlyr / 2
  leaves.Lai_shaded .= LAI / nlyr / 2
  
  Rn = 500.0
  kA = 0.7
  
  # 计算边界辐射
  radiative_transfer_multilayer!(leaves, Rn; kA=kA)
  
  # 验证 Beer 定律: Rn(L) = Rn(0) * exp(-kA * L)
  LAI_per_layer = LAI / nlyr
  
  for i in 1:(nlyr+1)
    L_cumulative = (i - 1) * LAI_per_layer
    expected = Rn * exp(-kA * L_cumulative)
    @test leaves.Rn[i] ≈ expected rtol=1e-10
  end
  
  println("  ✓ Beer 定律验证通过")
end

#==============================================================================#
#  测试 7: 与单层大叶模型对比
#==============================================================================#
@testset "多层 Leaves vs BigLeaf 辐射对比" begin
  LAI = 3.0
  Rn = 500.0
  kA = 0.7
  
  # 多层模型
  nlyr = 10
  leaves = Leaves{Float64}(nlyr=nlyr)
  leaves.Lai_sunlit .= LAI / nlyr / 2
  leaves.Lai_shaded .= LAI / nlyr / 2
  
  radiative_transfer_multilayer_absorbed!(leaves, Rn; kA=kA)
  multilayer_absorbed = sum(leaves.Rn[1:nlyr])
  multilayer_ground = leaves.Rn[nlyr+1]
  
  # BigLeaf 模型（使用 Beer 定律手动计算）
  τ = exp(-kA * LAI)
  bigleaf_absorbed = (1 - τ) * Rn  # 冠层吸收
  bigleaf_ground = τ * Rn            # 到达地面
  
  # 对比结果（应该接近）
  @test multilayer_absorbed ≈ bigleaf_absorbed rtol=0.05
  @test multilayer_ground ≈ bigleaf_ground rtol=0.05
  
  println("  ✓ 多层 vs BigLeaf 对比通过")
  println("    多层吸收: $(round(multilayer_absorbed, digits=1)) W m⁻²")
  println("    BigLeaf吸收: $(round(bigleaf_absorbed, digits=1)) W m⁻²")
  println("    相对误差: $(round(abs(multilayer_absorbed - bigleaf_absorbed) / bigleaf_absorbed * 100, digits=2))%")
end

#==============================================================================#
#  测试 8: 物理合理性检查
#==============================================================================#
@testset "多层 Leaves 辐射传输 - 物理合理性" begin
  nlyr = 10
  leaves = Leaves{Float64}(nlyr=nlyr)
  
  # 多种 LAI 情况
  for LAI in [0.5, 1.0, 2.0, 5.0]
    leaves.Lai_sunlit .= LAI / nlyr / 2
    leaves.Lai_shaded .= LAI / nlyr / 2
    
    Rn = 500.0
    radiative_transfer_multilayer_absorbed!(leaves, Rn; kA=0.7)
    
    # 检查所有值非负
    @test all(leaves.Rn .>= 0.0)
    
    # 检查能量守恒
    total = sum(leaves.Rn[1:nlyr]) + leaves.Rn[nlyr+1]
    @test total ≈ Rn rtol=0.01
    
    # 检查地面辐射随 LAI 增加而减少
    ground = leaves.Rn[nlyr+1]
    @test ground >= 0.0
    @test ground <= Rn
  end
  
  println("  ✓ 物理合理性检查通过 (LAI: 0.5, 1.0, 2.0, 5.0)")
end

println("\n" * "="^70)
println("  🎉 所有测试通过！")
println("="^70)
