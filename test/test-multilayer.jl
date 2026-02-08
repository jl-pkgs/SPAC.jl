using Test
using SPAC
using Statistics

"""
测试多层 Leaves 模型的基本功能
"""

@testset "多层 Leaves 模型测试" begin

  @testset "多层初始化" begin
    # 测试默认 10 层
    leaves = Leaves{Float64}()
    @test leaves.nlyr == 10

    # 测试自定义层数
    leaves_5 = Leaves{Float64}(nlyr=5)
    @test leaves_5.nlyr == 5

    leaves_20 = Leaves{Float64}(nlyr=20)
    @test leaves_20.nlyr == 20
  end

  @testset "多层 LAI 和 VCmax 分配" begin
    # 创建多层冠层
    nlyr = 10
    leaves = Leaves{Float64}(nlyr=nlyr)

    # 测试参数
    LAI = 5.0
    Ω = 1.0
    CosZs = 0.866  # 30° 太阳天顶角
    VCmax25 = 62.5
    N_leaf = 3.10
    slope = 20.72 / 62.5

    # 初始化多层
    VCmax_sunlit = zeros(Float64, nlyr)
    VCmax_shaded = zeros(Float64, nlyr)
    leaves, VCmax_sunlit, VCmax_shaded = initialize_multilayer_with_vcmax!(
      leaves, LAI, Ω, CosZs, VCmax25, N_leaf, slope, VCmax_sunlit, VCmax_shaded
    )

    # 检查 LAI 守恒
    total_LAI_sunlit = sum(leaves.Lai_sunlit)
    total_LAI_shaded = sum(leaves.Lai_shaded)
    total_LAI = total_LAI_sunlit + total_LAI_shaded

    @test total_LAI ≈ LAI atol=0.01

    # 检查所有 LAI 非负
    @test all(leaves.Lai_sunlit .>= 0)
    @test all(leaves.Lai_shaded .>= 0)

    # 检查 VCmax 非负
    @test all(VCmax_sunlit .>= 0)
    @test all(VCmax_shaded .>= 0)

    # 检查 VCmax 垂向分布（顶部应该更高）
    # 由于氮素从顶部到底部衰减，VCmax 也应该衰减
    if total_LAI_sunlit > 0.1
      @test VCmax_sunlit[1] > VCmax_sunlit[nlyr] * 0.5  # 顶部至少是底部的 50%
    end
  end

  @testset "边界条件：零 LAI" begin
    nlyr = 10
    leaves = Leaves{Float64}(nlyr=nlyr)

    LAI = 0.0
    Ω = 1.0
    CosZs = 0.866
    VCmax25 = 62.5
    N_leaf = 3.10
    slope = 20.72 / 62.5

    leaves, VCmax_sunlit, VCmax_shaded = initialize_multilayer!(
      leaves, LAI, Ω, CosZs, VCmax25, N_leaf, slope
    )

    @test sum(leaves.Lai_sunlit) ≈ 0.0 atol=1e-6
    @test sum(leaves.Lai_shaded) ≈ 0.0 atol=1e-6
    @test all(leaves.Lai_sunlit .≈ 0.0)
    @test all(leaves.Lai_shaded .≈ 0.0)
  end

  @testset "边界条件：夜间（CosZs ≤ 0）" begin
    nlyr = 10
    leaves = Leaves{Float64}(nlyr=nlyr)

    LAI = 5.0
    Ω = 1.0
    CosZs = 0.0  # 夜间
    VCmax25 = 62.5
    N_leaf = 3.10
    slope = 20.72 / 62.5

    leaves, VCmax_sunlit, VCmax_shaded = initialize_multilayer!(
      leaves, LAI, Ω, CosZs, VCmax25, N_leaf, slope
    )

    # 夜间所有叶片都是背阴叶
    @test sum(leaves.Lai_sunlit) ≈ 0.0 atol=1e-6
    @test sum(leaves.Lai_shaded) ≈ LAI atol=0.01
  end

  @testset "多层光合作用计算" begin
    nlyr = 10
    leaves = Leaves{Float64}(nlyr=nlyr)

    LAI = 5.0
    Ω = 1.0
    CosZs = 0.866
    VCmax25 = 62.5
    N_leaf = 3.10
    slope = 20.72 / 62.5

    # 初始化
    leaves, VCmax_sunlit, VCmax_shaded = initialize_multilayer!(
      leaves, LAI, Ω, CosZs, VCmax25, N_leaf, slope
    )

    # 创建光合作用模型和大气数据
    photo = Photosynthesis_Rong2018{Float64}()
    air = AirLayer{Float64}(
      Tavg=25.0,
      Rs=200.0,
      VPD=2.0,
      Pa=101.325,
      Ca=380.0,
      PC=1.0
    )

    # 计算多层光合作用
    GPP, Rd = photosynthesis_multilayer!(leaves, photo, air, VCmax_sunlit, VCmax_shaded)

    # 检查结果合理性
    @test GPP >= 0  # 白天 GPP 应该 >= 0
    @test Rd >= 0   # 呼吸应该 >= 0

    # 检查每层的 GPP
    @test all(leaves.GPP_sunlit .>= 0)
    @test all(leaves.GPP_shaded .>= 0)

    # 检查 GPP 总和
    GPP_sum = sum(leaves.GPP_sunlit) + sum(leaves.GPP_shaded)
    @test GPP_sum ≈ GPP atol=0.01
  end

  @testset "多层气孔导度计算" begin
    nlyr = 10
    leaves = Leaves{Float64}(nlyr=nlyr)

    LAI = 5.0
    Ω = 1.0
    CosZs = 0.866
    VCmax25 = 62.5
    N_leaf = 3.10
    slope = 20.72 / 62.5

    # 初始化
    leaves, VCmax_sunlit, VCmax_shaded = initialize_multilayer!(
      leaves, LAI, Ω, CosZs, VCmax25, N_leaf, slope
    )

    # 创建模型
    photo = Photosynthesis_Rong2018{Float64}()
    stomatal = Stomatal_BallBerry1987{Float64}()
    air = AirLayer{Float64}(
      Tavg=25.0,
      Rs=200.0,
      VPD=2.0,
      Pa=101.325,
      Ca=380.0,
      PC=1.0
    )

    # 计算多层气孔导度
    GPP, rs = stomatal_conductance_multilayer!(leaves, photo, stomatal, air, VCmax_sunlit, VCmax_shaded)

    # 检查结果合理性
    @test GPP >= 0
    @test rs >= 0  # 阻力应该非负
    @test isfinite(rs)  # 阻力应该是有限值

    # 检查每层的导度
    @test all(leaves.gs_sunlit .>= 0)
    @test all(leaves.gs_shaded .>= 0)

    # 检查 GPP 和导度的对应关系
    # 有 GPP 的层应该有导度
    for i in 1:nlyr
      if leaves.GPP_sunlit[i] > 0
        @test leaves.gs_sunlit[i] > 0
      end
      if leaves.GPP_shaded[i] > 0
        @test leaves.gs_shaded[i] > 0
      end
    end
  end

  @testset "多层辐射传输" begin
    nlyr = 10
    leaves = Leaves{Float64}(nlyr=nlyr)

    LAI = 5.0
    Ω = 1.0
    CosZs = 0.866
    VCmax25 = 62.5
    N_leaf = 3.10
    slope = 20.72 / 62.5

    # 初始化 LAI
    initialize_multilayer!(leaves, LAI, Ω, CosZs, VCmax25, N_leaf, slope)

    Rn = 150.0  # W m-2
    kA = 0.7

    # 计算辐射传输
    radiative_transfer_multilayer!(leaves, Rn; kA=kA)

    # 检查 Rn 数组长度
    @test length(leaves.Rn) == nlyr + 1

    # 检查辐射递减
    for i in 1:nlyr
      @test leaves.Rn[i] >= leaves.Rn[i+1]  # 辐射应该递减
    end

    # 检查地面辐射非负
    @test leaves.Rn[nlyr+1] >= 0

    # 检查冠层顶部辐射等于输入
    @test leaves.Rn[1] ≈ Rn atol=0.01

    # 检查吸收的辐射版本
    radiative_transfer_multilayer_absorbed!(leaves, Rn; kA=kA)

    # 检查吸收的辐射非负
    @test all(leaves.Rn[1:nlyr] .>= 0)

    # 检查吸收的辐射总和
    absorbed_sum = sum(leaves.Rn[1:nlyr])
    @test absorbed_sum <= Rn  # 吸收的辐射不应该超过总辐射
  end

  @testset "层数可配置性" begin
    # 测试不同层数
    for n in [5, 10, 20, 50]
      leaves = Leaves{Float64}(nlyr=n)

      LAI = 5.0
      Ω = 1.0
      CosZs = 0.866
      VCmax25 = 62.5
      N_leaf = 3.10
      slope = 20.72 / 62.5

      leaves, VCmax_sunlit, VCmax_shaded = initialize_multilayer!(
        leaves, LAI, Ω, CosZs, VCmax25, N_leaf, slope
      )

      # 检查数组长度
      @test length(leaves.Lai_sunlit) == n
      @test length(leaves.Lai_shaded) == n
      @test length(VCmax_sunlit) == n
      @test length(VCmax_shaded) == n

      # 检查 LAI 守恒
      total_LAI = sum(leaves.Lai_sunlit) + sum(leaves.Lai_shaded)
      @test total_LAI ≈ LAI atol=0.01
    end
  end

  @testset "VCmax 垂向分布验证" begin
    nlyr = 20
    leaves = Leaves{Float64}(nlyr=nlyr)

    LAI = 5.0
    Ω = 1.0
    CosZs = 0.866
    VCmax25 = 62.5
    N_leaf = 3.10
    slope = 20.72 / 62.5

    leaves, VCmax_sunlit, VCmax_shaded = initialize_multilayer!(
      leaves, LAI, Ω, CosZs, VCmax25, N_leaf, slope
    )

    # 检查 VCmax 从顶部到底部递减
    # 计算加权平均 VCmax
    total_LAI_sunlit = sum(leaves.Lai_sunlit)
    total_LAI_shaded = sum(leaves.Lai_shaded)

    if total_LAI_sunlit > 0.1
      VCmax_sunlit_mean = sum(VCmax_sunlit .* leaves.Lai_sunlit) / total_LAI_sunlit

      # 顶部（前 3 层）应该高于平均
      top_3_mean = mean(VCmax_sunlit[1:min(3, nlyr)])
      @test top_3_mean > VCmax_sunlit_mean * 0.9  # 顶部应该接近或高于平均

      # 底部（后 3 层）应该低于平均
      bottom_3_mean = mean(VCmax_sunlit[max(1, nlyr-2):nlyr])
      @test bottom_3_mean < VCmax_sunlit_mean * 1.1  # 底部应该接近或低于平均
    end
  end

  @testset "完整流程测试" begin
    # 测试完整的初始化 -> 光合 -> 导度流程
    nlyr = 10
    leaves = Leaves{Float64}(nlyr=nlyr)

    LAI = 5.0
    Ω = 1.0
    CosZs = 0.866
    VCmax25 = 62.5
    N_leaf = 3.10
    slope = 20.72 / 62.5

    # 1. 初始化
    leaves, VCmax_sunlit, VCmax_shaded = initialize_multilayer!(
      leaves, LAI, Ω, CosZs, VCmax25, N_leaf, slope
    )

    # 2. 辐射传输
    Rn = 150.0
    radiative_transfer_multilayer!(leaves, Rn; kA=0.7)

    # 3. 光合作用
    photo = Photosynthesis_Rong2018{Float64}()
    air = AirLayer{Float64}(
      Tavg=25.0,
      Rs=200.0,
      VPD=2.0,
      Pa=101.325,
      Ca=380.0,
      PC=1.0
    )
    GPP1, Rd1 = photosynthesis_multilayer!(leaves, photo, air, VCmax_sunlit, VCmax_shaded)

    # 4. 气孔导度
    stomatal = Stomatal_BallBerry1987{Float64}()
    GPP2, rs = stomatal_conductance_multilayer!(leaves, photo, stomatal, air, VCmax_sunlit, VCmax_shaded)

    # 检查一致性
    @test GPP1 ≈ GPP2 atol=0.1  # 两次计算的 GPP 应该接近

    # 检查数值稳定性
    @test isfinite(GPP1)
    @test isfinite(Rd1)
    @test isfinite(rs)
  end

end
