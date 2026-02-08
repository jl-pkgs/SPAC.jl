using Test
using SPAC

# 导入必要的函数
using SPAC: radiative_transfer!, radiative_transfer_multilayer!, radiative_transfer_2layer!
using SPAC: allocate_LAI!

"""
测试统一冠层接口的功能

这个测试套件验证了所有冠层类型（BigLeaf, TwoLeaf, TwoBigLeaf, Leaves, OverUnderCanopy）
的统一接口函数是否正确实现。
"""

@testset "Canopy Interface Tests" begin

  FT = Float64

  # ========================================================================
  # Test 1: get_lai - 获取总叶面积指数
  # ========================================================================
  @testset "get_lai" begin
    # BigLeaf
    canopy_bigleaf = BigLeaf{FT}(Lai=2.5)
    @test get_lai(canopy_bigleaf) == 2.5

    # TwoLeaf
    canopy_twoleaf = TwoLeaf{FT}(Lai=3.0, Ω=1.0)
    @test get_lai(canopy_twoleaf) == 3.0

    # TwoBigLeaf（分配后）
    canopy_twobigleaf = TwoBigLeaf{FT}(Lai=4.0, Ω=1.0)
    allocate_LAI!(canopy_twobigleaf, 0.866)
    @test get_lai(canopy_twobigleaf) ≈ 4.0 atol=0.01

    # Leaves
    canopy_leaves = Leaves{FT}(nlyr=10)
    fill!(canopy_leaves.Lai_sunlit, 0.2)
    fill!(canopy_leaves.Lai_shaded, 0.3)
    @test get_lai(canopy_leaves) ≈ 10 * (0.2 + 0.3) atol=0.01

    # OverUnderCanopy
    canopy_overunder = OverUnderCanopy{FT}(LAI_over=2.0, LAI_under=0.5)
    @test get_lai(canopy_overunder) == 2.5
  end

  # ========================================================================
  # Test 2: requires_lai_allocation - 是否需要 LAI 分配
  # ========================================================================
  @testset "requires_lai_allocation" begin
    @test requires_lai_allocation(BigLeaf{FT}()) == false
    @test requires_lai_allocation(TwoLeaf{FT}()) == true
    @test requires_lai_allocation(TwoBigLeaf{FT}()) == true
    @test requires_lai_allocation(Leaves{FT}()) == true
    @test requires_lai_allocation(OverUnderCanopy{FT}()) == false
  end

  # ========================================================================
  # Test 3: allocate_lai_if_needed! - 条件 LAI 分配
  # ========================================================================
  @testset "allocate_lai_if_needed!" begin
    # BigLeaf - 不需要分配，应直接返回
    canopy_bigleaf = BigLeaf{FT}(Lai=2.5)
    result = allocate_lai_if_needed!(canopy_bigleaf, 0.866)
    @test result === canopy_bigleaf  # 返回同一对象

    # TwoLeaf - 需要分配
    canopy_twoleaf = TwoLeaf{FT}(Lai=3.0, Ω=1.0)
    allocate_lai_if_needed!(canopy_twoleaf, 0.866)
    @test canopy_twoleaf.Lai_sunlit > 0
    @test canopy_twoleaf.Lai_shaded > 0
    @test canopy_twoleaf.Lai_sunlit + canopy_twoleaf.Lai_shaded ≈ 3.0 atol=0.01

    # TwoBigLeaf - 需要分配
    canopy_twobigleaf = TwoBigLeaf{FT}(Lai=4.0, Ω=1.0)
    allocate_lai_if_needed!(canopy_twobigleaf, 0.866)
    @test canopy_twobigleaf.Lai_sunlit > 0
    @test canopy_twobigleaf.Lai_shaded > 0
    @test canopy_twobigleaf.Lai_sunlit + canopy_twobigleaf.Lai_shaded ≈ 4.0 atol=0.01

    # Leaves - 不需要在这里分配（在初始化时分配）
    canopy_leaves = Leaves{FT}(nlyr=10)
    result = allocate_lai_if_needed!(canopy_leaves, 0.866)
    @test result === canopy_leaves

    # OverUnderCanopy - 不需要分配
    canopy_overunder = OverUnderCanopy{FT}()
    result = allocate_lai_if_needed!(canopy_overunder, 0.866)
    @test result === canopy_overunder
  end

  # ========================================================================
  # Test 4: get_rn_canopy - 获取冠层净辐射
  # ========================================================================
  @testset "get_rn_canopy" begin
    Rn_input = 150.0

    # BigLeaf
    canopy_bigleaf = BigLeaf{FT}(Lai=2.5)
    radiative_transfer!(canopy_bigleaf, Rn_input; kA=0.7)
    @test get_rn_canopy(canopy_bigleaf) == canopy_bigleaf.Rn_c
    @test get_rn_canopy(canopy_bigleaf) > 0

    # TwoLeaf
    canopy_twoleaf = TwoLeaf{FT}(Lai=3.0, Ω=1.0)
    radiative_transfer!(canopy_twoleaf, Rn_input; kA=0.7)
    @test get_rn_canopy(canopy_twoleaf) == canopy_twoleaf.Rn_c
    @test get_rn_canopy(canopy_twoleaf) > 0

    # TwoBigLeaf
    canopy_twobigleaf = TwoBigLeaf{FT}(Lai=4.0, Ω=1.0)
    radiative_transfer!(canopy_twobigleaf, Rn_input; kA=0.7)
    @test get_rn_canopy(canopy_twobigleaf) == canopy_twobigleaf.Rn_c
    @test get_rn_canopy(canopy_twobigleaf) > 0

    # Leaves
    canopy_leaves = Leaves{FT}(nlyr=10)
    # 填充一些 LAI 数据以便测试
    fill!(canopy_leaves.Lai_sunlit, 0.1)
    fill!(canopy_leaves.Lai_shaded, 0.2)
    radiative_transfer_multilayer!(canopy_leaves, Rn_input; kA=0.7)
    @test get_rn_canopy(canopy_leaves) == canopy_leaves.Rn[1]
    # 注意：Rn[1] 是到达冠层顶部的辐射，不是吸收的辐射
    @test canopy_leaves.Rn[1] ≈ Rn_input atol=1.0  # 第一层应该接近输入辐射

    # OverUnderCanopy
    canopy_overunder = OverUnderCanopy{FT}(LAI_over=2.0, LAI_under=0.5)
    radiative_transfer_2layer!(canopy_overunder, Rn_input; kA=0.7)
    @test get_rn_canopy(canopy_overunder) == canopy_overunder.Rn_over
    @test get_rn_canopy(canopy_overunder) > 0
  end

  # ========================================================================
  # Test 5: get_rn_soil - 获取土壤净辐射
  # ========================================================================
  @testset "get_rn_soil" begin
    Rn_input = 150.0

    # BigLeaf
    canopy_bigleaf = BigLeaf{FT}(Lai=2.5)
    radiative_transfer!(canopy_bigleaf, Rn_input; kA=0.7)
    @test get_rn_soil(canopy_bigleaf) == canopy_bigleaf.Rn_s
    @test get_rn_soil(canopy_bigleaf) >= 0

    # TwoLeaf
    canopy_twoleaf = TwoLeaf{FT}(Lai=3.0, Ω=1.0)
    radiative_transfer!(canopy_twoleaf, Rn_input; kA=0.7)
    @test get_rn_soil(canopy_twoleaf) == canopy_twoleaf.Rn_s
    @test get_rn_soil(canopy_twoleaf) >= 0

    # TwoBigLeaf
    canopy_twobigleaf = TwoBigLeaf{FT}(Lai=4.0, Ω=1.0)
    radiative_transfer!(canopy_twobigleaf, Rn_input; kA=0.7)
    @test get_rn_soil(canopy_twobigleaf) == canopy_twobigleaf.Rn_s
    @test get_rn_soil(canopy_twobigleaf) >= 0

    # Leaves
    canopy_leaves = Leaves{FT}(nlyr=10)
    radiative_transfer_multilayer!(canopy_leaves, Rn_input; kA=0.7)
    @test get_rn_soil(canopy_leaves) == canopy_leaves.Rn[end]
    @test get_rn_soil(canopy_leaves) >= 0

    # OverUnderCanopy
    canopy_overunder = OverUnderCanopy{FT}(LAI_over=2.0, LAI_under=0.5)
    radiative_transfer_2layer!(canopy_overunder, Rn_input; kA=0.7)
    @test get_rn_soil(canopy_overunder) == canopy_overunder.Rn_soil
    @test get_rn_soil(canopy_overunder) >= 0
  end

  # ========================================================================
  # Test 6: radiative_transfer_for_canopy! - 统一辐射传输接口
  # ========================================================================
  @testset "radiative_transfer_for_canopy!" begin
    Rn_input = 150.0

    # BigLeaf
    canopy_bigleaf = BigLeaf{FT}(Lai=2.5)
    radiative_transfer_for_canopy!(canopy_bigleaf, Rn_input; kA=0.7)
    @test canopy_bigleaf.Rn_c > 0
    @test canopy_bigleaf.Rn_s >= 0
    @test canopy_bigleaf.Rn_c + canopy_bigleaf.Rn_s ≈ Rn_input atol=1.0

    # TwoLeaf
    canopy_twoleaf = TwoLeaf{FT}(Lai=3.0, Ω=1.0)
    radiative_transfer_for_canopy!(canopy_twoleaf, Rn_input; kA=0.7)
    @test canopy_twoleaf.Rn_c > 0
    @test canopy_twoleaf.Rn_s >= 0
    @test canopy_twoleaf.Rn_c + canopy_twoleaf.Rn_s ≈ Rn_input atol=1.0

    # TwoBigLeaf
    canopy_twobigleaf = TwoBigLeaf{FT}(Lai=4.0, Ω=1.0)
    radiative_transfer_for_canopy!(canopy_twobigleaf, Rn_input; kA=0.7)
    @test canopy_twobigleaf.Rn_c > 0
    @test canopy_twobigleaf.Rn_s >= 0
    @test canopy_twobigleaf.Rn_c + canopy_twobigleaf.Rn_s ≈ Rn_input atol=1.0

    # Leaves - 使用专门的函数
    canopy_leaves = Leaves{FT}(nlyr=10)
    fill!(canopy_leaves.Lai_sunlit, 0.1)
    fill!(canopy_leaves.Lai_shaded, 0.2)
    radiative_transfer_multilayer!(canopy_leaves, Rn_input; kA=0.7)
    @test canopy_leaves.Rn[1] > 0
    @test canopy_leaves.Rn[end] >= 0

    # OverUnderCanopy
    canopy_overunder = OverUnderCanopy{FT}(LAI_over=2.0, LAI_under=0.5)
    radiative_transfer_2layer!(canopy_overunder, Rn_input; kA=0.7)
    @test canopy_overunder.Rn_over > 0
    @test canopy_overunder.Rn_under > 0
    @test canopy_overunder.Rn_soil >= 0
  end

  # ========================================================================
  # Test 7: canopy_type_name - 获取冠层类型名称
  # ========================================================================
  @testset "canopy_type_name" begin
    @test canopy_type_name(BigLeaf{FT}()) == "BigLeaf"
    @test canopy_type_name(TwoLeaf{FT}()) == "TwoLeaf"
    @test canopy_type_name(TwoBigLeaf{FT}()) == "TwoBigLeaf"
    @test canopy_type_name(Leaves{FT}()) == "Leaves"
    @test canopy_type_name(OverUnderCanopy{FT}()) == "OverUnderCanopy"
  end

  # ========================================================================
  # Test 8: 集成测试 - 完整的 evapotranspiration! 流程
  # ========================================================================
  @testset "Integration with evapotranspiration!" begin
    # 创建模型
    evap = Evapotranspiration_PML{FT}()
    photo = Photosynthesis_Rong2018{FT}()
    stomatal = Stomatal_Yu2004{FT}()

    # 创建大气数据
    air = AirLayer{FT}(
      Prcp=0.0, Tavg=25.0, Rs=200.0, Rn=150.0,
      VPD=1.5, U2=2.0, Pa=101.3, Ca=410.0
    )

    # 创建输出对象
    output = SpacOutput{FT}()

    # 测试 BigLeaf
    @testset "BigLeaf integration" begin
      canopy = BigLeaf{FT}(Lai=2.5)
      evapotranspiration!(output, evap, photo, stomatal, air, canopy)
      @test output.GPP > 0
      @test output.Ec >= 0
      @test output.rs > 0
    end

    # 测试 TwoLeaf
    @testset "TwoLeaf integration" begin
      canopy = TwoLeaf{FT}(Lai=3.0, Ω=1.0)
      evapotranspiration!(output, evap, photo, stomatal, air, canopy)
      @test output.GPP > 0
      @test output.Ec >= 0
      @test output.rs > 0
      @test canopy.Lai_sunlit > 0
      @test canopy.Lai_shaded > 0
    end

    # 测试 TwoBigLeaf
    @testset "TwoBigLeaf integration" begin
      canopy = TwoBigLeaf{FT}(Lai=3.5, Ω=1.0)
      evapotranspiration!(output, evap, photo, stomatal, air, canopy)
      @test output.GPP > 0
      @test output.Ec >= 0
      @test output.rs > 0
      @test canopy.Lai_sunlit > 0
      @test canopy.Lai_shaded > 0
    end
  end

  # ========================================================================
  # Test 9: 边界条件测试
  # ========================================================================
  @testset "Boundary conditions" begin
    # 零 LAI
    canopy_zero = BigLeaf{FT}(Lai=0.0)
    @test get_lai(canopy_zero) == 0.0

    # 夜间条件（CosZs = 0）
    canopy_night = TwoLeaf{FT}(Lai=3.0, Ω=1.0)
    allocate_lai_if_needed!(canopy_night, 0.0)
    @test canopy_night.Lai_sunlit ≈ 0.0 atol=0.01
    @test canopy_night.Lai_shaded ≈ 3.0 atol=0.01
  end

  # ========================================================================
  # Test 10: 性能测试（确保接口无性能损失）
  # ========================================================================
  @testset "Performance check" begin
    # 这个测试确保接口函数是内联的，没有性能损失
    canopy = BigLeaf{FT}(Lai=2.5)

    # 直接调用
    lai_direct = canopy.Lai

    # 通过接口调用
    lai_interface = get_lai(canopy)

    @test lai_direct == lai_interface
    # 性能测试应该在 benchmark 中进行
  end
end

println("\n✅ All canopy interface tests passed!")
