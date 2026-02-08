using SPAC
using Test
import SPAC: update!

@testset "LandModel-Canopy Integration" begin
  # 1. 默认 BigLeaf
  model1 = LandModel{FT}(
    evap = Evapotranspiration_PML{FT}(),
    photo = Photosynthesis_Rong2018{FT}(),
    stomatal = Stomatal_Yu2004{FT}()
  )
  @test model1.canopy isa BigLeaf{FT}
  
  # 2. 自定义 TwoLeaf
  model2 = LandModel{FT}(
    evap = Evapotranspiration_PML{FT}(),
    photo = Photosynthesis_Rong2018{FT}(),
    stomatal = Stomatal_Yu2004{FT}(),
    canopy = TwoLeaf{FT}(Lai=3.5, Ω=0.9)
  )
  @test model2.canopy isa TwoLeaf{FT}
  @test model2.canopy.Lai == 3.5
  @test model2.canopy.Ω == 0.9
  
  # 3. 参数收集（Lai 不是参数，是逐时刻输入）
  params = Params(model2)
  @test !any(params.name .== :Lai)  # Lai 不应该在参数列表中
  @test any(params.name .== :Ω)  # Ω 是可优化参数
  
  # 4. 参数更新（只更新 Ω，不更新 Lai）
  update!(model2, [:Ω], [0.8]; params=params)
  @test model2.canopy.Ω == 0.8
  # Lai 可以直接设置
  model2.canopy.Lai = 4.0
  @test model2.canopy.Lai == 4.0
  
  # 5. 多冠层对比
  models = Dict(
    "BigLeaf" => LandModel{FT}(
      evap = Evapotranspiration_PML{FT}(),
      photo = Photosynthesis_Rong2018{FT}(),
      stomatal = Stomatal_Yu2004{FT}(),
      canopy = BigLeaf{FT}()
    ),
    "TwoLeaf" => LandModel{FT}(
      evap = Evapotranspiration_PML{FT}(),
      photo = Photosynthesis_Rong2018{FT}(),
      stomatal = Stomatal_Yu2004{FT}(),
      canopy = TwoLeaf{FT}(Lai=3.5)
    )
  )
  
  for (name, model) in models
    @test model.canopy isa AbstractLeaf{FT}
  end
  
  # 6. 向后兼容性测试（三参数构造）
  model3 = LandModel{FT}(
    Evapotranspiration_PML{FT}(),
    Photosynthesis_Rong2018{FT}(),
    Stomatal_Yu2004{FT}()
  )
  @test model3.canopy isa BigLeaf{FT}
  
  # 7. 测试 evapotranspiration 函数使用 model.canopy
  air = AirLayer{FT}()
  canopy_bare = BigLeaf{FT}(; Lai = 0.0)
  canopy_veg = BigLeaf{FT}(; Lai = 2.0)
  evap = Evapotranspiration_PML{FT}()
  photo = Photosynthesis_Rong2018{FT}()
  stomatal = Stomatal_Yu2004{FT}()
  
  r_bare = evapotranspiration(evap, photo, stomatal, air, canopy_bare)
  r_veg = evapotranspiration(evap, photo, stomatal, air, canopy_veg)
  
  @test r_bare.ET_water == r_veg.ET_water
  @test r_bare[[:GPP, :Ec, :Ecr, :Eca, :Ei, :Pi]] == zeros(6)
  
  # 8. 测试 LandModel 的 Params 包含冠层参数
  model4 = LandModel{FT}(
    evap = Evapotranspiration_PML{FT}(),
    photo = Photosynthesis_Rong2018{FT}(),
    stomatal = Stomatal_Yu2004{FT}(),
    canopy = TwoLeaf{FT}(Lai=3.5, Ω=0.9)
  )
  params4 = Params(model4)
  @test !any(params4.name .== :Lai)  # Lai 不是参数
  @test any(params4.name .== :Ω)  # Ω 是参数
  @test any(params4.name .== :α)  # evap parameter
  @test any(params4.name .== :VCmax25)  # photo parameter
  @test any(params4.name .== :g1)  # stomatal parameter
end
