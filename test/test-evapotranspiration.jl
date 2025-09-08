using SPAC, Test
using Serialization

@testset "evapotranspiration One moment" begin
  FT = Float64

  air = AirLayer{FT}()
  canopy_bare = BigLeaf{FT}(; Lai = 0.0)
  canopy_veg = BigLeaf{FT}(; Lai = 2.0)
  evap = Evapotranspiration_PML{FT}()
  photo = Photosynthesis_Rong2018{FT}()
  stomatal = Stomatal_Yu2004{FT}()
  stomatal = Stomatal_Medlyn2011{FT}()

  r_bare = evapotranspiration(evap, photo, stomatal, air, canopy_bare)
  r_veg = evapotranspiration(evap, photo, stomatal, air, canopy_veg)

  r_bare.ET_water == r_veg.ET_water
  @test r_bare[[:GPP, :Ec, :Ecr, :Eca, :Ei, :Pi]] == zeros(6)

  model = LandModel{FT}(; evap, photo, stomatal)
  Params(model) |> DataFrame
end


## 测试单站点的运行与模拟
# @testset "evapotranspiration One Site" 
begin
  FT = Float64
  evap = Evapotranspiration_PML{FT}()
  photo = Photosynthesis_Rong2018{FT}()
  stomatal = Stomatal_Yu2004{FT}()

  df_out, df, _par = deserialize(file_FLUXNET_CRO_USTwt)
  res = evapotranspiration(evap, photo, stomatal, df) |> DataFrame

  # df_out.ET_obs
  GOF(df.GPPobs, df_out.GPP_sim)
  GOF(df.ETobs, df_out.ET_sim)

  GOF(df.GPPobs, res.GPP)
  GOF(df.ETobs, res.ET)
end
