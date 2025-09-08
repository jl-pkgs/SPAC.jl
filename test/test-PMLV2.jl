using SPAC, Test, Ipaper, DataFrames


begin
  FT = Float64
  evap = Evapotranspiration_PML{FT}()
  photo = Photosynthesis_Rong2018{FT}()
  stomatal = Stomatal_Yu2004{FT}()

  df_out, df, _par = deserialize(file_FLUXNET_CRO_USTwt)
  model = LandModel{FT}(evap, photo, stomatal)
  println(model)
  model_gof(model, df)

  params = Params(model)
  params[params.name.==:hc, :bound] = [(0.1, 2.0)] # 修改参数bounds
end

begin
  df_out, df, _par = deserialize(file_FLUXNET_CRO_USTwt)
  (; date, Prcp, LAI, Rn, Rs, Tavg, U2, VPD, Ca, Pa) = df
  forcing = DataFrame(; date, Prcp, LAI, Rn, Rs, Tavg, U2, VPD, Ca, Pa)

  _par = (α=0.03265625, η=0.069296875, g1=9.552734375,
    VCmax25=17.671875, VPDmin=1.21515625, VPDmax=3.5, D0=0.6541015625,
    kQ=0.10114375, kA=0.89921875, S_sls=0.01015625, fER0=0.152734375, hc=0.5)

  parNames = keys(_par) |> collect
  theta = values(_par) |> collect

  SPAC.update!(model, parNames, theta; params=Params(model)) # updating params
  r = evapotranspiration(model, forcing) |> DataFrame
end

# df.GPP_obs = df.GPPobs
# df.ET_obs = df.ETobs
# par = Param_PMLV2(; _par..., hc=0.5)
# r = PMLV2_sites(df; par)

@testset "CHECK PMLV2 RESULT" begin
  # @test GOF(df.Eeq, r.Es_eq).MAE <= 0.01
  @test GOF(df_out.Es_eq, r.Es_eq).MAE <= 0.002
  # @test GOF(df_out.Ga, r.Ga).MAE <= 1e-10
  # @test GOF(df_out.Gc, r.Gc_w).MAE <= 1e-3
  @test GOF(df_out.Ei, r.Ei).MAE <= 1e-8 # Ei passed Test
  @test GOF(df_out.Es, r.Es).MAE <= 0.01
  
  @test GOF(df_out.GPP_sim, r.GPP).MAE <= 1E-8
  
  @test GOF(df_out.ET_sim, r.ET).MAE <= 0.004
  @test GOF(df_out.Ec, r.Ec).MAE <= 0.004
  @test GOF(df_out.Ecr, r.Ecr).MAE <= 0.01 # 这两处有较大误差
  @test GOF(df_out.Eca, r.Eca).MAE <= 0.01
end

# DataFrame(; previous=df_out.Es_eq, current = r.Es_eq)

# @testset "PMLV2 scalar" begin
#   Prcp = 2.0 # mm
#   Tavg = 20.0
#   Rs = 200.0
#   Rn = 50.0
#   VPD = 2.0
#   U2 = 2.0
#   LAI = 2.0
#   @test_nowarn PMLV2(Prcp, Tavg, Rs, Rn, VPD, U2, LAI;)
# end
