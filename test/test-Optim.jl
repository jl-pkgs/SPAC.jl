using SPAC

begin
  FT = Float64
  evap = Evapotranspiration_PML{FT}()
  photo = Photosynthesis_Rong2018{FT}()
  stomatal = Stomatal_Yu2004{FT}()
  model = LandModel{FT}(evap, photo, stomatal)
  println(model)


  df_out, df, _par = deserialize(file_FLUXNET_CRO_USTwt)
  model_gof(model, df)

  params = Params(model)
  params[params.name.==:hc, :bound] = [(0.1, 2.0)] # 修改参数bounds
end

@testset "model_goal" begin
  parNames = [:kQ, :VCmax25, :VPDmin]
  parValues = [0.6, 10., 0.8]
  @test model_goal(parValues, model, df; parNames, params) <= -0.15
end


@testset "Optims" begin
  parNames = setdiff(params.name, [:d_PC])
  d_par = default_params(model; parNames)

  theta, gof = optim(model, df; parNames, params, maxn=5000, fun_gof=of_KGE)
  # DataFrame(; name=d_par.name, default=d_par.value, optim=theta)
  @test gof.NSE[1] >= 0.54
  @test gof.NSE[2] >= 0.60

  theta, gof = optim(model, df; parNames, params, maxn=1500, fun_gof=of_NSE)
  @test gof.NSE[1] >= 0.50
  @test gof.NSE[2] >= 0.60
end
