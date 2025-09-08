using SPAC

begin
  FT = Float64
  evap = Evapotranspiration_PML{FT}()
  photo = Photosynthesis_Rong2018{FT}()
  stomatal = Stomatal_Yu2004{FT}()

  df_out, df, _par = deserialize(file_FLUXNET_CRO_USTwt)
  model = LandModel{FT}(evap, photo, stomatal)
  model_gof(model, df)
end

@testset "model_goal" begin
  parnames = [:kQ, :VCmax25, :VPDmin]
  parvalues = [0.6, 10., 0.8]
  @test model_goal(parvalues, model, df; parnames, params) <= -0.15
end


@testset "Optims" begin
  params = Params(model) |> DataFrame
  params[params.name.==:hc, :bound] = [(0.1, 2.0)] # 修改参数bounds

  parnames = setdiff(params.name, [:d_PC])
  d_par = default_params(model; parnames)

  theta, gof = optim(model, df; parnames, params, maxn=1000, fun_gof=of_KGE)
  # DataFrame(; name=d_par.name, default=d_par.value, optim=theta)
  @test gof.NSE[1] >= 0.50
  @test gof.NSE[2] >= 0.55

  theta, gof = optim(model, df; parnames, params, maxn=1000, fun_gof=of_NSE)
  @test gof.NSE[1] >= 0.50
  @test gof.NSE[2] >= 0.60
end
