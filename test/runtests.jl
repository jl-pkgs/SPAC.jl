using SPAC, Test


# include("test-Ipaper.jl")
include("test-PMLV2.jl")
include("test-Optim.jl")
include("test-stomatal_conductance.jl")
include("test-photosynthesis.jl")
include("test-PET.jl")

include("test-radiation.jl")
include("test-evapotranspiration.jl")
include("test-VCmax.jl")

# New canopy type tests
include("test-twoleaf.jl")
include("test-twobigleaf.jl")

# Multilayer canopy tests
include("test-multilayer.jl")

# Integration tests for all canopy types
include("test-integration.jl")
include("test-twoleaf-integration.jl")


@testset "Model Params update!" begin
  FT = Float64
  model = Photosynthesis_Rong2018{FT}()

  params = Params(model)

  parnames = [:kQ, :VCmax25, :VPDmin]
  parvalues = [0.6, 10., 0.8]
  @time SPAC.update!(model, parnames, parvalues; params)

  @test model.kQ == 0.6
  @test model.VCmax25 == 10.
  @test model.watercons.VPDmin == 0.8
end
