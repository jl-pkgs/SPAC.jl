include("Canopy.jl")
include("Interface.jl")

include("PET.jl")
include("utilize.jl")

include("water_constrain.jl")
include("photosynthesis.jl")
include("stomatal_conductance.jl")
include("leaf_conductance.jl")

# Include multilayer functions after photosynthesis is defined
include("multilayer.jl")

include("evaporation_interception.jl")
# include("evapotranspiration.jl")


function transpiration()
end

function evaporation()
end
