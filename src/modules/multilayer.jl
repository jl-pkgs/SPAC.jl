# Multilayer canopy functions
# These must be included after photosynthesis.jl because they depend on Photosynthesis_Rong2018

include("../Canopy/multilayer_init.jl")
include("../Canopy/multilayer_photosynthesis.jl")
include("../Canopy/multilayer_conductance.jl")
include("../Canopy/multilayer_radiation.jl")
