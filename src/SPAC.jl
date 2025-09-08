module SPAC

export PMLV2, PMLV2_sites,
  photosynthesis, cal_Ei_Dijk2021,
  T_adjust_Vm25, f_VPD_Zhang2019

export file_FLUXNET_CRO, file_FLUXNET_CRO_USTwt
export DataFrame, GOF, sceua
export fread, fwrite, melt_list, deserialize


using UnPack, Printf
using Parameters, FieldMetadata
using Statistics
using DataFrames, RTableTools
using DocStringExtensions


import ModelParams: GOF, sceua
import Ipaper: par_map, deserialize


## global data
dir_proj = "$(@__DIR__)/.."
file_FLUXNET_CRO = "$dir_proj/data/CRO/FLUXNET_CRO" |> abspath
file_FLUXNET_CRO_USTwt = "$dir_proj/data/CRO/FLUXNET_CRO_US-Twt" |> abspath

include("ModelParam.jl")
include("SpacOutput.jl")

include("modules/modules.jl")
include("evapotranspiration.jl")
include("ET_PMLV2.jl")

include("Optim.jl")

include("Radiation/tridiagonal_solver.jl")
include("Radiation/helper.jl")
include("Radiation/Norman_Longwave.jl")
include("Radiation/Norman_Shortwave.jl")

include("tools_Ipaper.jl")
# include("ModelCalib.jl")

end # module SPACmodels
