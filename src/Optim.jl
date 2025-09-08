export default_params, model_gof, model_goal, optim
export of_KGE, of_NSE, GOF
export DataFrame

import DataFrames: DataFrame
using ModelParams: of_KGE, of_NSE, GOF


function default_params(model; parnames=nothing)
  params = Params(model)
  isnothing(parnames) && (parnames = params.name)
  inds = indexin(parnames, params.name) # target params index
  params[inds, :]
end


function model_gof(model, df; parvalues=nothing, parnames=nothing)
  if !isnothing(parvalues) && !isnothing(parnames)
    SPAC.update!(model, parnames, parvalues; params=Params(model)) # updating params
  end

  res = evapotranspiration(model, df) |> DataFrame
  of_GPP = of_KGE(df.GPPobs, res.GPP)
  of_ET = of_KGE(df.ETobs, res.ET)
  goal = -(of_GPP + of_ET) / 2

  DataFrame([
    (; type=:GPP, GOF(df.GPPobs, res.GPP)...),
    (; type=:ET, GOF(df.ETobs, res.ET)...)
  ])
end


function model_goal(theta::AbstractVector{T}, model::LandModel{T}, df::DataFrame;
  fun_gof=of_KGE,
  parnames::Vector, params::Union{Nothing,DataFrame}=nothing) where {T<:Real}

  SPAC.update!(model, parnames, theta; params)
  res = evapotranspiration(model, df) |> DataFrame

  of_GPP = fun_gof(df.GPPobs, res.GPP)
  of_ET = fun_gof(df.ETobs, res.ET)
  goal = -(of_GPP + of_ET) / 2
  goal
end


"""
# Arguments
- `parnames`: 需要指定，需要率定的模型参数名
- `params`: `name`, `value`, `bound`, `path`
"""
function optim(model::LandModel{T}, df::DataFrame;
  parnames::Vector,
  fun_gof=of_KGE,
  theta0::Union{Nothing,Vector{T}}=nothing, maxn::Int=1000,
  params::Union{Nothing,DataFrame}=nothing) where {T<:Real}

  isnothing(params) && (params = Params(model))

  inds = indexin(parnames, params.name) # target params index
  lower = map(x -> x[1], params.bound[inds])
  upper = map(x -> x[2], params.bound[inds])

  isnothing(theta0) && (theta0 = params.value[inds])

  _goal(theta) = model_goal(theta, model, df; parnames, params, fun_gof)

  theta, feval, exitflag = sceua(_goal, theta0, lower, upper; maxn)
  SPAC.update!(model, params.name, params.value; params) # backup

  gof = model_gof(model, df; parvalues=theta, parnames)
  theta, gof
end
