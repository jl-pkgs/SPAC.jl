export default_params, model_gof, model_goal, optim
export of_KGE, of_NSE, GOF
export DataFrame

import DataFrames: DataFrame
using ModelParams: of_KGE, of_NSE, GOF


function default_params(model; parNames=nothing)
  params = Params(model)
  isnothing(parNames) && (parNames = params.name)
  inds = indexin(parNames, params.name) # target params index
  params[inds, :]
end


"""
- `var_obs`: `[:ETobs, :GPPobs]`
"""
function model_gof(model, df; parValues=nothing, parNames=nothing,
  fun_gof=of_KGE, var_obs=[:ETobs, :GPPobs])

  if !isnothing(parValues) && !isnothing(parNames)
    SPAC.update!(model, parNames, parValues; params=Params(model)) # updating params
  end

  ET_obs = df[:, var_obs[1]]
  GPP_obs = df[:, var_obs[2]]

  res = evapotranspiration(model, df) |> DataFrame
  of_ET = fun_gof(ET_obs, res.ET)
  of_GPP = fun_gof(GPP_obs, res.GPP)
  goal = -(of_GPP + of_ET) / 2

  DataFrame([
    (; type=:GPP, GOF(GPP_obs, res.GPP)...),
    (; type=:ET, GOF(ET_obs, res.ET)...)
  ])
end


function model_goal(theta::AbstractVector{T}, model::LandModel{T}, df::DataFrame;
  fun_gof=of_KGE, var_obs=[:ETobs, :GPPobs],
  parNames::Vector, params::Union{Nothing,DataFrame}=nothing) where {T<:Real}

  SPAC.update!(model, parNames, theta; params)
  res = evapotranspiration(model, df) |> DataFrame

  of_ET = fun_gof(df[:, var_obs[1]], res.ET)
  of_GPP = fun_gof(df[:, var_obs[2]], res.GPP)
  goal = -(of_GPP + of_ET) / 2
  goal
end


"""
# Arguments
- `parNames`: 需要指定，需要率定的模型参数名
- `params`: `name`, `value`, `bound`, `path`
"""
function optim(model::LandModel{T}, df::DataFrame;
  parNames::Vector,
  fun_gof=of_KGE, var_obs=[:ETobs, :GPPobs],
  theta0::Union{Nothing,Vector{T}}=nothing, maxn::Int=1000,
  params::Union{Nothing,DataFrame}=nothing) where {T<:Real}

  isnothing(params) && (params = Params(model))

  inds = indexin(parNames, params.name) # target params index
  lower = map(x -> x[1], params.bound[inds])
  upper = map(x -> x[2], params.bound[inds])

  isnothing(theta0) && (theta0 = params.value[inds])

  _goal(theta) = model_goal(theta, model, df; parNames, params, fun_gof, var_obs)

  theta, feval, exitflag = sceua(_goal, theta0, lower, upper; maxn)
  SPAC.update!(model, params.name, params.value; params) # backup

  gof = model_gof(model, df; parValues=theta, parNames, var_obs)
  theta, gof
end
