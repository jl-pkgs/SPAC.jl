export Params, AbstractModel, update!


import FieldMetadata: @bounds, bounds, @units, units

abstract type AbstractModel{FT} end

function unlist(list::Vector)
  res = []
  for x in list
    isa(x, Vector) ? append!(res, unlist(x)) : push!(res, x)
  end
  return map(x -> x, res)
end


_fieldname(x, field::Symbol) = field

"递归获取模型参数属性"
function get_ModelParamRecur(x::T, path=Vector{Symbol}(); fun=bounds) where {FT,T<:AbstractModel{FT}}
  fileds = fieldnames(T) |> collect
  map(field -> begin
      value = getfield(x, field)
      _path = [path..., field]
      r = fun(x, field)
      TYPE = typeof(value)
      if TYPE == FT 
        return (_path => r) 
      elseif TYPE <: AbstractModel
        get_ModelParamRecur(value, _path; fun)
      else
        return []
      end
    end, fileds) |> unlist
end


function Params(model::AbstractModel)
  _names = get_ModelParamRecur(model; fun=_fieldname)
  _bounds = get_ModelParamRecur(model; fun=bounds)
  _values = get_ModelParamRecur(model; fun=getfield)
  _units = get_ModelParamRecur(model; fun=units)

  # 返回NamedTuple向量，符合Tables接口
  [(name=last(n), value=last(v), bound=last(b),
    unit=last(u), path=first(n))
   for (n, v, u, b) in zip(_names, _values, _units, _bounds)] |> DataFrame
end


function update!(model::AbstractModel{FT}, names::Vector{Symbol}, values::Vector{FT},
  ; params::Union{Nothing,DataFrame}=nothing) where {FT}
  isnothing(params) && (params = ModelParams(model))

  for (name, value) in zip(names, values)
    rows = filter(row -> row.name == name, params)
    @assert size(rows, 1) == 1 "Duplicated parameters are not allowed!"
    update!(model, rows.path[1], value)
  end
end

function update!(model::AbstractModel{FT}, path::Vector{Symbol}, value::FT) where {FT}
  if length(path) == 1
    setfield!(model, path[1], value)
  elseif length(path) > 1
    submodel = getfield(model, path[1])
    update!(submodel, path[2:end], value)
  end
end
