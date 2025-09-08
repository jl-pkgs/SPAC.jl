export Evapotranspiration_PML
export evapotranspiration, evapotranspiration_multi


# 1t: struct forcing
function evapotranspiration(
  evap::AbstractEvapotranspirationModel{T},
  photo::AbstractPhotosynthesisModel{T},
  stomatal::AbstractStomatalModel{T},
  air::AirLayer{T}, canopy::BigLeaf{T}) where {T<:Real}

  output = SpacOutput{T}()
  evapotranspiration!(output, evap, photo, stomatal, air, canopy)
end


# 1site: DataFrame forcing
function evapotranspiration(
  evap::AbstractEvapotranspirationModel{T},
  photo::AbstractPhotosynthesisModel{T},
  stomatal::AbstractStomatalModel{T},
  d::AbstractDataFrame; kw...) where {T}

  forcing = get_forcing(d, evap)
  evapotranspiration(evap, photo, stomatal, forcing...; kw...)
end

function evapotranspiration(model::LandModel{T}, d::AbstractDataFrame; kw...) where {T}
  (; evap, photo, stomatal) = model
  evapotranspiration(evap, photo, stomatal, d; kw...)
end


# 相同植被类型多个站点一起的计算
function evapotranspiration_multi(
  evap::AbstractEvapotranspirationModel{T},
  photo::AbstractPhotosynthesisModel{T},
  stomatal::AbstractStomatalModel{T},
  df::AbstractDataFrame; kw...) where {T}

  "site" ∉ names(df) && return evapotranspiration(evap, photo, stomatal, df; kw...)
  sites = df.site
  grps = unique(sites)

  res = []
  for grp in grps # 可以加并行
    inds = sites .== grp
    d = df[inds, :]
    r = evapotranspiration(evap, photo, stomatal, d; kw...) |> DataFrame
    push!(res, r)
  end
  vcat(res...)
end
