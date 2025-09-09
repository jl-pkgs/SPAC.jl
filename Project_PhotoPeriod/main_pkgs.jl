using RTableTools, Ipaper, SPAC
using Plots
import Base: NamedTuple
using DataFrames: Cols
NamedTuple(names::AbstractVector, values::AbstractVector) =
  NamedTuple{tuple(names...)}(values)


FT = Float64
function build_model(; FT=Float64,
  evap=Evapotranspiration_PML{FT}(),
  photo=Photosynthesis_Rong2018{FT}(),
  stomatal=Stomatal_Yu2004{FT}(),
  PC_photo::Bool=false, PC_g0::Bool=false, PC_g1::Bool=false
)
  model = LandModel{FT}(evap, photo, stomatal)
  model.photo.PC_photo = PC_photo
  model.stomatal.PC_g1 = PC_g1
  if typeof(model.stomatal) <: Stomatal_Medlyn2011
    model.stomatal.PC_g0 = PC_g0
  end
  model
end


function run_model(; site="", maxn=5000, of_gof=:NSE,
  optim_PC=false,
  PC_photo::Bool=false, PC_g0::Bool=false, PC_g1::Bool=false,
  stomatal=Stomatal_Yu2004,
  type_lai="whit", kw...)

  model = build_model(; stomatal=stomatal{Float64}(), PC_photo, PC_g0, PC_g1)
  params = Params(model)

  vars_rm = optim_PC ? Vector{Symbol}() : [:d_PC]
  parNames = setdiff(params.name, vars_rm)

  ## run Models
  d = df[df.name.==site, :]
  (; date, GPP, ET, prcp, LAI_whit, LAI_glass, Rnl, Rns, Rs, Tavg, U2, VPD, Ca, Pa, PC) = d

  Rn = Rns + Rnl
  forcing = (; date, GPP_obs=GPP, ET_obs=ET,
    Prcp=prcp, Tavg, Rs, Rn, VPD, U2,
    PC, Pa, Ca) |> DataFrame
  type_lai == "whit" && (forcing.LAI = LAI_whit)
  type_lai == "glass" && (forcing.LAI = LAI_glass)

  ## 模型参数率定
  theta, gof = optim(model, forcing; parNames, params, var_obs=[:ET_obs, :GPP_obs], maxn)

  SPAC.update!(model, parNames, theta; params=Params(model)) # updating params
  out = evapotranspiration(model, forcing) |> DataFrame
  cbind(out; GPP_obs=GPP, ET_obs=ET)

  gof = [
    (; site, kw..., var="ET", GOF(ET, out.ET)...),
    (; site, kw..., var="GPP", GOF(GPP, out.GPP)...)] |> DataFrame
  (; out=out[:, Cols(:GPP_obs, :ET_obs, 1:end)],
    par=NamedTuple(parNames, theta), gof)
end


function get_prefix(; PC_photo, optim_PC, type_lai)
  prefix = PC_photo ? "ConstPC" : "NonPC"
  if PC_photo && optim_PC
    prefix = "WithPC"
  end
  "LAI_$type_lai" * ",$prefix"
end


function process(; optim_PC=true,
  PC_photo::Bool=false, PC_g0::Bool=false, PC_g1::Bool=false,
  type_lai="glass",
  outdir, stomatal, kw...)

  mkpath(outdir)

  prefix = get_prefix(; PC_photo, optim_PC, type_lai)
  N = length(sites)
  res = Vector{Any}(undef, N)
  @par for i in 1:N
    site = sites[i]
    printstyled("[$i] $site\n", color=:green, bold=true)
    res[i] = run_model(; site, optim_PC,
      PC_photo, PC_g0, PC_g1,
      type_lai, of_gof=:NSE, stomatal, kw...)
  end

  df_gof = vcat(map(x -> x.gof, res)...)
  df_out = vcat(map(x -> x.out, res)...)

  mat_par = cat(map(x -> collect(x.par), res)..., dims=2) |> transpose |> collect
  parNames = collect(keys(res[1].par))
  df_par = cbind(DataFrame(mat_par, parNames); site=sites)[:, Cols(:site, 1:end)]

  fwrite(df_gof, "$outdir/PMLV2China_flux37_$(prefix)_gof.csv")
  fwrite(df_par, "$outdir/PMLV2China_flux37_$(prefix)_par.csv")
  fwrite(df_out, "$outdir/PMLV2China_flux37_$(prefix)_OUTPUT.csv")

  [(; var="ET", GOF(df_out.ET_obs, df_out.ET)...),
    (; var="GPP", GOF(df_out.GPP_obs, df_out.GPP)...)] |> DataFrame
end
