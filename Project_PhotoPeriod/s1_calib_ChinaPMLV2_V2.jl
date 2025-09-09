using RTableTools, Ipaper, SPAC
using Plots
import Base: NamedTuple
using DataFrames: Cols
NamedTuple(names::AbstractVector, values::AbstractVector) =
  NamedTuple{tuple(names...)}(values)


f = "Z:/Researches/ET_ModelDev/data-raw/backup/Forcing_PMLV2_China_8day_2003-2022_flux37_v20250108.csv"
df = fread(f) |> replace_missing!

sites = df.name |> unique_sort
site = "三江源"
# sites = setdiff(sites, ["元江"])
# site = "那曲"
# site = "哀牢山"
# d = df[df.name.==site, :]

# parNames = [
#   :α, :η, :g1, :VCmax25, :VPDmin, :VPDmax, :D0, :kQ, :kA, :S_sls, :fER0#, :d_pc # :hc
# ]

FT = Float64
function build_model(; FT=Float64, 
  evap = Evapotranspiration_PML{FT}(),
  photo=Photosynthesis_Rong2018{FT}(),
  stomatal=Stomatal_Yu2004{FT}())
  LandModel{FT}(evap, photo, stomatal)
end

# site = "哀牢山"
# out, par, gof1_WithPC = run_model(; site, PC_photo=true, optim_PC=true)
# out, par, gof1_NonPC = run_model(; site, PC_photo=false)
# out, par, gof2_WithPC = run_model(; site, PC_photo=true, optim_PC=true, stomatal)
# out, par, gof2_NonPC = run_model(; site, PC_photo=false, stomatal)

function run_model(; site="",
  PC_photo=true, optim_PC=false, type_lai="whit",
  of_gof=:NSE, maxn=5000, 
  stomatal=Stomatal_Yu2004{FT}(),
  kw...)

  ## Initial model
  model = build_model(; stomatal)
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
  PC_photo && (forcing.PC = PC)
  type_lai == "whit" && (forcing.LAI = LAI_whit)
  type_lai == "glass" && (forcing.LAI = LAI_glass)

  ## 模型参数率定
  model.photo.PC_photo = PC_photo
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

# r_PC = run_model(; site, PC_photo=true)
# r_NonPC = run_model(; site, PC_photo=false)

# run_model(; site, PC_photo=true, type_lai="whit", of_gof=:NSE)#.gof
# run_model(; site, PC_photo=false, type_lai="whit", of_gof=:NSE).gof
function get_prefix(; PC_photo, optim_PC, type_lai)
  prefix = PC_photo ? "ConstPC" : "NonPC"
  if PC_photo && optim_PC
    prefix = "WithPC"
  end
  "LAI_$type_lai" * ",$prefix"
end

function process(; PC_photo=false, optim_PC=true, type_lai="glass", 
  outdir, stomatal, kw...)
  prefix = get_prefix(; PC_photo, optim_PC, type_lai)
  N = length(sites)
  res = Vector{Any}(undef, N)
  @par for i in 1:N
    site = sites[i]
    printstyled("[$i] $site\n", color=:green, bold=true)
    res[i] = run_model(; site, PC_photo, optim_PC, type_lai, of_gof=:NSE, stomatal, kw...)
  end

  df_gof = vcat(map(x -> x.gof, res)...)
  df_out = vcat(map(x -> x.out, res)...)

  mat_par = cat(map(x -> collect(x.par), res)..., dims=2) |> transpose |> collect
  parNames = collect(keys(res[1].par))
  df_par = cbind(DataFrame(mat_par, parNames); site=sites)[:, Cols(:site, 1:end)]

  mkpath(outdir)

  fwrite(df_gof, "$outdir/PMLV2China_flux37_$(prefix)_gof.csv")
  fwrite(df_par, "$outdir/PMLV2China_flux37_$(prefix)_par.csv")
  fwrite(df_out, "$outdir/PMLV2China_flux37_$(prefix)_OUTPUT.csv")

  [(; var="ET", GOF(df_out.ET_obs, df_out.ET)...),
    (; var="GPP", GOF(df_out.GPP_obs, df_out.GPP)...)] |> DataFrame
end

stomatal_Medlyn2011 = Stomatal_Medlyn2011{Float64}()
stomatal_Yu2024 = Stomatal_Yu2004{Float64}()

# outdir = "./Project_PhotoPeriod/OUTPUT"
outdir = "./Project_PhotoPeriod/OUTPUT/Medlyn2011_V2"
kw = (; outdir, stomatal = Stomatal_Medlyn2011{Float64}())

begin
  r_whit_WithPC = process(; type_lai="whit", PC_photo=true, optim_PC=true, kw...)
  r_whit_ConstPC = process(; type_lai="whit", PC_photo=true, optim_PC=false, kw...)
  r_whit_NonPC = process(; type_lai="whit", PC_photo=false, kw...)

  r_glass_WithPC = process(; type_lai="glass", PC_photo=true, optim_PC=true, kw...)
  r_glass_ConstPC = process(; type_lai="glass", PC_photo=true, optim_PC=false, kw...)
  r_glass_NonPC = process(; type_lai="glass", PC_photo=false, kw...)

  cbind(r_whit_NonPC; type_lai="WHIT", type_pc="NonPC")
  cbind(r_whit_ConstPC; type_lai="WHIT", type_pc="ConstPC")
  cbind(r_whit_WithPC; type_lai="WHIT", type_pc="WithPC")

  cbind(r_glass_NonPC; type_lai="GLASS", type_pc="NonPC")
  cbind(r_glass_ConstPC; type_lai="GLASS", type_pc="ConstPC")
  cbind(r_glass_WithPC; type_lai="GLASS", type_pc="WithPC")

  R = vcat(r_whit_NonPC, r_whit_ConstPC, r_whit_WithPC,
    r_glass_NonPC, r_glass_ConstPC, r_glass_WithPC)[:, Cols(:type_lai, :type_pc, 1:end)]
end

# # LAI_whit & Non_PC
# NSE = 0.4901658678314159, R2 = 0.5713317663570666, KGE = 0.7435630375012554

# # LAI_whit & With_PC
# NSE = 0.5080601316516402, R2 = 0.5889340168907004, KGE = 0.7519358746627769

# # LAI_glass & Non_PC
# NSE = 0.4935561078799632, R2 = 0.5797932490279245, KGE = 0.7532742010794355

# # LAI_glass & With_PC
# NSE = 0.5343449820736267, R2 = 0.6058902430946861, KGE = 0.7663923308482895
