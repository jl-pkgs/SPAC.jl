using RTableTools, Ipaper, SPAC
using Plots
import Base: NamedTuple
using DataFrames: Cols
NamedTuple(names::AbstractVector, values::AbstractVector) = 
  NamedTuple{tuple(names...)}(values)

# include("../examples/main_pkgs.jl")
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
begin
  FT = Float64
  evap = Evapotranspiration_PML{FT}()
  photo = Photosynthesis_Rong2018{FT}()
  stomatal = Stomatal_Yu2004{FT}()

  model = LandModel{FT}(evap, photo, stomatal)
  model.photo.use_PC = true

  params = Params(model)
  parnames = setdiff(params.name, [:d_PC])
end

# site, var, NSE, R2, KGE, R, RMSE, MAE, bias, bias_perc, n_valid
# 三江源, ET, 0.6314213848002423, 0.6717941683525882, 0.7677636119861359, 0.8196305072119926, 0.7201708085099746, 0.5147756659383849, -0.17945556206709948, -13.81443029817728, 229
# 三江源, GPP, 0.7512404633610329, 0.7542890674960799, 0.7805236341912695, 0.8684981678138878, 1.0984843194102694, 0.6990723794209167, 0.10099157842344733, 6.744850722759586, 229

function run_model(; site="",
  use_PC=true, type_lai="whit", 
  of_gof=:NSE, kw...)

  d = df[df.name.==site, :]
  (; date, GPP, ET, prcp, LAI_whit, LAI_glass, Rnl, Rns, Rs, Tavg, U2, VPD, Ca, Pa, PC) = d

  Rn = Rns + Rnl
  forcing = (; date, GPP_obs=GPP, ET_obs=ET,
    Prcp=prcp, Tavg, Rs, Rn, VPD, U2,
    PC,
    # LAI=LAI_glass, 
    # LAI=LAI_whit,
    Pa, Ca) |> DataFrame
  use_PC && (forcing.PC = PC)
  type_lai == "whit" && (forcing.LAI = LAI_whit)
  type_lai == "glass" && (forcing.LAI = LAI_glass)

  ## 模型参数率定
  model.photo.use_PC = use_PC
  theta, gof = optim(model, forcing; parNames, params, var_obs=[:ET_obs, :GPP_obs], maxn=5000)
  # model.photo.use_PC = false
  # theta, gof = optim(model, forcing; parNames, params, var_obs=[:ET_obs, :GPP_obs], maxn=5000)
  
  SPAC.update!(model, parNames, theta; params=Params(model)) # updating params
  out = evapotranspiration(model, forcing) |> DataFrame
  cbind(out; GPP_obs=GPP, ET_obs=ET)

  gof = [
    (; site, kw..., var="ET", GOF(ET, out.ET)...),
    (; site, kw..., var="GPP", GOF(GPP, out.GPP)...)] |> DataFrame
  (; out=out[:, Cols(:GPP_obs, :ET_obs, 1:end)], 
    par = NamedTuple(parNames, theta), gof)
end

# r_PC = run_model(; site, use_PC=true)
# r_NonPC = run_model(; site, use_PC=false)

# run_model(; site, use_PC=true, type_lai="whit", of_gof=:NSE)#.gof
# run_model(; site, use_PC=false, type_lai="whit", of_gof=:NSE).gof
function get_prefix(use_PC, type_lai)
  prefix = use_PC ? "WithPC" : "NonPC"
  if use_PC && :d_pc ∉ parNames
    prefix = "ConstPC"
  end
  "LAI_$type_lai" * ",$prefix"
end

function process(; use_PC=false, type_lai="glass")
  prefix = get_prefix(use_PC, type_lai)
  N = length(sites)
  res = Vector{Any}(undef, N)
  @par for i in 1:N
    site = sites[i]
    printstyled("[$i] $site\n", color=:green, bold=true)
    res[i] = run_model(; site, use_PC, type_lai, of_gof=:NSE)
  end

  df_gof = vcat(map(x -> x.gof, res)...)
  df_out = vcat(map(x -> x.out, res)...)
  
  mat_par = cat(map(x -> collect(x.par), res)..., dims=2) |> transpose |> collect
  df_par = cbind(DataFrame(mat_par, parNames); site = sites)[:, Cols(:site, 1:end)]

  outdir = "./Project_PhotoPeriod/OUTPUT"
  mkpath(outdir)
  
  fwrite(df_gof, "$outdir/PMLV2China_flux37_$(prefix)_gof.csv")
  fwrite(df_par, "$outdir/PMLV2China_flux37_$(prefix)_par.csv")
  fwrite(df_out, "$outdir/PMLV2China_flux37_$(prefix)_OUTPUT.csv")

  [(; var="ET", GOF(df_out.ET_obs, df_out.ET)...),
   (; var="GPP", GOF(df_out.GPP_obs, df_out.GPP)...)] |> DataFrame
end




begin
  r_whit_NonPC = process(; type_lai="whit", use_PC=false)
  r_whit_WithPC = process(; type_lai="whit", use_PC=true)
  r_glass_NonPC = process(; type_lai="glass", use_PC=false)
  r_glass_WithPC = process(; type_lai="glass", use_PC=true)

  cbind(r_glass_NonPC; type_lai="GLASS", type_pc="NonPC")
  cbind(r_glass_WithPC; type_lai="GLASS", type_pc="WithPC")

  cbind(r_whit_NonPC; type_lai="WHIT", type_pc="NonPC")
  cbind(r_whit_WithPC; type_lai="WHIT", type_pc="WithPC")

  R = vcat(r_whit_NonPC, r_whit_WithPC, r_glass_NonPC, r_glass_WithPC)[:, Cols(:type_lai, :type_pc, 1:end)]
end


# # LAI_whit & Non_PC
# NSE = 0.4901658678314159, R2 = 0.5713317663570666, KGE = 0.7435630375012554

# # LAI_whit & With_PC
# NSE = 0.5080601316516402, R2 = 0.5889340168907004, KGE = 0.7519358746627769

# # LAI_glass & Non_PC
# NSE = 0.4935561078799632, R2 = 0.5797932490279245, KGE = 0.7532742010794355

# # LAI_glass & With_PC
# NSE = 0.5343449820736267, R2 = 0.6058902430946861, KGE = 0.7663923308482895
