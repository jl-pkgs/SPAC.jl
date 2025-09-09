include("main_pkgs.jl")

function runtests(config; prefix="Medlyn2011")
  for i in eachindex(config)
    kw = config[i]
    (; PC_photo, PC_g1, PC_g0) = kw
    _prefix_PC = "(PC_photo=$PC_photo,g1=$PC_g1,g0=$PC_g0)"
    outdir = "./Project_PhotoPeriod/OUTPUT/$(prefix)_$(_prefix_PC)"
    println(outdir)
    r = process(; type_lai="whit", optim_PC=false, outdir, kw...)
  end
end

f = "Z:/Researches/ET_ModelDev/data-raw/backup/Forcing_PMLV2_China_8day_2003-2022_flux37_v20250108.csv"
df = fread(f) |> replace_missing!

sites = df.name |> unique_sort
site = "三江源"

# build_model(; PC_photo=true, PC_g1=true, PC_g0=true)
# r_PC = run_model(; site, PC_photo=true)
# r_NonPC = run_model(; site, PC_photo=false)
# run_model(; site, PC_photo=true, type_lai="whit", of_gof=:NSE)#.gof
# run_model(; site, PC_photo=false, type_lai="whit", of_gof=:NSE).gof

types_lai = ["whit", "glass"]

## 方案1: Stomatal_Medlyn2011
stomatal = Stomatal_Medlyn2011
configs = [
  (; stomatal, PC_photo=false, PC_g1=false, PC_g0=false), # 对照组
  (; stomatal, PC_photo=true, PC_g1=true, PC_g0=true),    # A1
  (; stomatal, PC_photo=true, PC_g1=true, PC_g0=false),
  (; stomatal, PC_photo=true, PC_g1=false, PC_g0=true),
]
runtests(configs)


## 方案2: 
stomatal = Stomatal_Yu2004
configs = [
  (; stomatal, PC_photo=false, PC_g1=false, PC_g0=false), # 对照组
  (; stomatal, PC_photo=true, PC_g1=true, PC_g0=false),   # A1
  (; stomatal, PC_photo=true, PC_g1=false, PC_g0=false),
]
runtests(configs; prefix="Yu2004")
