pacman::p_load(
  Ipaper, data.table, dplyr, lubridate, 
  gg.layers
)

df_org = fread("OUTPUT/PMLV2China_flux37_LAI_whit,ConstPC_OUTPUT.csv")
df_new = fread("Project_PhotoPeriod/OUTPUT/PMLV2China_flux37_LAI_whit,ConstPC_OUTPUT.csv")

list(
  org = df_org[, GOF(GPP_obs, GPP)],
  new = df_new[, GOF(GPP_obs, GPP)]
)

list(
  org = df_org[, GOF(ET_obs, ET)],
  new = df_new[, GOF(ET_obs, ET)]
)
