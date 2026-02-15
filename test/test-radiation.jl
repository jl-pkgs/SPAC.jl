using SPAC, Test


@testset "Norman_Longwave" begin
  n = 50
  ϵ = [1.0; fill(0.98, n - 1)]
  T_leaf = [20.0; fill(25.0, n - 1)]
  τd = ones(n) .* 0.915     # transmittance of diffuse radiation through each leaf layer
  L_up, L_dn, Rln, Rln_soil, Rln_veg = Norman_Longwave(T_leaf, ϵ, τd)
  @test Rln_soil ≈ 28.382474037679856
  @test Rln_veg ≈ -75.5489070386997
end


@testset "Beer sunlit/shaded partition" begin
  Lai_sunlit = 1.2
  Lai_shaded = 2.8
  Rs = 600.0
  Rs_sun, Rs_sha = partition_sunshade_radiation_beer(Lai_sunlit, Lai_shaded, Rs; f_dif=0.2, enable_partition=true)

  @test Rs_sun > Rs_sha

  weighted_flux = Rs_sun * Lai_sunlit + Rs_sha * Lai_shaded
  @test weighted_flux ≈ Rs rtol=1e-6

  # 零 LAI 边界
  Rs_sun0, Rs_sha0 = partition_sunshade_radiation_beer(0.0, 0.0, Rs)
  @test Rs_sun0 == 0.0
  @test Rs_sha0 == 0.0
end

@testset "Norman_Shortwave" begin
  # 从下到上
  dLAI = [0.1, 0.2, 0.3] |> reverse
  Kb = 1.0
  Ω = 1.0
  ## 
  fsun = cal_fsun(vcat(NaN, dLAI), Kb, Ω)
  @test fsun[2:end] == [0.6376281516217733, 0.8187307530779818, 0.9512294245007141]
    
  PAR_sun, PAR_sha, frac_sha, frac_sun = Norman_Shortwave(dLAI)
  @test PAR_sha ≈ [168.89332828782798, 159.20858966542008, 143.23329123213492]
end
