@testset "stomatal_conductance" begin
  FT = Float64
  Ag = 10. # [umol m-2 s-1]
  Rd = 0.15Ag
  VPD = 2.0 # kPa
  Ca = 380. # ppm, [umol-1 mol]

  Tavg = 25.0  # degC
  Pa = 101.325 # kPa

  stomatal_Yu2004 = Stomatal_Yu2004{FT}(D0=0.7, g1=10)
  stomatal_Medlyn2011 = Stomatal_Medlyn2011{FT}(; g0=0.0001, g1=2.0)

  # Gs在0.1 ~ 0.4 [mol m-2 s-1]是合理的值域，对应的阻力在 100~400 [s m-1]
  gs_yu2004 = stomatal_conductance(stomatal_Yu2004, Ag, Rd, VPD, Ca)
  gs_medlyn2011 = stomatal_conductance(stomatal_Medlyn2011, Ag, Rd, VPD, Ca)

  @test gs_yu2004 ≈ 0.10916179337231968
  @test gs_medlyn2011 ≈ 0.08650343275861604

  # 1 [mol m⁻² s⁻¹] ≈ 0.0245 [m s-1]
  @test mol2m(1., Tavg, Pa) ≈ 0.02446540217236615
  @test mol2m_rong2018(1., Tavg, Pa) ≈ 0.024468739156082532
end
