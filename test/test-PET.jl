@testset "ET0 models" begin
  Rn = 200.0    # W/m²
  Ta = 25.0     # °C
  VPD = 1.5     # kPa
  Uz = 3.0      # m/s

  ra_15 = aerodynamic_resistance(1.0; z_obs=15.0)  # s/m
  ra_2 = aerodynamic_resistance(1.0; z_obs=2.0)

  @test ra_15 ≈ 267.1155033258225
  @test ra_2 ≈ 207.66407000788683

  @test ET0_PT72(Rn, Ta) ≈ 6.564860193238437
  @test ET0_Penman48(Rn, Ta, VPD, Uz)[1] ≈ 7.9265544809634365
  @test ET0_FAO98(Rn, Ta, VPD, Uz)[1] ≈ 6.928646397419433

  @test ET0_Monteith65(250.0, 25.0, 1.0, 2.0; hc=0.12, z_obs=15.0)[1] ≈ 6.801658192680292
end
