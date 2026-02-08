using Test
using SPAC

@testset "Ipaper" begin
  # test nanmean2
  @test nanmean2(1, NaN) == 1
  @test nanmean2(NaN, 1) == 1
  r = nanmean2.([1, 2, 3], [1, 4, NaN])
  @test r == [1.0, 3.0, 3]
  @test eltype(r) == Float64
end
