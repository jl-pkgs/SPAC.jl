using SPAC, Test, Ipaper, DataFrames, Statistics
using Printf

println("\n" * "="^70)
println("  不同冠层方案 ET/GPP 模拟精度评估")
println("="^70)

FT = Float64  # Local variable

# 模型组件
evap = Evapotranspiration_PML{Float64}()
photo = Photosynthesis_Rong2018{Float64}()
stomatal = Stomatal_Yu2004{Float64}()

# 加载数据
df_out, df, _par = deserialize(file_FLUXNET_CRO_USTwt)
(; Prcp, LAI, Rn, Rs, Tavg, U2, VPD, Ca, Pa) = df
GPP_obs, ET_obs = df.GPPobs, df.ETobs

println("\n📊 数据概况: US-Twt (CRO)")
println("   样本数: $(length(df.date)), LAI: $(round(minimum(LAI), digits=1))-$(round(maximum(LAI), digits=1))")

#==============================================================================#
#  模拟函数：支持不同冠层方案
#==============================================================================#
function simulate_canopy(evap, photo, stomatal, df; canopy_type=:BigLeaf, Ω=1.0)
  (; Prcp, LAI, Rn, Rs, Tavg, U2, VPD, Ca, Pa) = df
  ntime = length(Prcp)

  GPP = zeros(Float64, ntime)
  ET = zeros(Float64, ntime)
  Ec = zeros(Float64, ntime)
  Es = zeros(Float64, ntime)
  Ei = zeros(Float64, ntime)

  air = AirLayer{Float64}()
  output = SpacOutput{Float64}()
  frame = 3  # 8天滑动窗口

  # 第一遍：计算 Pi 和 Es_eq
  Pi_arr = zeros(Float64, ntime)
  Es_eq_arr = zeros(Float64, ntime)

  @inbounds for i in 1:ntime
    # 创建对应的冠层对象
    canopy = if canopy_type == :BigLeaf
      BigLeaf{FT}(Lai=LAI[i])
    elseif canopy_type == :TwoLeaf
      TwoLeaf{FT}(Lai=LAI[i], Ω=Ω)
    elseif canopy_type == :TwoBigLeaf
      TwoBigLeaf{FT}(Lai=LAI[i], Ω=Ω)
    end

    SPAC.update!(air, Prcp[i], Tavg[i], Rs[i], Rn[i], VPD[i], U2[i], Pa[i]; Ca=Ca[i])
    evapotranspiration!(output, evap, photo, stomatal, air, canopy)

    GPP[i] = output.GPP
    Ec[i] = output.Ec
    Ei[i] = output.Ei
    Pi_arr[i] = output.Pi
    Es_eq_arr[i] = output.Es_eq
  end

  # 土壤蒸发（使用滑动窗口）
  β_Es = SPAC.movmean2(Pi_arr, frame, 0) ./ SPAC.movmean2(Es_eq_arr, frame, 0)
  clamp!(β_Es, 0, 1)
  Es .= β_Es .* Es_eq_arr
  ET .= Ec .+ Ei .+ Es

  DataFrame(; GPP, ET, Ec, Es, Ei)
end

#==============================================================================#
#  运行不同冠层方案模拟
#==============================================================================#
println("\n⏳ 运行模拟...")

r_big = simulate_canopy(evap, photo, stomatal, df; canopy_type=:BigLeaf)
r_two = simulate_canopy(evap, photo, stomatal, df; canopy_type=:TwoLeaf, Ω=0.9)
r_twobig = simulate_canopy(evap, photo, stomatal, df; canopy_type=:TwoBigLeaf, Ω=0.9)

println("   ✓ BigLeaf, TwoLeaf, TwoBigLeaf 完成")

#==============================================================================#
#  精度评估
#==============================================================================#
println("\n" * "="^70)
println("  📈 精度评估结果")
println("="^70)

# GPP 精度
gof_gpp = [
  ("BigLeaf", GOF(GPP_obs, r_big.GPP)),
  ("TwoLeaf(Ω=0.9)", GOF(GPP_obs, r_two.GPP)),
  ("TwoBigLeaf(Ω=0.9)", GOF(GPP_obs, r_twobig.GPP))
]

println("\n▸ GPP 模拟精度 (gC m⁻² d⁻¹):")
println("  " * "-"^66)
println(@sprintf("  %-18s %6s  %6s  %6s  %6s  %6s",
  "Model", "R²", "NSE", "RMSE", "MAE", "Bias"))
println("  " * "-"^66)
for (name, g) in gof_gpp
  println(@sprintf("  %-18s %6.3f  %6.3f  %6.3f  %6.3f  %6.3f",
    name, g.R2, g.NSE, g.RMSE, g.MAE, g.bias))
end

# ET 精度
gof_et = [
  ("BigLeaf", GOF(ET_obs, r_big.ET)),
  ("TwoLeaf(Ω=0.9)", GOF(ET_obs, r_two.ET)),
  ("TwoBigLeaf(Ω=0.9)", GOF(ET_obs, r_twobig.ET))
]

println("\n▸ ET 模拟精度 (mm d⁻¹):")
println("  " * "-"^66)
println(@sprintf("  %-18s %6s  %6s  %6s  %6s  %6s",
  "Model", "R²", "NSE", "RMSE", "MAE", "Bias"))
println("  " * "-"^66)
for (name, g) in gof_et
  println(@sprintf("  %-18s %6.3f  %6.3f  %6.3f  %6.3f  %6.3f",
    name, g.R2, g.NSE, g.RMSE, g.MAE, g.bias))
end

#==============================================================================#
#  测试断言
#==============================================================================#
println("\n" * "="^70)
println("  🧪 单元测试")
println("="^70)

@testset "冠层方案精度测试" begin
  # BigLeaf 基准测试
  mask = r_big.Es .> 0
  @test GOF(df_out.Es_eq[mask], r_big.Es[mask]).MAE <= 0.5

  # GPP 模拟 R² > 0.5
  @test gof_gpp[1][2].R2 > 0.5  # BigLeaf
  @test gof_gpp[2][2].R2 > 0.5  # TwoLeaf
  @test gof_gpp[3][2].R2 > 0.5  # TwoBigLeaf

  # ET 模拟：BigLeaf NSE > 0，其他方案 R² > 0.3
  @test gof_et[1][2].NSE > 0  # BigLeaf
  @test gof_et[2][2].R2 > 0.3  # TwoLeaf
  @test gof_et[3][2].R2 > 0.3  # TwoBigLeaf

  # 物理一致性
  @test mean(abs.(r_big.GPP .- r_two.GPP)) / mean(r_big.GPP) < 0.3
  @test mean(abs.(r_two.GPP .- r_twobig.GPP)) / mean(r_two.GPP) < 0.1
end

#==============================================================================#
#  综合评价
#==============================================================================#
println("\n" * "="^70)
println("  📋 综合评价")
println("="^70)

best_gpp_idx = argmax([g[2].NSE for g in gof_gpp])
best_et_idx = argmax([g[2].NSE for g in gof_et])

println(@sprintf("\n  🏆 GPP 最佳方案: %-15s (NSE=%.3f)",
  gof_gpp[best_gpp_idx][1], gof_gpp[best_gpp_idx][2].NSE))
println(@sprintf("  🏆 ET  最佳方案: %-15s (NSE=%.3f)",
  gof_et[best_et_idx][1], gof_et[best_et_idx][2].NSE))

println("\n  💡 结论:")
println("     - BigLeaf: 计算简单，适合大尺度应用")
println("     - TwoLeaf: 区分向阳叶/背阴叶，更精确的光合计算")
println("     - TwoBigLeaf: CLM风格导度积分，适合复杂冠层")
println()
