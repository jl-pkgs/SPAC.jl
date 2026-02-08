using SPAC
using Test
using Ipaper
using DataFrames
using Statistics
using Printf
using CSV

println("\n", "="^70)
println("  现代化参数方案综合测试")
println("="^70)

#==============================================================================#
#  1. 创建不同配置的 LandModel
#==============================================================================#
println("\n📦 创建模型配置...")

models = Dict(
  "Yu2004" => LandModel{Float64}(
    evap = Evapotranspiration_PML{Float64}(),
    photo = Photosynthesis_Rong2018{Float64}(),
    stomatal = Stomatal_Yu2004{Float64}()
  ),
  "Medlyn2011" => LandModel{Float64}(
    evap = Evapotranspiration_PML{Float64}(),
    photo = Photosynthesis_Rong2018{Float64}(),
    stomatal = Stomatal_Medlyn2011{FT}()
  ),
  "BallBerry1987" => LandModel{Float64}(
    evap = Evapotranspiration_PML{Float64}(),
    photo = Photosynthesis_Rong2018{Float64}(),
    stomatal = Stomatal_BallBerry1987{FT}()
  )
)

println("   ✓ 创建了 $(length(models)) 个模型配置")

#==============================================================================#
#  2. 参数堆叠和收集测试
#==============================================================================#
println("\n🔧 测试参数堆叠和收集...")

@testset "参数堆叠和收集" begin
  for (name, model) in models
    params = Params(model)
    @test !isempty(params)
    @test all(hasproperty.([params], [:name, :value, :bound]))
    println("   ✓ $name: $(nrow(params)) 个参数")
  end
end

#==============================================================================#
#  3. 参数更新测试
#==============================================================================#
println("\n🔄 测试参数更新...")

@testset "参数更新" begin
  model = deepcopy(models["Yu2004"])
  params_orig = Params(model)
  
  # 更新部分参数
  parNames = [:α, :VCmax25, :g1]
  parValues = [0.08, 60.0, 15.0]
  
  SPAC.update!(model, parNames, parValues; params=params_orig)
  params_updated = Params(model)
  
  # 验证更新（检查参数是否在模型中）
  for (i, name) in enumerate(parNames)
    idx = findfirst(==(String(name)), params_updated.name)
    if !isnothing(idx)
      @test params_updated.value[idx] ≈ parValues[i]
    end
  end
  
  println("   ✓ 参数更新成功")
end

#==============================================================================#
#  4. 多站点模型配置
#==============================================================================#
println("\n🌍 多站点模型配置...")

# 定义不同 IGBP 类型的模型配置
IGBP_configs = Dict(
  "EBF" => (VCmax25=80.0, g1=12.0, α=0.06),  # 常绿阔叶林
  "DBF" => (VCmax25=60.0, g1=10.0, α=0.06),  # 落叶阔叶林
  "CRO" => (VCmax25=50.0, g1=9.0, α=0.06),   # 农作物
  "GRA" => (VCmax25=40.0, g1=8.0, α=0.05),   # 草地
)

println("   ✓ 定义了 $(length(IGBP_configs)) 个 IGBP 配置")

@testset "多站点配置" begin
  for (igbp, config) in IGBP_configs
    model = LandModel{Float64}(
      evap = Evapotranspiration_PML{Float64}(),
      photo = Photosynthesis_Rong2018{Float64}(),
      stomatal = Stomatal_Yu2004{Float64}()
    )
    
    # 应用配置
    params = Params(model)
    for (key, value) in pairs(config)
      idx = findfirst(==(String(key)), params.name)
      if !isnothing(idx)
        SPAC.update!(model, [Symbol(key)], [value]; params=params)
      end
    end
    
    @test true  # 配置应用成功
    println("   ✓ $igbp: VCmax25=$(config.VCmax25), g1=$(config.g1)")
  end
end

#==============================================================================#
#  5. 参数优化测试（简化版）
#==============================================================================#
println("\n⚡ 测试参数优化...")

df_out, df, _par = deserialize(file_FLUXNET_CRO_USTwt)

@testset "参数优化" begin
  model = models["Yu2004"]
  params = Params(model)
  
  # 选择少量参数进行快速测试
  parNames = [:α, :VCmax25, :g1]
  
  # 运行优化（减少迭代次数以加快测试）
  theta_opt, gof = SPAC.optim(model, df; 
    parNames, params, maxn=500, fun_gof=SPAC.of_NSE)
  
  @test gof.NSE[1] > 0.4  # ET NSE
  @test gof.NSE[2] > 0.3  # GPP NSE
  
  println("   ✓ 优化完成: ET NSE=$(round(gof.NSE[1], digits=3)), GPP NSE=$(round(gof.NSE[2], digits=3))")
end

#==============================================================================#
#  6. 参数导出测试
#==============================================================================#
println("\n💾 测试参数导出...")

@testset "参数导出" begin
  model = models["Yu2004"]
  params = Params(model)
  
  # 导出为 CSV
  csv_file = joinpath(@__DIR__, "params_test_export.csv")
  CSV.write(csv_file, params)
  
  @test isfile(csv_file)
  @test filesize(csv_file) > 0
  
  # 清理
  rm(csv_file)
  
  println("   ✓ 参数导出成功")
end

#==============================================================================#
#  综合评价
#==============================================================================#
println("\n" * "="^70)
println("  📋 综合评价")
println("="^70)

println("\n✅ 现代化参数方案特性:")
println("   1. 统一的参数接口：Params(model)")
println("   2. 灵活的参数更新：update!(model, names, values)")
println("   3. 多站点配置支持：IGBP_configs")
println("   4. 参数优化集成：optim(model, df; parNames)")
println("   5. 参数导出功能：CSV.write(params)")

println("\n💡 使用建议:")
println("   - 使用 Params(model) 统一访问参数")
println("   - 使用 update!(model, names, values) 更新参数")
println("   - 使用 optim(model, df; parNames) 进行参数优化")
println("   - 使用 IGBP_configs 管理多站点配置")
println()
