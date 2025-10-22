# VCmax垂向分布示例
#
# 本示例展示如何使用VCmax垂向分布函数来计算冠层的最大羧化速率分布
#
# 作者: Claude
# 日期: 2025

using SPAC

# ============================================================================
# 示例1: 基本的VCmax垂向分布计算
# ============================================================================

println("=".^80)
println("示例1: 计算10层冠层的VCmax垂向分布")
println("=".^80)

# 设置参数（基于Chen et al., 2012, GBC）
nlyr = 10           # 冠层层数
lai = 5.0           # 总叶面积指数 [m² m⁻²]
Ω = 1.0             # 遮挡指数（clumping index）[-]
CosZs = 0.866       # cos(30°) - 太阳天顶角30度
VCmax25 = 62.5      # 25°C下的最大羧化速率 [μmol m⁻² s⁻¹]
N_leaf = 3.10       # 冠层平均叶片氮含量 [g m⁻²]
slope = 20.72 / 62.5 # Vcmax-N曲线的斜率（归一化）

# 计算垂向分布
VCmax_sunlit, VCmax_shaded, LAI_sunlit, LAI_shaded =
    VCmax_profile(nlyr, lai, Ω, CosZs, VCmax25, N_leaf, slope)

# 输出结果
println("\n各层VCmax分布（从冠层顶部到底部）:")
println("-".^80)
println("层号 | 累积LAI | 向阳LAI | 背阴LAI | VCmax_sunlit | VCmax_shaded")
println("-".^80)

for i in 1:nlyr
    L_mid = (i - 0.5) * lai / nlyr  # 层中心的累积LAI
    @printf("%3d  | %7.2f | %7.3f | %7.3f | %12.2f | %12.2f\n",
            i, L_mid, LAI_sunlit[i], LAI_shaded[i],
            VCmax_sunlit[i], VCmax_shaded[i])
end

println("\n总LAI验证:")
println("  输入总LAI: $lai")
println("  计算总LAI: $(round(sum(LAI_sunlit) + sum(LAI_shaded), digits=4))")
println("  向阳总LAI: $(round(sum(LAI_sunlit), digits=3))")
println("  背阴总LAI: $(round(sum(LAI_shaded), digits=3))")

# ============================================================================
# 示例2: 计算冠层平均VCmax（用于大叶模型）
# ============================================================================

println("\n" * "=".^80)
println("示例2: 计算考虑垂向分布的冠层平均VCmax")
println("=".^80)

# 使用相同参数
VCmax_mean, VCmax_sunlit_mean, VCmax_shaded_mean =
    VCmax_profile_mean(nlyr, lai, Ω, CosZs, VCmax25, N_leaf, slope)

println("\n冠层平均VCmax（按LAI加权）:")
println("  总体平均VCmax:    $(round(VCmax_mean, digits=2)) μmol m⁻² s⁻¹")
println("  向阳叶平均VCmax:  $(round(VCmax_sunlit_mean, digits=2)) μmol m⁻² s⁻¹")
println("  背阴叶平均VCmax:  $(round(VCmax_shaded_mean, digits=2)) μmol m⁻² s⁻¹")

# 与原始方法比较
VCmax_sunlit_old, VCmax_shaded_old =
    VCmax(lai, Ω, CosZs, VCmax25, N_leaf, slope)

println("\n与原始二分法（向阳/背阴）比较:")
println("  原始方法 - 向阳叶VCmax: $(round(VCmax_sunlit_old, digits=2)) μmol m⁻² s⁻¹")
println("  原始方法 - 背阴叶VCmax: $(round(VCmax_shaded_old, digits=2)) μmol m⁻² s⁻¹")

# ============================================================================
# 示例3: 不同太阳天顶角的影响
# ============================================================================

println("\n" * "=".^80)
println("示例3: 太阳天顶角对VCmax分布的影响")
println("=".^80)

zenith_angles = [0, 30, 45, 60, 75]  # 天顶角（度）
println("\n天顶角 | CosZs | VCmax_mean | VCmax_sunlit_mean | VCmax_shaded_mean")
println("-".^80)

for zenith in zenith_angles
    cos_zenith = cosd(zenith)
    if cos_zenith > 0
        vcmax_m, vcmax_sun_m, vcmax_sha_m =
            VCmax_profile_mean(nlyr, lai, Ω, cos_zenith, VCmax25, N_leaf, slope)
        @printf("%4d°   | %.3f | %10.2f | %17.2f | %17.2f\n",
                zenith, cos_zenith, vcmax_m, vcmax_sun_m, vcmax_sha_m)
    end
end

# ============================================================================
# 示例4: 不同LAI的影响
# ============================================================================

println("\n" * "=".^80)
println("示例4: 叶面积指数对VCmax分布的影响")
println("=".^80)

LAI_values = [1.0, 2.0, 3.0, 5.0, 7.0]
println("\nLAI   | VCmax_mean | 向阳LAI | 背阴LAI")
println("-".^80)

for lai_test in LAI_values
    vcmax_m, _, _ = VCmax_profile_mean(nlyr, lai_test, Ω, CosZs, VCmax25, N_leaf, slope)
    _, _, lai_sun, lai_sha = VCmax_profile(nlyr, lai_test, Ω, CosZs, VCmax25, N_leaf, slope)
    @printf("%.1f   | %10.2f | %7.3f | %7.3f\n",
            lai_test, vcmax_m, sum(lai_sun), sum(lai_sha))
end

# ============================================================================
# 示例5: 可视化VCmax垂向剖面
# ============================================================================

println("\n" * "=".^80)
println("示例5: VCmax垂向剖面（ASCII图）")
println("=".^80)

# 计算20层的高分辨率分布
nlyr_high = 20
VCmax_sunlit_hr, VCmax_shaded_hr, _, _ =
    VCmax_profile(nlyr_high, lai, Ω, CosZs, VCmax25, N_leaf, slope)

println("\n深度(LAI) |" * " ".^15 * "VCmax (μmol m⁻² s⁻¹)")
println("-".^10 * "+" * "-".^60)

max_vcmax = maximum([VCmax_sunlit_hr; VCmax_shaded_hr])
for i in 1:nlyr_high
    L_mid = (i - 0.5) * lai / nlyr_high
    # 归一化到60个字符
    n_chars_sun = Int(round(VCmax_sunlit_hr[i] / max_vcmax * 30))
    n_chars_sha = Int(round(VCmax_shaded_hr[i] / max_vcmax * 30))

    bar_sun = "█"^n_chars_sun
    bar_sha = "░"^n_chars_sha

    @printf("%6.2f    | %s%s (%.1f/%.1f)\n",
            L_mid, bar_sun, bar_sha, VCmax_sunlit_hr[i], VCmax_shaded_hr[i])
end

println("\n图例: █ = 向阳叶, ░ = 背阴叶")

# ============================================================================
# 总结
# ============================================================================

println("\n" * "=".^80)
println("总结")
println("=".^80)
println("""
VCmax垂向分布功能说明:

1. VCmax_profile(): 计算多层垂向分布
   - 输入: 层数, LAI, 遮挡指数, 太阳角度, VCmax25, 叶片氮含量, 斜率
   - 输出: 每层的向阳/背阴叶的VCmax和LAI

2. VCmax_profile_mean(): 计算冠层平均VCmax
   - 基于垂向分布计算LAI加权平均值
   - 适用于大叶模型

3. 应用场景:
   - 多层冠层模型的光合作用模拟
   - 大叶模型的精确化
   - 辐射传输模型的耦合
   - 碳通量的垂向分布研究

4. 理论基础:
   - 基于Chen et al., 2012, GBC的氮素垂向分布模型
   - 考虑光和氮素的消光效应
   - 区分向阳和背阴叶

参考文献:
Chen, J.M., et al. (2012). Daily canopy photosynthesis model through
temporal and spatial scaling for remote sensing applications.
Ecological Modelling, 124, 99-119.
""")

println("=".^80)
