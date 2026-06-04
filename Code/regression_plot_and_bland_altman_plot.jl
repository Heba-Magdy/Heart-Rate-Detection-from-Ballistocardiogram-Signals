# ===========================================================================
# REGRESSION PLOT & BLAND-ALTMAN PLOT
# ===========================================================================

using Plots

function plot_regression_and_bland_altman(patient_ids::Vector{String},
                                           bcg_rates::Dict{String, Vector{Float64}},
                                           ref_rates::Dict{String, Vector{Float64}};
                                           save_path::String=".")

    # --- Build per-patient scalar vectors ---
    predicted = Float64[]
    reference = Float64[]
    valid_ids = String[]

    for pid in patient_ids
        has_bcg = haskey(bcg_rates, pid) && !isempty(bcg_rates[pid])
        has_ref = haskey(ref_rates,  pid) && !isempty(ref_rates[pid])
        if has_bcg && has_ref
            push!(predicted, mean(bcg_rates[pid]))
            push!(reference, mean(ref_rates[pid]))
            push!(valid_ids, pid)
        end
    end

    n = length(predicted)
    if n < 2
        println("⚠ Not enough patients for plots.")
        return
    end

    # ===========================================================
    # 1. REGRESSION PLOT
    # ===========================================================
    # Manually compute linear regression coefficients
    # y = slope * x + intercept
    # slope     = Σ(xi - x̄)(yi - ȳ) / Σ(xi - x̄)²
    # intercept = ȳ - slope * x̄

    mean_ref  = sum(reference) / n
    mean_pred = sum(predicted) / n

    num   = 0.0
    denom = 0.0
    for i in 1:n
        num   += (reference[i] - mean_ref) * (predicted[i] - mean_pred)
        denom += (reference[i] - mean_ref)^2
    end
    slope     = num / denom
    intercept = mean_pred - slope * mean_ref

    # Regression line points
    x_min  = minimum(reference) - 5
    x_max  = maximum(reference) + 5
    x_line = [x_min, x_max]
    y_line = slope .* x_line .+ intercept

    # Identity line (perfect agreement): y = x
    id_line = [x_min, x_max]

    p1 = scatter(reference, predicted,
                 xlabel     = "Reference HR (bpm)",
                 ylabel     = "BCG Estimated HR (bpm)",
                 title      = "Regression Plot\n(slope=$(round(slope,digits=3)), intercept=$(round(intercept,digits=3)))",
                 label      = "Patients",
                 color      = :steelblue,
                 markersize = 7,
                 legend     = :topleft)

    plot!(p1, x_line, y_line,
          label     = "Regression line",
          color     = :red,
          linewidth = 2)

    plot!(p1, id_line, id_line,
          label     = "Identity (y=x)",
          color     = :black,
          linewidth = 1,
          linestyle = :dash)

    # Annotate each point with patient ID
    for i in 1:n
        annotate!(p1, reference[i], predicted[i],
                  text("  $(valid_ids[i])", 7, :left, :gray))
    end

    # ===========================================================
    # 2. BLAND-ALTMAN PLOT
    # ===========================================================
    # x-axis: mean of BCG and reference   = (predicted + reference) / 2
    # y-axis: difference                  = predicted - reference
    # Bias   = mean(differences)
    # Limits of agreement = bias ± 1.96 * std(differences)

    means = [(predicted[i] + reference[i]) / 2.0 for i in 1:n]
    diffs = [predicted[i] - reference[i]          for i in 1:n]

    # Bias (mean difference) — manual
    bias = sum(diffs) / n

    # Std of differences — manual
    sq_sum = sum((d - bias)^2 for d in diffs)
    std_diff = sqrt(sq_sum / (n - 1))

    loa_upper = bias + 1.96 * std_diff
    loa_lower = bias - 1.96 * std_diff

    x_ba_min = minimum(means) - 5
    x_ba_max = maximum(means) + 5

    p2 = scatter(means, diffs,
                 xlabel     = "Mean of BCG & Reference HR (bpm)",
                 ylabel     = "Difference BCG - Reference (bpm)",
                 title      = "Bland-Altman Plot",
                 label      = "Patients",
                 color      = :steelblue,
                 markersize = 7,
                 legend     = :topright)

    # Bias line
    hline!(p2, [bias],
           label     = "Bias = $(round(bias, digits=2)) bpm",
           color     = :red,
           linewidth = 2)

    # Upper limit of agreement
    hline!(p2, [loa_upper],
           label     = "+1.96 SD = $(round(loa_upper, digits=2)) bpm",
           color     = :orange,
           linewidth = 1,
           linestyle = :dash)

    # Lower limit of agreement
    hline!(p2, [loa_lower],
           label     = "-1.96 SD = $(round(loa_lower, digits=2)) bpm",
           color     = :orange,
           linewidth = 1,
           linestyle = :dash)

    # Zero reference line
    hline!(p2, [0.0],
           label     = "Zero",
           color     = :black,
           linewidth = 1,
           linestyle = :dot)

    # Annotate each point with patient ID
    for i in 1:n
        annotate!(p2, means[i], diffs[i],
                  text("  $(valid_ids[i])", 7, :left, :gray))
    end

    # ===========================================================
    # 3. PRINT BLAND-ALTMAN STATISTICS
    # ===========================================================
    println("\n" * "="^50)
    println("BLAND-ALTMAN STATISTICS")
    println("="^50)
    println("  Bias (mean diff)     = $(round(bias,      digits=3)) bpm")
    println("  Std of differences   = $(round(std_diff,  digits=3)) bpm")
    println("  Upper LoA (+1.96 SD) = $(round(loa_upper, digits=3)) bpm")
    println("  Lower LoA (-1.96 SD) = $(round(loa_lower, digits=3)) bpm")
    println("\n  Interpretation:")
    println("  Bias close to 0     → $(abs(bias) < 5 ? "✓ Good" : "✗ Systematic offset exists")")
    println("  Narrow LoA band     → $(abs(loa_upper - loa_lower) < 20 ? "✓ Good agreement" : "✗ Wide limits — poor agreement")")

    # ===========================================================
    # 4. COMBINE AND SAVE
    # ===========================================================
    combined = plot(p1, p2,
                    layout  = (1, 2),
                    size    = (1200, 500),
                    margin  = 5Plots.mm)

    savefig(combined, joinpath(save_path, "regression_bland_altman.png"))
    println("\n  Plots saved to: $(joinpath(save_path, "regression_bland_altman.png"))")

    return p1, p2
end