include("data_preprocessing.jl")
include("metrics_computation.jl")
include("pearson_relation_coeff.jl")
include("p_value_calculations.jl")
include("regression_plot_and_bland_altman_plot.jl")

# ===========================================================================
# MAIN 
# ===========================================================================

base_dir = raw"C:\Users\Heba\Desktop\PhD_Courses\AS\project\dataset\dataset\modified_signals"
all_hr_rates, all_br_rates, all_hr_ref = process_patients_data(base_dir)

# ===========================================================================
# SUMMARY — BCG vs Reference
# ===========================================================================
display_summary_for_data_analysis(all_hr_rates, all_hr_ref)


# ===========================================================================
# SUMMARY — metrics across all patients
# ===========================================================================

# Mean Absolute Error (MAE) && Root Mean Square Error (RMSE) && Mean Absolute Percentage Error (MAPE)
# ---------------------------------------------------------------------------
all_pids = sort(collect(keys(all_hr_rates)))
metrics  = compute_metrics(all_pids, all_hr_rates, all_hr_ref)

if !isnothing(metrics)
    println("\n" * "="^65)
    println("OVERALL METRICS ACROSS $(metrics.n) PATIENTS")
    println("="^65)
    println("  MAE  = $(round(metrics.mae,  digits=3)) bpm   ← avg absolute error per patient")
    println("  RMSE = $(round(metrics.rmse, digits=3)) bpm   ← penalises large errors more")
    println("  MAPE = $(round(metrics.mape, digits=3)) %     ← scale-independent percentage error")
end

# Pearson relation coeff
# ---------------------------------------------------------------------------
r = compute_pearson(all_pids, all_hr_rates, all_hr_ref)

# P-value of the correlation
# ---------------------------------------------------------------------------
n_pts = length(filter(pid -> haskey(all_hr_rates, pid) &&
                                           haskey(all_hr_ref,   pid), all_pids))
p, tstat, tail_label = compute_pvalue(r, n_pts, tail="right")

println("\n  Pearson r  = $(round(r,      digits=4))")
println("  t-stat     = $(round(tstat,   digits=4))  (df = $(n_pts - 2))")
println("  Test       = $tail_label")
println("  p-value    = $(round(p,       digits=6))")
println("  Significant (p < 0.05): $(p < 0.05 ? "YES ✓" : "NO ✗")")

# Regression Plot And Bland–Altman Plot
# ---------------------------------------------------------------------------
p1, p2 = plot_regression_and_bland_altman(
    all_pids,
    all_hr_rates,
    all_hr_ref;
    save_path = raw"C:\Users\Heba\Desktop\PhD_Courses\AS\\projectConvertedCode\results"
)

display(p1)   # show regression plot
display(p2)   # show bland-altman plot
