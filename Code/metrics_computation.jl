# ===========================================================================
# COMPARISON METRICS — one value per patient, across all patients
# ===========================================================================

function compute_metrics(patient_ids::Vector{String},
                         bcg_rates::Dict{String, Vector{Float64}},
                         ref_rates::Dict{String, Vector{Float64}})

    predicted = Float64[]
    reference = Float64[]
    valid_ids = String[]

    # Build one scalar per patient
    for pid in patient_ids
        has_bcg = haskey(bcg_rates, pid) && !isempty(bcg_rates[pid])
        has_ref = haskey(ref_rates, pid)  && !isempty(ref_rates[pid])

        if has_bcg && has_ref
            push!(predicted, mean(bcg_rates[pid]))
            push!(reference, mean(ref_rates[pid]))
            push!(valid_ids, pid)
        else
            println("  ⚠ Patient $pid skipped — missing BCG or reference data")
        end
    end

    n = length(predicted)
    if n == 0
        println("⚠ No patients with both BCG and reference data.")
        return nothing
    end

    # Print per-patient breakdown
    println("\n" * "="^65)
    println("PER-PATIENT AVERAGES")
    println("="^65)
    println("Patient | BCG avg (bpm) | Ref avg (bpm) | Error (bpm)")
    println("-"^65)
    for i in 1:n
        err = predicted[i] - reference[i]
        println("  $(valid_ids[i])    |   $(round(predicted[i], digits=1))        |   $(round(reference[i], digits=1))        |   $(round(err, digits=1))")
    end

    # --- MAE ---
    # MAE = (1/n) * Σ |pred_i - ref_i|
    mae = 0.0
    for i in 1:n
        mae += abs(predicted[i] - reference[i])
    end
    mae /= n

    # --- RMSE ---
    # RMSE = sqrt( (1/n) * Σ (pred_i - ref_i)^2 )
    mse = 0.0
    for i in 1:n
        mse += (predicted[i] - reference[i])^2
    end
    rmse = sqrt(mse / n)

    # --- MAPE ---
    # MAPE = (100/n) * Σ |pred_i - ref_i| / |ref_i|
    mape       = 0.0
    mape_count = 0
    for i in 1:n
        if reference[i] != 0.0
            mape += abs(predicted[i] - reference[i]) / abs(reference[i])
            mape_count += 1
        end
    end
    mape = mape_count > 0 ? (mape / mape_count) * 100.0 : NaN

    return (n=n, mae=mae, rmse=rmse, mape=mape, ids=valid_ids,
            predicted=predicted, reference=reference)
end