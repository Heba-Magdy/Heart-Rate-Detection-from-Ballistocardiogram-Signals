# ===========================================================================
# PEARSON CORRELATION COEFFICIENT (manual implementation)
# ===========================================================================

function compute_pearson(patient_ids::Vector{String},
                         bcg_rates::Dict{String, Vector{Float64}},
                         ref_rates::Dict{String, Vector{Float64}})

    predicted = Float64[]
    reference = Float64[]

    for pid in patient_ids
        has_bcg = haskey(bcg_rates, pid) && !isempty(bcg_rates[pid])
        has_ref = haskey(ref_rates, pid)  && !isempty(ref_rates[pid])
        if has_bcg && has_ref
            push!(predicted, mean(bcg_rates[pid]))
            push!(reference, mean(ref_rates[pid]))
        end
    end

    n = length(predicted)
    if n < 2
        println("⚠ Need at least 2 patients to compute correlation.")
        return NaN
    end

    # Step 1: compute means
    mean_pred = 0.0
    mean_ref  = 0.0
    for i in 1:n
        mean_pred += predicted[i]
        mean_ref  += reference[i]
    end
    mean_pred /= n
    mean_ref  /= n

    # Step 2: compute deviations from mean
    # numerator   = Σ (pred_i - mean_pred) * (ref_i - mean_ref)
    # denominator = sqrt( Σ(pred_i - mean_pred)² * Σ(ref_i - mean_ref)² )
    numerator    = 0.0
    sum_sq_pred  = 0.0
    sum_sq_ref   = 0.0
    for i in 1:n
        dev_pred  = predicted[i] - mean_pred
        dev_ref   = reference[i] - mean_ref
        numerator   += dev_pred * dev_ref
        sum_sq_pred += dev_pred^2
        sum_sq_ref  += dev_ref^2
    end

    denominator = sqrt(sum_sq_pred * sum_sq_ref)

    # Step 3: guard against zero denominator (all values identical)
    if denominator == 0.0
        println("⚠ Denominator is zero — all predicted or reference values are identical.")
        return NaN
    end

    r = numerator / denominator

    # Step 4: clamp to [-1, 1] to handle floating point drift
    r = max(-1.0, min(1.0, r))
	
	println("\n  Pearson r = $(round(r, digits=4))")
	println("  Interpretation:")
	if abs(r) >= 0.9
		println("    → Very strong correlation")
	elseif abs(r) >= 0.7
		println("    → Strong correlation")
	elseif abs(r) >= 0.5
		println("    → Moderate correlation")
	elseif abs(r) >= 0.3
		println("    → Weak correlation")
	else
		println("    → Little to no correlation")
	end
	
    return r
end