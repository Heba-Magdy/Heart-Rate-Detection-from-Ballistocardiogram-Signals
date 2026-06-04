"""
Complete Julia implementation of sleep/vitals analysis pipeline.
Multi-patient version — iterates over all patient BCG folders.
"""

using FFTW
using DSP
using Statistics
using LinearAlgebra
using CSV
using DataFrames
using Dates

# ===========================================================================
# PROCESS REFERENCE RR FILE
# ===========================================================================
function process_patients_data(base_dir::String)

    patient_folders = sort(filter(isdir, [joinpath(base_dir, d) for d in readdir(base_dir)]))

    all_hr_rates = Dict{String, Vector{Float64}}()
    all_br_rates = Dict{String, Vector{Float64}}()
    all_hr_ref   = Dict{String, Vector{Float64}}()

    for patient_dir in patient_folders
        patient_id = basename(patient_dir)
        bcg_dir    = joinpath(patient_dir, "BCG")

        if !isdir(bcg_dir)
            println("⚠ No BCG folder for patient $patient_id, skipping.")
            continue
        end

        bcg_files = sort(filter(f -> endswith(lowercase(f), ".csv"), readdir(bcg_dir, join=true)))
        if isempty(bcg_files)
            println("⚠ No CSV files in BCG folder for patient $patient_id, skipping.")
            continue
        end

        bcg_file           = bcg_files[1]
        hr_rates, br_rates = process_patient(bcg_file)

        if !isnothing(hr_rates)
            all_hr_rates[patient_id] = hr_rates
            all_br_rates[patient_id] = br_rates
        end

        rr_file = find_rr_file(patient_dir, bcg_file)
        if !isnothing(rr_file)
            all_hr_ref[patient_id] = process_patient_reference_data(rr_file, 10)
        end

    end

    return (all_hr_rates, all_br_rates, all_hr_ref)
end


function display_summary_for_data_analysis(all_hr_rates::Dict{String, Vector{Float64}},
                         all_hr_ref::Dict{String, Vector{Float64}})
    println("\n" * "="^70)
    println("COMPARISON — BCG Heart Rate vs Reference")
    println("="^70)
    println("Patient | BCG avg | Ref avg | BCG min | BCG max | Ref min | Ref max")
    println("-"^70)
    
    for pid in sort(collect(keys(all_hr_rates)))
        hr_bcg = all_hr_rates[pid]
        hr_ref = get(all_hr_ref, pid, Float64[])
    
        bcg_str = isempty(hr_bcg) ? "  N/A  " : "$(round(mean(hr_bcg), digits=1))"
        ref_str = isempty(hr_ref) ? "  N/A  " : "$(round(mean(hr_ref), digits=1))"
        bcg_min = isempty(hr_bcg) ? "N/A" : "$(round(minimum(hr_bcg)))"
        bcg_max = isempty(hr_bcg) ? "N/A" : "$(round(maximum(hr_bcg)))"
        ref_min = isempty(hr_ref) ? "N/A" : "$(round(minimum(hr_ref)))"
        ref_max = isempty(hr_ref) ? "N/A" : "$(round(maximum(hr_ref)))"
    
        println("  $pid    |  $bcg_str  |  $ref_str  |  $bcg_min  |  $bcg_max  |  $ref_min  |  $ref_max")
    end
    
    println("\nDone. Processed $(length(all_hr_rates)) patients.") 
end


function process_patient_reference_data(rr_file::String, window_sec::Int=10)
    println("  Reference file: $rr_file")

    # --- Parse RR file ---
    # Format: Timestamp,Heart Rate,RR Interval in seconds
    df = CSV.read(rr_file, DataFrame, header=1)

    # Parse timestamps
    fmt         = dateformat"yyyy/mm/d HH:MM:SS"
    timestamps  = DateTime.(String.(df[:, 1]), fmt)
    heart_rates = Float64.(df[:, 2])

    t_start = timestamps[1]
    t_end   = timestamps[end]
    total_sec = round(Int, (t_end - t_start).value / 1000)  # ms → seconds

    println("  RR duration: $(total_sec) seconds  Rows: $(nrow(df))")

    # --- Compute per-window average HR ---
    # For each 10-second window, collect all HR readings that fall in that window
    # then average them
    window_hr = Float64[]

    win_start = t_start
    win_end   = win_start + Second(window_sec)

    while win_start < t_end
        # Find all rows whose timestamp falls in [win_start, win_end)
        mask = (timestamps .>= win_start) .& (timestamps .< win_end)
        if any(mask)
            push!(window_hr, mean(heart_rates[mask]))
        else
            # No readings in this window — interpolate using nearest value
            diffs = abs.((timestamps .- win_start) ./ Millisecond(1))
            push!(window_hr, heart_rates[argmin(diffs)])
        end
        win_start = win_end
        win_end   = win_start + Second(window_sec)
    end

    println("  Reference windows: $(length(window_hr))")
    return window_hr
end

# ===========================================================================
# HELPER — find matching RR file for a BCG file
# ===========================================================================

function find_rr_file(patient_dir::String, bcg_file::String)
    rr_dir = joinpath(patient_dir, "RR")
    if !isdir(rr_dir)
        println("  ⚠ No RR folder found at: $rr_dir")
        return nothing
    end

    rr_files = sort(filter(f -> endswith(lowercase(f), ".csv"), readdir(rr_dir, join=true)))
    if isempty(rr_files)
        println("  ⚠ No RR CSV files found in: $rr_dir")
        return nothing
    end

    # Match by date: extract date string from BCG filename
    # BCG filename format: 01_20231104_BCG.csv → date = "20231104"
    bcg_name  = basename(bcg_file)
    date_part = match(r"\d{8}", bcg_name)   # find 8-digit date in filename

    if !isnothing(date_part)
        matched = filter(f -> occursin(date_part.match, basename(f)), rr_files)
        if !isempty(matched)
            return matched[1]
        end
    end

    # Fallback: just return first RR file
    println("  ⚠ No date-matched RR file found, using first available.")
    return rr_files[1]
end

# ===========================================================================
# WAVELET FILTERS (bior3.9)
# ===========================================================================

function get_bior39_filters()
    rec_lo = [0.0, 0.0, 0.0, 0.007782754936752098, -0.006935194945129213,
              -0.05765991509500517, 0.09039613683044469, 0.09787642690172737,
              -0.4266740994852978,  0.8563591981812583, -0.4266740994852978,
               0.09787642690172737,  0.09039613683044469, -0.05765991509500517,
              -0.006935194945129213, 0.007782754936752098, 0.0, 0.0, 0.0]

    dec_lo = [0.0, 0.0, 0.0, -0.007782754936752098, -0.006935194945129213,
               0.05765991509500517,  0.09039613683044469, -0.09787642690172737,
              -0.4266740994852978,  -0.8563591981812583, -0.4266740994852978,
              -0.09787642690172737,   0.09039613683044469,  0.05765991509500517,
              -0.006935194945129213, -0.007782754936752098, 0.0, 0.0, 0.0]

    L      = length(rec_lo)
    rec_hi = [(-1)^(k+1) * dec_lo[L - k] for k in 0:L-1]

    return rec_lo ./ sqrt(2), rec_hi ./ sqrt(2)
end

# ===========================================================================
# MODWT
# ===========================================================================

function modwt(x::Vector{Float64}, wname::String, J::Int)
    x = vec(x)
    datalength = length(x)

    Jmax = floor(Int, log2(datalength))
    if J <= 0 || J > Jmax
        error("Wavelet:modwt:MRALevel - J must be between 1 and floor(log2(N))")
    end

    siglen     = length(x)
    Lo, Hi     = get_bior39_filters()
    filter_len = length(Lo)
    Nrep       = siglen

    if siglen < filter_len
        reps = div(filter_len - siglen, siglen) + 1
        xp   = repeat(x, reps + 1)
        x    = xp[1 : siglen + (filter_len - siglen)]
        Nrep = length(x)
    end

    N    = Nrep
    G    = fft([Lo; zeros(N - filter_len)])
    H    = fft([Hi; zeros(N - filter_len)])
    Vhat = fft(x)

    function modwtdec(X, G, H, J)
        N        = length(X)
        upfactor = 2^J
        indices  = mod.(upfactor .* (0:N-1), N) .+ 1
        return G[indices] .* X, H[indices] .* X
    end

    w = Vector{Vector{Float64}}()
    for jj in 0:J-1
        Vhat, What = modwtdec(Vhat, G, H, jj)
        push!(w, real.(ifft(What)))
    end
    push!(w, real.(ifft(Vhat)))

    W = hcat(w...)'
    return W[:, 1:siglen]
end

# ===========================================================================
# MODWTMRA
# ===========================================================================

function modwtmra(w::Matrix{Float64}, wname::String)
    cfslength  = size(w, 2)
    J0         = size(w, 1) - 1
    nullinput  = zeros(Float64, cfslength)
    N          = cfslength
    Lo, Hi     = get_bior39_filters()
    filter_len = length(Lo)

    if cfslength < filter_len
        reps = div(filter_len - cfslength, cfslength) + 1
        wp   = repeat(w, 1, reps)
        w    = hcat(w, wp[:, 1:(filter_len - cfslength)])
        cfslength = size(w, 2)
        nullinput = zeros(Float64, cfslength)
    end

    G = fft([Lo; zeros(cfslength - filter_len)])
    H = fft([Hi; zeros(cfslength - filter_len)])

    function imodwtrec(Vin, Win, G, H, J)
        N_local  = length(Vin)
        upfactor = 2^J
        indices  = mod.(upfactor .* (0:N_local-1), N_local) .+ 1
        Gup = conj.(G[indices]); Hup = conj.(H[indices])
        return real.(ifft(Gup .* fft(Vin) .+ Hup .* fft(Win)))
    end

    function imodwtDetails(coefs, nullinput, lev, G, H, N)
        v = nullinput; w_inner = coefs
        for jj in (lev+1):-1:1
            vout    = imodwtrec(v, w_inner, G, H, jj-1)
            w_inner = copy(nullinput); v = vout
        end
        return v[1:N]
    end

    function imodwtSmooth(sc, nullinput, G, H, N, J0)
        v = sc
        for J in (J0+1):-1:1
            v = imodwtrec(v, copy(nullinput), G, H, J-1)
        end
        return v[1:N]
    end

    mra = Vector{Vector{Float64}}()
    for J in J0:-1:1
        wcfs = vec(w[J, :])
        push!(mra, imodwtDetails(wcfs, nullinput, J-1, G, H, cfslength)[1:N])
    end

    smooth = imodwtSmooth(vec(w[J0+1, :]), copy(nullinput), G, H, cfslength, J0-1)
    reverse!(mra)
    push!(mra, smooth[1:N])

    return hcat(mra...)'
end

# ===========================================================================
# DETECT PATTERNS
# ===========================================================================

function detect_patterns(pt1::Int, pt2::Int, win_size::Int,
                         data::Vector{<:Real}, time::Vector{<:Real};
                         plot_flag::Int=0)
    data = Float64.(data)
    time = Float64.(time)
    pt1_ = pt1; pt2_ = pt2

    limit       = floor(Int, length(data) / win_size)
    flag        = zeros(Int, length(data))
    event_flags = zeros(Int, limit)
    segments_sd = Float64[]

    for i in 1:limit
        push!(segments_sd, std(data[pt1+1:pt2+1], corrected=true))
        pt1 = pt2; pt2 += win_size
    end

    mad     = sum(abs.(segments_sd .- mean(segments_sd))) / length(segments_sd)
    thresh1 = sort(segments_sd)[max(1, floor(Int, 0.10 * length(segments_sd)))]

    thresh2 = 2.0 * mad
    pt1 = pt1_; pt2 = pt2_

    for j in 1:limit
        std_fos = round(segments_sd[j])
        if std_fos < thresh1
            flag[pt1+1:pt2+1] .= 3; event_flags[j] = 3
        elseif std_fos > thresh2
            flag[pt1+1:pt2+1] .= 2; event_flags[j] = 2
        else
            flag[pt1+1:pt2+1] .= 1; event_flags[j] = 1
        end
        pt1 = pt2; pt2 += win_size
    end

    mask = flag .== 1
    return data[mask], time[mask]
end

# ===========================================================================
# BAND PASS FILTERING
# ===========================================================================

function scipy_filtfilt(filter_coefs, x::Vector{Float64})
    b  = coefb(filter_coefs)
    a  = coefa(filter_coefs)
    n  = max(length(a), length(b))
    b2 = vcat(b, zeros(n - length(b)))
    a2 = vcat(a, zeros(n - length(a)))

    pad_len   = n * 3
    left_pad  = 2 * x[1]   .- x[pad_len+1:-1:2]
    right_pad = 2 * x[end] .- x[end-1:-1:end-pad_len]
    x_padded  = [left_pad; x; right_pad]

    comp = zeros(n-1, n-1)
    for i in 1:n-2; comp[i+1, i] = 1.0; end
    comp[1, :] = -a2[2:n] ./ a2[1]
    zi = (I - comp') \ (b2[2:n] .- a2[2:n] .* b2[1])

    y_fwd = filt(b, a, x_padded, zi .* x_padded[1])
    y_bwd = reverse(filt(b, a, reverse(y_fwd), zi .* y_fwd[end]))

    return y_bwd[pad_len+1 : pad_len+length(x)]
end

function band_pass_filtering(data::Vector{<:Real}, fs::Int, filter_type::String)
    data = Float64.(data)
    if filter_type == "bcg"
        bcg_          = scipy_filtfilt(digitalfilter(Highpass(0.8 / (fs/2)), Chebyshev1(2, 0.5)), data)
		filtered_data = scipy_filtfilt(digitalfilter(Lowpass(15.0 / (fs/2)), Chebyshev1(4, 0.5)), data)
    elseif filter_type == "breath"
        bcg_          = scipy_filtfilt(digitalfilter(Highpass(0.01 / (fs/2)), Chebyshev1(2, 0.5)), data)
        filtered_data = scipy_filtfilt(digitalfilter(Lowpass(0.4   / (fs/2)), Chebyshev1(4, 0.5)), bcg_)
    else
        filtered_data = data
    end
    return filtered_data
end

# ===========================================================================
# DETECT PEAKS
# ===========================================================================

function detect_peaks(x::Vector{Float64};
                      mph::Union{Nothing,Float64}=nothing,
                      mpd::Int=1, threshold::Float64=0.0,
                      edge::Union{Nothing,String}="rising",
                      kpsh::Bool=false, valley::Bool=false)
    x = copy(x)
    length(x) < 3 && return Int[]
    valley && (x = -x)

    dx           = x[2:end] .- x[1:end-1]
    indnan       = findall(isnan, x)
    if !isempty(indnan)
        x[indnan] .= Inf
        dx[findall(isnan, dx)] .= Inf
    end

    dx_r = [dx; 0.0]; dx_l = [0.0; dx]
    ine = Int[]; ire = Int[]; ife = Int[]
    if isnothing(edge)
        ine = findall((dx_r .< 0) .& (dx_l .> 0))
    else
        lowercase(edge) in ["rising",  "both"] && (ire = findall((dx_r .<= 0) .& (dx_l .> 0)))
        lowercase(edge) in ["falling", "both"] && (ife = findall((dx_r .<  0) .& (dx_l .>= 0)))
    end
    ind = sort(unique([ine; ire; ife]))

    if !isempty(ind) && !isempty(indnan)
        nb = filter(i -> 1 <= i <= length(x), unique([indnan; indnan.-1; indnan.+1]))
        ind = filter(i -> i ∉ nb, ind)
    end
    !isempty(ind) && ind[1]   == 1          && (ind = ind[2:end])
    !isempty(ind) && ind[end] == length(x)  && (ind = ind[1:end-1])
    !isempty(ind) && !isnothing(mph)        && (ind = filter(i -> x[i] >= mph, ind))
    if !isempty(ind) && threshold > 0
        ind = ind[min.(x[ind].-x[ind.-1], x[ind].-x[ind.+1]) .>= threshold]
    end
    if !isempty(ind) && mpd > 1
        ind  = ind[sortperm(x[ind], rev=true)]
        idel = falses(length(ind))
        for i in 1:length(ind)
            if !idel[i]
                tc = (ind .>= ind[i]-mpd) .& (ind .<= ind[i]+mpd)
                kpsh && (tc .&= x[ind[i]] .> x[ind])
                idel .|= tc; idel[i] = false
            end
        end
        ind = sort(ind[.!idel])
    end
    !isempty(indnan) && (x[indnan] .= NaN)
    return ind .- 1   # 0-based to match Python
end

# ===========================================================================
# COMPUTE RATE  (fs-based, no utc_time needed)
# ===========================================================================

function compute_rate(beats::Vector{Float64}, fs::Int, mpd::Int)
    indices = detect_peaks(beats, mpd=mpd)
    println("  peaks: $(length(indices))  spacings(samples): $(diff(indices)[1:min(5,end)])")
    if length(indices) > 1
        mean_samples = mean(diff(indices))
        mean_ms      = mean_samples * (1000.0 / fs)
        bpm_avg      = 1000.0 * (60.0 / mean_ms)
        return round(bpm_avg, digits=2), indices
    else
        return 0.0, Int[]
    end
end

# ===========================================================================
# VITALS
# ===========================================================================

function vitals(t1::Int, t2::Int, win_size::Int, window_limit::Int,
                sig::Vector{Float64}, fs::Int, mpd::Int)
    all_rate = Float64[]
    for j in 0:window_limit-1
        rate, _ = compute_rate(sig[t1+1:t2], fs, mpd)
        push!(all_rate, rate)
        t1 = t2; t2 += win_size
    end
    return all_rate
end

# ===========================================================================
# PROCESS ONE PATIENT FILE
# ===========================================================================

function process_patient(file::String)
    println("\n>>> Processing: $file")

    # --- Parse new BCG file format ---
    lines     = readlines(file)
    row2      = split(strip(lines[2]), ",")
    fs        = parse(Int,     row2[3])
    first_val = parse(Float64, row2[1])
    remaining = [parse(Float64, strip(l)) for l in lines[3:end] if !isempty(strip(l))]
    data_stream = vcat(first_val, remaining)

    # Dummy time vector (sample indices) — no UTC needed
    time_vec = collect(0:length(data_stream)-1)

    println("  Samples: $(length(data_stream))  fs: $(fs) Hz  Duration: $(round(length(data_stream)/fs/3600, digits=2))h")

    # Window size = 10 seconds
    win_size     = fs * 10
    start_point  = 0
    end_point    = win_size

    # Step 1: Detect patterns
    data_stream, time_vec = detect_patterns(
        start_point, end_point, win_size,
        data_stream, time_vec; plot_flag=0)
    println("  After detect_patterns: $(length(data_stream)) samples")

    length(data_stream) < win_size && (println("  ⚠ Too few samples after filtering, skipping."); return nothing, nothing)

    # Step 2: Bandpass filters
    movement  = band_pass_filtering(data_stream, fs, "bcg")
    breathing = band_pass_filtering(data_stream, fs, "breath")

    # Step 3: Wavelet transform
    wavelet_cycle = movement

    # Step 4: Vital signs windows
    t1            = 0
    t2            = win_size
    window_shift2 = win_size
    limit         = floor(Int, length(breathing) / window_shift2)

    # mpd scaled to fs
    # Heart rate: max ~200bpm → min spacing = fs*60/200 samples
    # Breathing:  max ~60bpm  → min spacing = fs*60/60  samples
    mpd_hr = floor(Int, fs * 60 / 120) 
    mpd_br = max(1, floor(Int, fs * 60 / 60))

    hr_beats = vitals(t1, t2, window_shift2, limit, wavelet_cycle, fs, mpd_hr)
    br_beats = vitals(t1, t2, window_shift2, limit, breathing,     fs, mpd_br)

    hr_nz = filter(x -> x > 0, hr_beats)
    br_nz = filter(x -> x > 0, br_beats)

    if !isempty(hr_nz)
        println("  Heart Rate     — min: $(round(minimum(hr_nz)))  max: $(round(maximum(hr_nz)))  avg: $(round(mean(hr_nz)))")
    else
        println("  Heart Rate     — no valid windows")
    end
    if !isempty(br_nz)
        println("  Breathing Rate — min: $(round(minimum(br_nz)))  max: $(round(maximum(br_nz)))  avg: $(round(mean(br_nz)))")
    else
        println("  Breathing Rate — no valid windows")
    end

    return hr_nz, br_nz
end
