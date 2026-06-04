
# ===========================================================================
# P-VALUE FOR PEARSON CORRELATION (manual implementation)
# ===========================================================================

function compute_pvalue(r::Float64, n::Int;
                        tail::String="right")   # "two", "right", "left"

    if n < 3
        println("⚠ Need at least 3 patients to compute p-value.")
        return NaN, NaN
    end

    if abs(r) == 1.0
        return 0.0, r * Inf
    end

    # Step 1: t-statistic
    # t = r * sqrt(n-2) / sqrt(1 - r²)
    t_stat = r * sqrt(n - 2) / sqrt(1.0 - r^2)

    df = Float64(n - 2)

    # Step 2: compute P(T <= t) using regularised Beta
    # P(T <= t) = 1 - 0.5 * I_x(df/2, 0.5)  where x = df/(df + t²)
    x        = df / (df + t_stat^2)
    cum_prob = 1.0 - 0.5 * regularised_beta(x, df/2.0, 0.5)

    # Step 3: compute p-value based on tail
    if tail == "two"
        # H₁: r ≠ 0  →  p = 2 * P(T > |t|)
        p_value = 2.0 * (1.0 - cum_prob)
        tail_label = "two-tailed (H₁: r ≠ 0)"

    elseif tail == "right"
        # H₁: r > 0  →  p = P(T > t)
        p_value = 1.0 - cum_prob
        tail_label = "right-tailed (H₁: r > 0)"

    elseif tail == "left"
        # H₁: r < 0  →  p = P(T < t)
        p_value = cum_prob
        tail_label = "left-tailed (H₁: r < 0)"

    else
        error("tail must be 'two', 'right', or 'left'")
    end

    p_value = max(0.0, min(1.0, p_value))

    return p_value, t_stat, tail_label
end

# ---------------------------------------------------------------------------
# Regularised incomplete beta function I_x(a, b)
# Uses continued fraction expansion (Lentz algorithm) for numerical stability
# ---------------------------------------------------------------------------

function regularised_beta(x::Float64, a::Float64, b::Float64)

    if x < 0.0 || x > 1.0
        error("x must be in [0, 1]")
    end
    x == 0.0 && return 0.0
    x == 1.0 && return 1.0

    # Use symmetry relation for better convergence:
    # I_x(a,b) = 1 - I_(1-x)(b,a)  when x > (a+1)/(a+b+2)
    if x > (a + 1.0) / (a + b + 2.0)
        return 1.0 - regularised_beta(1.0 - x, b, a)
    end

    # Log of the beta function normalisation factor
    # ln B(a,b) = ln Γ(a) + ln Γ(b) - ln Γ(a+b)
    log_beta_ab = log_gamma(a) + log_gamma(b) - log_gamma(a + b)

    # Front factor: x^a * (1-x)^b / B(a,b)
    front = exp(a * log(x) + b * log(1.0 - x) - log_beta_ab) / a

    # Continued fraction via Lentz algorithm
    cf = beta_continued_fraction(x, a, b)

    return front * cf
end

function beta_continued_fraction(x::Float64, a::Float64, b::Float64;
                                  max_iter::Int=200, tol::Float64=1e-12)
    # Lentz algorithm for the continued fraction representation
    qab = a + b
    qap = a + 1.0
    qam = a - 1.0

    c = 1.0
    d = 1.0 - qab * x / qap
    abs(d) < 1e-30 && (d = 1e-30)
    d = 1.0 / d
    h = d

    for m in 1:max_iter
        m2 = 2 * m

        # Even step
        aa =  m * (b - m) * x / ((qam + m2) * (a + m2))
        d  = 1.0 + aa * d;  abs(d) < 1e-30 && (d = 1e-30)
        c  = 1.0 + aa / c;  abs(c) < 1e-30 && (c = 1e-30)
        d  = 1.0 / d
        h *= d * c

        # Odd step
        aa = -(a + m) * (qab + m) * x / ((a + m2) * (qap + m2))
        d  = 1.0 + aa * d;  abs(d) < 1e-30 && (d = 1e-30)
        c  = 1.0 + aa / c;  abs(c) < 1e-30 && (c = 1e-30)
        d  = 1.0 / d
        delta = d * c
        h    *= delta

        abs(delta - 1.0) < tol && break
    end

    return h
end

function log_gamma(x::Float64)
    # Lanczos approximation for ln Γ(x), accurate to ~15 significant figures
    if x < 0.5
        # Reflection formula: Γ(x)Γ(1-x) = π/sin(πx)
        return log(π / sin(π * x)) - log_gamma(1.0 - x)
    end

    x -= 1.0
    coeffs = [0.99999999999980993, 676.5203681218851, -1259.1392167224028,
              771.32342877765313, -176.61502916214059, 12.507343278686905,
              -0.13857109526572012, 9.9843695780195716e-6, 1.5056327351493116e-7]
    g = 7.0
    t = x + g + 0.5
    s = coeffs[1]
    for i in 2:9
        s += coeffs[i] / (x + i - 1)
    end
    return 0.5 * log(2π) + (x + 0.5) * log(t) - t + log(s)
end


