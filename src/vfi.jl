struct Model
    wageV :: Vector{Float64}
    R :: Float64
    beta :: Float64
    u :: AbstractUtility
    # Number of capital grid points
    nk :: Int64
    # Max k grid point
    kMax :: Float64
end

## ----------  Helpers

c_min(m :: Model) = 1e-4;
betar(m :: Model) = m.beta * m.R;
n_periods(m :: Model) = length(m.wageV);
wage(m :: Model, t) = m.wageV[t];

# g + g^2 + ... + g^(T-1)
pvfactor(g, T) = (g ^ T - 1.0) / (g - 1.0);

present_value(xV, R) = sum(xV ./ (R .^ (0 : (length(xV)-1))));

euler_dev(m :: Model, ctV) = euler_dev(m.u, ctV, betar(m));

budget_kprime(w, R, k, c) = w .+ R .* k .- c;
budget_c(w, R, k, kPrime) = w .+ R .* k .- kPrime;

make_k_grid(m :: Model, t) = LinRange(0.0, m.kMax, m.nk);
kprime_max(m :: Model, t) = m.kMax;
# Max kPrime consistent with c > c_min
kprime_max(m :: Model, t, k) = budget_kprime(wage(m, t), m.R, k, c_min(m));
kprime_min(m :: Model, t) = 0.0;


function init_test_model()
    T = 5;
    wageV = collect(LinRange(1.0, 2.0, T));
    kMax = 20.0;
    u = UtilityCRRA(2.0);
    return Model(wageV, 1.04, 0.98, u, 20, kMax)
end


## ---------- Solve by Policy Fct Iteration

"""
	$(SIGNATURES)

Solve a model with policy function iteration.
Returns `Vector`s with interpolated policy functions k'(k) and c(k).
"""
function solve(m :: Model)
    T = n_periods(m);

    # Results are stored as interpolation objects
    kPrime_tV = Vector{Any}(undef, T);
    c_tV = Vector{Any}(undef, T);

    # Last period: easy
    kPrime_tV[T], c_tV[T] = solve_last_period(m);

    # Backward induction
    for t = (T-1) : -1 : 1
        kGridV, kPrime_tV[t], c_tV[t] = solve_one_period(m, t, kPrime_tV[t+1]);
    end

    isValid = check_solution(m, kPrime_tV, c_tV);
    if !isValid
        @warn "Solution fails checks"
    end

    return kPrime_tV, c_tV
end

function solve_last_period(m :: Model)
    T = n_periods(m);
    kGridV = make_k_grid(m, T);
    kPrimeFct = interpolate_kprime(kGridV, zeros(length(kGridV)));
    cV = budget_c(wage(m, T), m.R, kGridV, 0.0);
    cFct = interpolate_kprime(kGridV, cV);
    return kPrimeFct, cFct
end

function solve_one_period(m :: Model, t, kPrimeTomorrow)
    kGridV = make_k_grid(m, t);
    kPrimeV = zeros(length(kGridV));
    for (ik, k) in enumerate(kGridV)
        kPrimeV[ik] = solve_one_point(m, t, k, kPrimeTomorrow);
    end
    cV = budget_c(wage(m, t), m.R, kGridV, kPrimeV);
    kPrimeFct = interpolate_kprime(kGridV, kPrimeV);
    cFct = interpolate_kprime(kGridV, cV);
    return kGridV, kPrimeFct, cFct
end

function solve_one_point(m :: Model, t, k, kPrimeTomorrow)
    e_dev(kPrime) = euler_dev_one_point(m, t, k, kPrime, kPrimeTomorrow);
    kPrimeMin = 0.0;
    kPrimeMax = min(kprime_max(m, t, k), kprime_max(m, t));
    @assert kPrimeMax > kPrimeMin

    # Try corners
    eDev0 = e_dev(kPrimeMin);
    if eDev0 <= 0.0
        kPrime = kPrimeMin;
    else
        eDevMax = e_dev(kPrimeMax);
        if eDevMax >= 0.0
            kPrime = kPrimeMax
        else
            kPrime = find_zero(e_dev, (kPrimeMin, kPrimeMax, Bisection()));
        end
    end
    @assert kPrimeMin <= kPrime <= kprime_max(m, t)
    return kPrime
end

# Euler equation deviation.
# `dev > 0` means that `kPrime` implies a `c` today that is too high.
function euler_dev_one_point(m :: Model, t, k, kPrime, kPrimeTomorrow)
    cPrime = budget_c(wage(m, t+1), m.R, kPrime, kPrimeTomorrow(kPrime));
    betaRUprime = betar(m) * marg_utility(m.u, cPrime);
    dev = budget_c(wage(m, t), m.R, k, kPrime) - inv_marg_utility(m.u, betaRUprime);
    return dev
end

# Alternative Euler equation deviation.
# Easier to understand, but more non-linear.
function euler_dev_one_point2(m :: Model, t, k, kPrime, kPrimeTomorrow)
    c = budget_c(wage(m, t), m.R, k, kPrime);
    kPrimePrime = kPrimeTomorrow(kPrime);
    cPrime = budget_c(wage(m, t+1), m.R, kPrime, kPrimePrime);
    eDev = euler_dev(m, [c, cPrime]);
    return only(eDev)
end


function interpolate_kprime(kGridV :: AbstractVector{F}, kPrimeV :: AbstractVector{F}) where F
    interp = LinearInterpolation(kGridV, kPrimeV);
    return interp
end


function check_solution(m :: Model, kPrime_tV, c_tV)
    bcValid = true;
    eulerValid = true;

    T = n_periods(m);
    for t = 1 : (T-1)
        kGridV = make_k_grid(m, t);
        kp_fct = kPrime_tV[t];
        c_fct = c_tV[t];
        wage_t = wage(m, t);
        for (ik, k) in enumerate(kGridV)
            kPrime = kp_fct(k);
            c = budget_c(wage_t, m.R, k, kPrime);
            bcValid = bcValid && isapprox(c, c_fct(k), atol = 1e-4);

            if kprime_min(m, t) < kPrime < kprime_max(m, t)
                kPrimePrime = kPrime_tV[t+1](kPrime);
                cPrime = budget_c(wage(m, t+1), m.R, kPrime, kPrimePrime);
                eDev = only(euler_dev(m, [c, cPrime]));
                eulerValid = eulerValid && abs(eDev < 1e-4);
            end
        end
    end

    if !bcValid
        @warn "Budget constraint violated"
    end
    if !eulerValid
        @warn "Euler equation violated"
    end

    return bcValid  &&  eulerValid
end


# ----------------