using LinearAlgebra
using Printf
using Plots
using StatsBase

if !isdir("10_multifractal/graphics")
    mkpath("10_multifractal/graphics")
end
## params
const r = 3.569945672          
const N_total = 3_000_000        
const N_discard = 1_000        
const x0 = 0.3   

epsilons = 2.0 .^ (-16:-1:-16)

qs = collect(-5:0.5:5)

logistic(x) = r * x * (1 - x)

function generate_trajectory(x0, N_total, N_discard)
    x = x0
    data = Float64[]
    for n in 1:N_total
        x = logistic(x)
        if n > N_discard
            push!(data, x)
        end
    end
    return data
end
function box_counting(data, epsilons)
    Ns = Float64[]
    for ε in epsilons
        nbins = Int(ceil(1 / ε))
        bins = falses(nbins)
        for x in data
            i = clamp(Int(floor(x / ε)) + 1, 1, nbins)
            bins[i] = true
        end
        push!(Ns, count(bins))
    end
    return Ns
end

function ex1()
    xdata = generate_trajectory(x0, N_total, N_discard)
    N = length(xdata)
    
    Ns = box_counting(xdata, epsilons)
    
    X = log.(1 ./ epsilons)
    Y = log.(Ns)
    
    A = [X ones(length(X))]
    coeff = A \ Y
    D0 = coeff[1]
    residuals = Y - A * coeff
    σD0 = std(residuals) / sqrt(sum((X .- mean(X)).^2))
    
    @printf("\nWymiar pudełkowy D0 = %.5f ± %.5f\n", D0, σD0)
    
    plot(X, Y, seriestype=:scatter, label="data",
         xlabel="ln(1/ε)", ylabel="ln N(ε)",
         title="wymiar pudełkowy", size = (800, 450))
    plot!(X, A * coeff, label="fit")
    savefig("10_multifractal/graphics/ex1.pdf")
end
# ex1()
## ===========================================================
# ex2
# ============================================================
function get_histogram()
    xdata = generate_trajectory(x0, N_total, N_discard)
    p_hist = plot(layout=(length(epsilons), 1), size=(600, 400))
    for (i, ε) in enumerate(epsilons)
        nbins = Int(ceil(1 / ε))
        h = fit(Histogram, xdata, nbins=nbins)
        plot!(p_hist[i], h.weights ./ sum(h.weights),
        label="ε=$(round(ε, sigdigits=2))", xlims = (13000, 14000, 2e04))
    end
    savefig("10_multifractal/graphics/ex2_hist3.pdf")
    # display(p_hist)
end
# get_histogram()
##

# ============================================================
# ex3
# ============================================================

function cantor_theoretical_tau(q, p)
    return -log(p^q + (1-p)^q) / log(3)
end

function cantor_theoretical_Dq(q, p)
    if abs(q - 1) < 1e-6  # q ≈ 1
        q_plus = q + 0.001
        q_minus = q - 0.001
        tau_plus = cantor_theoretical_tau(q_plus, p)
        tau_minus = cantor_theoretical_tau(q_minus, p)
        return (tau_plus - tau_minus) / (q_plus - q_minus)
    else
        return cantor_theoretical_tau(q, p) / (q - 1)
    end
end

function generate_cantor_measures(p, k_max=5)
    measures = Dict()
    epsilons_cantor = []
    
    for k in 1:k_max
        ε = 3.0^(-k)
        push!(epsilons_cantor, ε)
        
        n_intervals = 2^k
        μ = zeros(n_intervals)
        
        for i in 0:(n_intervals-1)
            binary_str = string(i, base=2, pad=k)
            μ_val = 1.0
            for bit_char in binary_str
                if bit_char == '0'
                    μ_val *= p
                else
                    μ_val *= (1-p)
                end
            end
            μ[i+1] = μ_val
        end
        
        measures[ε] = μ
    end
    
    return epsilons_cantor, measures
end

function calculate_Z_cantor(q, epsilons_cantor, measures)
    Z_values = Float64[]
    
    for ε in epsilons_cantor
        μ = measures[ε]
        Z_q_ε = sum(μ.^q)
        push!(Z_values, Z_q_ε)
        
        if abs(q - 1) < 1e-6  
            @printf("Z(1, ε=%.5f) = %.10f\n", ε, Z_q_ε)
        end
    end
    
    return Z_values
end

function fit_tau_cantor(q, epsilons_cantor, Z_values)
    ln_eps = log.(epsilons_cantor)
    ln_Z = log.(Z_values)
    
    A = [ln_eps ones(length(ln_eps))]
    coeff = A \ ln_Z
    tau_numerical = coeff[1]
    
    return tau_numerical
end

function test_cantor_multifractal()    
    p_values = [0.5, 0.25]
    q_test = [-1, 0, 1, 2]
    k_max = 6  
    
    for p in p_values
        epsilons_cantor, measures = generate_cantor_measures(p, k_max)
        
        results_τ = []
        results_Dq = []
        
        for q in q_test
            q_calc = (abs(q - 1) < 1e-6) ? 0.999 : q  
            
            Z_values = calculate_Z_cantor(q_calc, epsilons_cantor, measures)
            τ_numerical = fit_tau_cantor(q_calc, epsilons_cantor, Z_values)
            τ_theoretical = cantor_theoretical_tau(q, p)
            
            Dq_numerical = (abs(q - 1) < 1e-6) ? τ_numerical : τ_numerical / (q - 1)
            Dq_theoretical = cantor_theoretical_Dq(q, p)
            
            push!(results_τ, (q, τ_numerical, τ_theoretical))
            push!(results_Dq, (q, Dq_numerical, Dq_theoretical))
            
        end
        
        plot_cantor_validation(epsilons_cantor, measures, p, q_test)
    end

end

function plot_cantor_validation(epsilons_cantor, measures, p, q_test)
    ln_eps = log.(epsilons_cantor)

    p1 = plot(xlabel="ln ε", ylabel="ln Z(q,ε)",
              title="Cantor: ln Z vs ln ε (p=$p)",
              size=(800,500))

    colors = [:blue, :red, :green, :purple]

    τ_num = Float64[]
    τ_err = Float64[]
    τ_theo = Float64[]

    for (i,q) in enumerate(q_test)
        q_calc = abs(q-1)<1e-6 ? 0.999 : q
        Z = calculate_Z_cantor(q_calc, epsilons_cantor, measures)

        τn, στ = fit_tau_with_error(epsilons_cantor, Z)
        τt = cantor_theoretical_tau(q, p)

        push!(τ_num, τn)
        push!(τ_err, στ)
        push!(τ_theo, τt)

        plot!(p1, ln_eps, log.(Z), seriestype=:scatter,
              color=colors[i], label="q=$q (num)")

        plot!(p1, ln_eps, τt .* ln_eps .+ log(Z[1]) .- τt .* ln_eps[1],
              color=colors[i], linestyle=:dash, label="q=$q (teoria)")
    end

    savefig(p1, "10_multifractal/graphics/cantor_lnZ_p$(p).pdf")
    display(p1)

    p2 = plot(q_test, τ_theo, lw=3, label="τ(q) teoria",
              xlabel="q", ylabel="τ(q)",
              title="Cantor: τ(q) (p=$p)")

    scatter!(q_test, τ_num, yerror=τ_err,
             label="τ(q) numeryczne", ms=6)

    savefig(p2, "10_multifractal/graphics/cantor_tau_p$(p).pdf")
    display(p2)

    Dq_theo = Float64[]
    Dq_num = Float64[]
    
    for (i, q) in enumerate(q_test)
        if abs(q - 1) < 1e-6
            Dq_t = cantor_theoretical_Dq(q, p)
            Dq_n = τ_num[i]  
        else
            Dq_t = τ_theo[i] / (q - 1)
            Dq_n = τ_num[i] / (q - 1)
        end
        push!(Dq_theo, Dq_t)
        push!(Dq_num, Dq_n)
    end
    
    p3 = plot(q_test, Dq_theo, lw=3, label="D(q) teoria",
              xlabel="q", ylabel="D(q)",
              title="Cantor: D(q) (p=$p)")
              
    scatter!(p3, q_test, Dq_num, yerror=τ_err ./ abs.(q_test .- 1),
             label="D(q) numeryczne", ms=6)

    if abs(p-0.5) < 1e-8
        A = [q_test ones(length(q_test))]
        coeff = A \ τ_num
        τ_fit = A * coeff

        residuals = τ_num - τ_fit
        σ = std(residuals)

        @printf("a = %.6f, b = %.6f\n", coeff[1], coeff[2])
        @printf("σ_regresji = %.6e\n", σ)

        plot!(p2, q_test, τ_fit, linestyle=:dash, label="linear fit")
        savefig(p2, "10_multifractal/graphics/cantor_tau_linear_p05.pdf")

        Dq_fit = Float64[]
        for q in q_test
            if abs(q - 1) < 1e-6
                push!(Dq_fit, coeff[1])
            else
                tau_fitted = coeff[1] * q + coeff[2]
                push!(Dq_fit, tau_fitted / (q - 1))
            end
        end
        
        plot!(p3, q_test, Dq_fit, linestyle=:dash, label="linear fit")
        savefig(p3, "10_multifractal/graphics/cantor_Dq_linear_p05.pdf")
    end
    
    savefig(p3, "10_multifractal/graphics/cantor_Dq_p$(p).pdf")
    display(p3)
end


function fit_tau_with_error(epsilons, Z_values)
    x = log.(epsilons)
    y = log.(Z_values)

    A = [x ones(length(x))]
    coeff = A \ y
    τ = coeff[1]

    residuals = y - A * coeff
    σ2 = sum(residuals.^2) / (length(y) - 2)
    cov = σ2 * inv(transpose(A) * A)
    στ = sqrt(cov[1,1])

    return τ, στ
end


test_cantor_multifractal()