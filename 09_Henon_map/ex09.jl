using CairoMakie, DelimitedFiles

function henon_step(x::T, y::T, a::T, b::T) where {T<:Real}
    xn = one(T) - a*x^2 + y
    yn = b*x
    return xn, yn
end

function henon_iter(x0, y0, a, b, N)
    xs = Vector{typeof(x0)}(undef, N+1)
    ys = Vector{typeof(y0)}(undef, N+1)
    xs[1] = x0
    ys[1] = y0

    for n in 1:N
        xs[n+1], ys[n+1] = henon_step(xs[n], ys[n], a, b)
    end
    return xs, ys
end

function autocorr(x, maxlag)
    N = length(x)
    μ = mean(x)
    σ2 = var(x)

    C = zeros(Float64, maxlag+1)
    for k in 0:maxlag
        C[k+1] = sum((x[1:N-k] .- μ) .* (x[1+k:N] .- μ)) / ((N-k)*σ2)
    end
    return C
end

function bifurcation_diagram(a_min, a_max, n_a; b=0.3, Ntras=2000, Nsample=200)
    avals = range(a_min, a_max, length=n_a)
    xplot = Float64[]
    yplot = Float64[]

    for a in avals
        x = 0.1
        y = 0.0
        for i in 1:Ntras
            x, y = henon_step(x, y, a, b)
        end
        for i in 1:Nsample
            x, y = henon_step(x, y, a, b)
            push!(xplot, a)
            push!(yplot, x)
        end
    end

    return xplot, yplot
end

function poincare_section(a, b, x0, y0, N)
    xs, ys = henon_iter(x0, y0, a, b, N)
    return xs, ys
end

function zoom(xs, ys; 
                     xmin, xmax, ymin, ymax)
    mask = (xs .>= xmin) .& (xs .<= xmax) .&
           (ys .>= ymin) .& (ys .<= ymax)
    
    xs_zoom = xs[mask]
    ys_zoom = ys[mask]
    
    if length(xs_zoom) > 50000
        subsample = 1:max(1, length(xs_zoom) ÷ 50000):length(xs_zoom)
        xs_zoom = xs_zoom[subsample]
        ys_zoom = ys_zoom[subsample]
    end

    fig = Figure(size = (800, 450))
    ax = Axis(fig[1,1], xlabel = L"x", ylabel = L"y", title = L"\text{Zoom}", 
             ylabelsize = 25, xlabelsize = 25, titlesize = 25, 
             xticklabelsize = 18, yticklabelsize = 18, aspect = DataAspect())
    
    CairoMakie.scatter!(ax, xs_zoom, ys_zoom, markersize=1.5, color=:blue, alpha=0.6)
    
    xlims!(ax, xmin, xmax)
    ylims!(ax, ymin, ymax)
    
    return fig
end

# ================== ex7 and 8 ===============

function do_linear_regression(xs, ys; scales = 2:10)
    logeps = Float64[]
    logN   = Float64[]

    xmin, _ = extrema(xs)
    ymin, _ = extrema(ys)

    for k in scales
        ε = 2.0^-k

        boxes = Set{Tuple{Int,Int}}()

        for (x,y) in zip(xs, ys)
            ix = floor(Int,(x - xmin)/ε)
            iy = floor(Int,(y - ymin)/ε)# floor to change the cont into disc
            push!(boxes, (ix, iy))
        end

        push!(logeps, log(1/ε))
        push!(logN, log(length(boxes)))
    end

    X = [logeps ones(length(logeps))]
    β = X \ logN
    slope = β[1]

    residuals = logN .- (X*β)
    σ2 = sum(residuals.^2) / (length(logeps) - 2)
    cov = σ2 * inv(X'X)
    σ_slope = sqrt(cov[1,1])

    return slope, σ_slope, logeps, logN
end

function correlation_dimension(xs, ys; scales = 4:12)
    pts = collect(zip(xs, ys))
    N = length(pts)

    logr = Float64[]
    logC = Float64[]

    for k in scales
        r = 10.0^(-k/4)
        count = 0

        for i in 1:N
            for j in i+1:N
                dx = pts[i][1] - pts[j][1]
                dy = pts[i][2] - pts[j][2]
                if dx^2 + dy^2 < r^2
                    count += 1
                end
            end
        end

        C = 2*count / (N*(N-1))
        push!(logr, log(r))
        push!(logC, log(C))
    end

    X = [logr ones(length(logr))]
    β = X \ logC
    slope = β[1]

    residuals = logC .- (X*β)
    σ2 = sum(residuals.^2) / (length(logr) - 2)
    cov = σ2 * inv(X'X)
    σ_slope = sqrt(cov[1,1])

    return slope, σ_slope, logr, logC
end
