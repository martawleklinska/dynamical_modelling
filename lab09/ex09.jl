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
    println("  Filtering points in zoom region...")
    mask = (xs .>= xmin) .& (xs .<= xmax) .&
           (ys .>= ymin) .& (ys .<= ymax)
    
    xs_zoom = xs[mask]
    ys_zoom = ys[mask]
    
    println("  Found $(length(xs_zoom)) points in zoom region")
    
    if length(xs_zoom) > 50000
        subsample = 1:max(1, length(xs_zoom) ÷ 50000):length(xs_zoom)
        xs_zoom = xs_zoom[subsample]
        ys_zoom = ys_zoom[subsample]
        println("  Subsampled to $(length(xs_zoom)) points for plotting")
    end
    
    println("  Creating plot...")
    fig = Figure(size = (800, 450))
    ax = Axis(fig[1,1], xlabel = L"x", ylabel = L"y", title = L"\text{Zoom}", 
             ylabelsize = 25, xlabelsize = 25, titlesize = 25, 
             xticklabelsize = 18, yticklabelsize = 18, aspect = DataAspect())
    
    CairoMakie.scatter!(ax, xs_zoom, ys_zoom, markersize=1.5, color=:blue, alpha=0.6)
    
    xlims!(ax, xmin, xmax)
    ylims!(ax, ymin, ymax)
    
    return fig
end



as = [0.50, 1.10, 1.25, 1.40]
b = 0.3
N = 2000

###############################################################
# 1. TASK 
###############################################################
function run_task1()
    fig = Figure(size = (1200, 2000))
    for (i, a) in enumerate(as)
        xs, ys = henon_iter(0.0, 0.0, a, b, N)
        ax = Axis(fig[i,1], title="x_n for a=$a", xlabel=L"n", ylabel=L"x_n", xlabelsize = 35, ylabelsize = 35, xticklabelsize = 18, yticklabelsize = 28, titlesize = 35)
        CairoMakie.scatter!(ax, 1:length(xs), xs, color=:blue, markersize=3.)

        ax = Axis(fig[i,2], title="y_n for a=$a", xlabel=L"n", ylabel=L"y_n", xlabelsize = 35, ylabelsize = 35, xticklabelsize = 28, yticklabelsize = 28, titlesize = 35)
        CairoMakie.scatter!(ax, 1:length(ys), ys, color=:red, markersize=3.)
    end
    save("lab09/data/task1_xy_a.png", fig)
end
###############################################################
# 2. TASK 
###############################################################
function run_task2()
    fig = Figure(size = (1300, 2000))
    for (i, a) in enumerate(as)
        xs64, _ = henon_iter(0.0, 0.0, a, b, N)
        xs32, _ = henon_iter(0.0f0, 0.0f0, Float32(a), Float32(b), N)

        ax = Axis(fig[i,1], title="Float64 vs Float32, a=$a", xlabel=L"$n$", ylabel=L"$x_n$", 
        xlabelsize = 35, ylabelsize = 35, titlesize = 35, xticklabelsize = 28, yticklabelsize = 28)
        CairoMakie.scatter!(ax, 1:length(xs64), xs64, label="Float64", color=:blue, markersize=5.)
        CairoMakie.scatter!(ax, 1:length(xs32), Float64.(xs32), label="Float32", color=:red, markersize=5.)
        Legend(fig[5, 1], ax, labelsize = 35, orientation = :horizontal, framevisible = false)
    end
    save("lab09/data/task2_compare_a.png", fig)
end

###############################################################
# 3. TASK 
###############################################################
function run_task3()
    fig = Figure(size = (1300, 2000))
    for (i, a) in enumerate(as)
        xs, _ = henon_iter(0.0, 0.0, a, b, N)
        C = autocorr(xs, 200)

        ax = Axis(fig[i,1], xlabel=L"k", ylabel=L"C(k)", title="Autokorelacja przy a=$a (float64)", 
        xlabelsize = 35, ylabelsize = 35, titlesize = 35, xticklabelsize = 28, yticklabelsize = 28)
        if i < 3
            xlims!(ax, 0.0, 10.)
        elseif i == 3
            xlims!(ax, 0.0, 50)
        else
            xlims!(ax, 0.0, 200)
        end
        scatter!(ax, 0:length(C)-1, C, color=:green, markersize=10.)
        hlines!(ax, [0], color=:black, linestyle=:dash, alpha=0.5)
    end
    save("lab09/data/task3_autocorr_a.png", fig)

end

###############################################################
# 4. TASK 
###############################################################
function run_task4()

    xplot, yplot = bifurcation_diagram(0.0, 1.5, 2000)

    fig = Figure(size=(800, 450))
    ax = Axis(fig[1,1], xlabel=L"a", ylabel=L"x", title="Diagram bifurkacji", xticklabelsize = 18, 
    yticklabelsize = 18, xlabelsize = 25, ylabelsize = 25, titlesize = 25)
    ylims!(ax, -1.3, 1.8)
    CairoMakie.scatter!(ax, xplot, yplot, markersize=0.5, color=:black, alpha=0.7)
    save("lab09/data/task4_bifurcation.png", fig)
end

###############################################################
# 5. TASK 
###############################################################
function run_task5()
    fig = Figure(size=(1300, 2000))
    for (i, a) in enumerate(as)
        xs64, ys64 = henon_iter(0.0, 0.0, a, b, N)
        xs32, ys32 = henon_iter(0.0f0, 0.0f0, Float32(a), Float32(b), N)

        ax = Axis(fig[i,1], xlabel=L"x", ylabel=L"y", title="Poincaré for a=$a", aspect=DataAspect(), xticklabelsize = 28, yticklabelsize = 28, xlabelsize = 35, ylabelsize = 35, titlesize = 35)
        CairoMakie.scatter!(ax, xs64, ys64, markersize=4, label="Float64", color=:blue, alpha=0.7)
        CairoMakie.scatter!(ax, Float64.(xs32), Float64.(ys32), markersize=4, label="Float32", color=:red, alpha=0.7)
        Legend(fig[5,1], ax, framevisible = false, orientation = :horizontal, labelsize = 35)
    end
    save("lab09/data/task5_poincare_compare_a.png", fig)
end

###############################################################
# 6. TASK 
###############################################################
function run_task6()
    println("Running: Task 6 — self-similarity zoom…")

    a_zoom = 1.4
    # xs, ys = henon_iter(0.0, 0.0, a_zoom, b, 50000000)
    println("Loading data...")
    data = readdlm("lab09/data/task6_zooms.txt")
    println("Data loaded. Size: $(size(data))")
    
    xs = data[:, 1]
    ys = data[:, 2]
    
    N_total = length(xs)
    subsample_factor = max(1, N_total ÷ 200000)  
    indices = 1:subsample_factor:N_total
    
    xs_sub = xs[indices]
    ys_sub = ys[indices]
    
    println("Subsampled to $(length(xs_sub)) points (factor: $subsample_factor)")
    println("Plotting full attractor...")
    
    fig = Figure(size=(800, 450))
    ax = Axis(fig[1,1], xlabel=L"x", ylabel=L"y", title=L"a=1.4", aspect=DataAspect(), 
             xticklabelsize = 18, yticklabelsize = 18, xlabelsize = 25, ylabelsize = 25, titlesize = 25)
    CairoMakie.scatter!(ax, xs_sub, ys_sub, markersize=0.3, color=:blue, alpha=0.4)
    save("lab09/data/task6_attractor.png", fig)
    println("Full attractor saved.")

    fig1 = zoom(xs, ys; 
            xmin=1.0, xmax=1.35,
            ymin=-0.0, ymax=0.25)
    save("lab09/data/task6_zoom1.png", fig1)

    fig2 = zoom(xs, ys; 
            xmin=1.080, xmax=1.0807,
            ymin=0.0915, ymax=0.092)
    save("lab09/data/task6_zoom2.png", fig2)
    
end

function run_all_tasks()
    # run_task1()
    # run_task2()
    # run_task3()
    # run_task4()
    # run_task5()
    run_task6()
end
run_all_tasks()


## ================== ex7 and 8 ===============

xs, ys = henon_iter(0.0, 0.0, 1.4, b, 50000)

function box_count_dimension(xs, ys; scales = 3:10)
    logeps = Float64[]
    logN   = Float64[]

    xmin, xmax = extrema(xs)
    ymin, ymax = extrema(ys)

    for k in scales
        ε = 2.0^-k
        nx = ceil(Int, (xmax - xmin) / ε)
        ny = ceil(Int, (ymax - ymin) / ε)

        boxes = Set{Tuple{Int,Int}}()

        for (x,y) in zip(xs, ys)
            ix = floor(Int,(x - xmin)/ε)
            iy = floor(Int,(y - ymin)/ε)
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

D, err, lgε, lgN = box_count_dimension(xs, ys)
println("Box-counting dimension = $D ± $err")
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
D, err, _, _ = correlation_dimension(xs, ys)
println("Correlation dimension = $D ± $err")
