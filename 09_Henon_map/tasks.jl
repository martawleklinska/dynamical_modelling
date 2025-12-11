include("ex09.jl")

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
    save("la09_Henon_mapb09/data/task1_xy_a.png", fig)
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
    save("09_Henon_map/data/task2_compare_a.png", fig)
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
    save("09_Henon_map/data/task3_autocorr_a.png", fig)

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
    save("09_Henon_map/data/task4_bifurcation.png", fig)
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
    save("09_Henon_map/data/task5_poincare_compare_a.png", fig)
end

###############################################################
# 6. TASK 
###############################################################
function run_task6()
    println("Running: Task 6 — self-similarity zoom…")

    a_zoom = 1.4
    # xs, ys = henon_iter(0.0, 0.0, a_zoom, b, 50000000)
    println("Loading data...")
    data = readdlm("09_Henon_map/data/task6_zooms.txt")

    xs = data[:, 1]
    ys = data[:, 2]
    
    N_total = length(xs)
    subsample_factor = max(1, N_total ÷ 200000)  
    indices = 1:subsample_factor:N_total
    
    xs_sub = xs[indices]
    ys_sub = ys[indices]
    
    
    fig = Figure(size=(800, 450))
    ax = Axis(fig[1,1], xlabel=L"x", ylabel=L"y", title=L"a=1.4", aspect=DataAspect(), 
             xticklabelsize = 18, yticklabelsize = 18, xlabelsize = 25, ylabelsize = 25, titlesize = 25)
    CairoMakie.scatter!(ax, xs_sub, ys_sub, markersize=0.3, color=:blue, alpha=0.4)
    save("lab09_Henon_map09/data/task6_attractor.png", fig)

    fig1 = zoom(xs, ys; 
            xmin=1.0, xmax=1.35,
            ymin=-0.0, ymax=0.25)
    save("09_Henon_map/data/task6_zoom1.png", fig1)

    fig2 = zoom(xs, ys; 
            xmin=1.080, xmax=1.0807,
            ymin=0.0915, ymax=0.092)
    save("09_Henon_map/data/task6_zoom2.png", fig2)
    
end

# ======================= TASK 7 =============================
function run_ex7()
    D, err, lgε, lgN = do_linear_regression(xs, ys)
    
    X = [lgε ones(length(lgε))]
    β = X \ lgN
    slope = β[1]  
    intercept = β[2]
    
    lgε_line = range(minimum(lgε), maximum(lgε), length=100)
    lgN_fit = slope .* lgε_line .+ intercept
    
    fig = Figure(size = (800, 450))
    ax = Axis(fig[1,1], 
             xlabel = L"\log (1/\varepsilon)", 
             ylabel = L"\log N(\varepsilon)", 
             title = L"\text{Regresja liniowa do wymiaru pudełkowego}",
             xlabelsize = 25, 
             ylabelsize = 25,
             titlesize = 30,
             xticklabelsize = 18,
             yticklabelsize = 18)
    
    scatter!(ax, lgε, lgN, 
            markersize = 12, 
            color = :blue, 
            label = "punkty")
    
    lines!(ax, lgε_line, lgN_fit, 
           color = :red, 
           linewidth = 3,
           label = "fit")

    axislegend(ax, position = :lt, labelsize = 20)

    
    save("09_Henon_map/data/ex7_box_counting_dimension.png", fig)
    # display(fig)
    println("D = $D ± $err")
    return fig
end

# ======================= TASK 8 =============================
function run_ex8()
    D, err, logr, logC = correlation_dimension(xs, ys)
 
    X = [logr ones(length(logr))]
    β = X \ logC
    slope = β[1]      
    intercept = β[2]  
    
    logr_line = range(minimum(logr), maximum(logr), length=100)
    logC_fit = slope .* logr_line .+ intercept
    
    fig = Figure(size = (800, 450))
    ax = Axis(fig[1,1], 
             xlabel = L"\log \varepsilon", 
             ylabel = L"\log C(\varepsilon)", 
             title = L"\text{Regresja liniowa do wymiaru korelacyjnego}",
             xlabelsize = 25, 
             ylabelsize = 25,
             titlesize = 30,
             xticklabelsize = 18,
             yticklabelsize = 18)
    
    scatter!(ax, logr, logC, 
            markersize = 12, 
            color = :blue, 
            label = "punkty")
    
    lines!(ax, logr_line, logC_fit, 
           color = :red, 
           linewidth = 3,
           label = "fit")
    
    axislegend(ax, position = :lt, labelsize = 20)
    
    save("09_Henon_map/data/ex8_correlation_dimension.png", fig)
    # display(fig)
    println("corr: d = $D ± $err")

    return fig
end
