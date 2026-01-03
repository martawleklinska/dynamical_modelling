using CairoMakie
using DelimitedFiles
using FilePathsBase
## ====================== DUFFING OSCILLATOR PROJEKT =========================

struct Params
    zeta::Float64
    alpha::Float64
    beta::Float64
    gamma::Float64
    omega::Float64
end
p = [
    Params(0.1, 1.0, 5.0, 0.0, 0.0),
    Params(0.1, -1.0, 5.0, 0.0, 0.0),
    Params(.1, 1.0, -5.0, 0.0, 0.0),
    Params(.1, -1.0, -5.0, 0.0, 0.0),
    Params(3., -1.0, -5.0, 0.0, 0.0),
    Params(3., 1.0, 5.0, 0.0, 0.0),
    Params(-0.1, 1.0, 5.0, 0.0, 0.0),
    Params(-0.1, -1.0, -5.0, 0.0, 0.0),
    Params(-3., -1.0, -5.0, 0.0, 0.0)
    ]
filenames = [
    "ab_pos_zeta_small_",
            "b_pos_a_neg_zeta_small_",
            "a_pos_b_neg_zeta_small_",
            "ab_neg_zeta_small_",
            "ab_neg_zeta_big_",
            "ab_pos_zeta_big_",
            "ab_pos_zeta_neg_",
            "ab_neg_zeta_neg_x",
            "ab_neg_zeta_neg_big_"
            ]
titles = [  
          L"$(\zeta, \; \alpha, \; \beta) = (0.1,\; 1.0,\; 5.0)$",
          L"$(\zeta, \; \alpha, \; \beta) = (0.1,\; -1.0,\; 5.0)$",
          L"$(\zeta, \; \alpha, \; \beta) = (0.1, \;1.0, \;-5.0)$",
          L"$(\zeta, \; \alpha, \; \beta) = (0.1,\; -1.0,\; -5.0)$",
          L"$(\zeta, \; \alpha, \; \beta) = (3., \;-1.0, \;-5.0)$",
          L"$(\zeta, \; \alpha, \; \beta) = (3., \;1.0, \;5.0)$",
          L"$(\zeta, \; \alpha, \; \beta) = (-0.1, \;1.0, \;5.0)$",
          L"$(\zeta, \; \alpha, \; \beta) = (-0.1, \;-1.0, \;-5.0)$",
          L"$(\zeta, \; \alpha, \; \beta) = (-3., \;-1.0, \;-5.0)$"
]

## =======================================================================
# === TRAJECTORIES FOR 9 GIVEN COMBINATIONS OF PARAMS ζ, α, β with γ=0 ===
## =======================================================================
function get_trajectories(p::Params, filename, title, isdisplay)
    init_vals = [
        L"(x_0,\;v_0)=(-1.0, \;0.5)", 
        L"(x_0,\;v_0)=(-1., \;2.0)",
        L"(x_0,\;v_0)=(0.0, \;0.2)", 
        L"(x_0,\;v_0)=(0.0, \;0.5)",
        L"(x_0,\;v_0)=(0.5, \;-1.7)",
        L"(x_0,\;v_0)=(0.5, \;-2.0)", 
        ]
    points = [[-1.0, 0.5], [-1.0, 2.0], [0.0, 0.2], [0.0, 0.5], [0.5, -1.7], [0.5, -2.0]]
    xs = range(-1.2, 1.2, length = 20)
    vs = range(-2., 2., length = 20)

    X = [v for x in xs, v in vs]
    V = [-2.0 * p.zeta * v - p.alpha * x - p.beta * x^3 for x in xs, v in vs]
    mag = sqrt.(X.^2 .+ V.^2) .* 15
    X ./= mag
    V ./= mag
    
    # ============================================================
    fig = Figure(size = (1000, 500))
    ax = Axis(fig[1, 1], xlabel=L"x", ylabel=L"$v$", xlabelsize = 30, limits = ((-1.2, 1.2), (-2.2, 2.)),
    ylabelsize = 30, title=title, titlesize = 30, xticklabelsize = 20, yticklabelsize = 20)
    arrows!(ax, xs, vs, X, V, arrowsize=6, lengthscale = 1.2, linecolor=:gray, linewidth=8, alpha = 0.6)
    
    # ============================================================
    data_dir = "/home/marta/Documents/studia/dynamical_modelling/projekt_Duffing/build/data/"
    begin_file = "duff_$filename"

    files = filter(f -> occursin(begin_file, f) && endswith(f, ".txt"), readdir(data_dir; join=true))

    println("Found $(length(files)) data files:")
    foreach(println, files)

    cm = cgrad(:RdBu_6, 6)

    for (i, file) in enumerate(files)
        data = readdlm(file, '\t', skipstart=1)
        label = init_vals[i]
        if size(data, 2) >= 3
            x_vals = data[:, 2]  
            y_vals = data[:, 3]  
            lines!(ax, x_vals, y_vals, linewidth=4, label = label, alpha = 0.7, color = cm[i])
            scatter!(ax, [points[i][1]], [points[i][2]], marker = :circle, markersize = 20, color = cm[i])
        else
            @warn "cos jest nie tak z $(file) "
        end
    end

    scatter!(ax, 0.0, 0.0, color = :black, markersize = 20, marker = :star5)
    if p.alpha * p.beta < 0
        scatter!(ax, sqrt(-p.alpha/p.beta), 0.0, color = :blue, markersize = 20, marker = :star5)
        scatter!(ax, -sqrt(-p.alpha/p.beta), 0.0, color = :red, markersize = 20, marker = :star5)
    end
    Legend(fig[1,2], ax, labelsize = 25, framevisible = false)
    if isdisplay
        display(fig)
    else 
        save("projekt_Duffing/graphics/$filename.pdf", fig)
    end
    return fig
end

## ============================================================================
# === TIME-DEPENDENCIES FOR 9 GIVEN COMBINATIONS OF PARAMS ζ, α, β with γ=0 ===
## =======================================================================
function get_time_dependencies(p::Params, filename, title, isdisplay)
    init_vals = [
        L"(x_0,\;v_0)=(-1.0, \;0.5)", 
        L"(x_0,\;v_0)=(-1., \;2.0)",
        L"(x_0,\;v_0)=(0.0, \;0.2)", 
        L"(x_0,\;v_0)=(0.0, \;0.5)",
        L"(x_0,\;v_0)=(0.5, \;-1.7)",
        L"(x_0,\;v_0)=(0.5, \;-2.0)", 
        ]
    # ============================================================
    fig = Figure(size = (1000, 500))
    ax = Axis(fig[1, 1], xlabel=L"t", ylabel=L"$x$", xlabelsize = 30, limits = ((-0.1, 10.2), (-3.2, 3.2)),
    ylabelsize = 30, title=title, titlesize = 30, xticklabelsize = 20, yticklabelsize = 20)
    ax2 = Axis(fig[1, 2], xlabel=L"t", ylabel=L"$v$", xlabelsize = 30, limits = ((-0.1, 10.2), (-10.5, 10.7)),
    ylabelsize = 30, title=title, titlesize = 30, xticklabelsize = 20, yticklabelsize = 20)
    
    # ============================================================
    data_dir = "/home/marta/Documents/studia/dynamical_modelling/projekt_Duffing/build/data/"
    begin_file = "duff_$filename"

    files = filter(f -> occursin(begin_file, f) && endswith(f, "00000.txt"), readdir(data_dir; join=true))

    println("Found $(length(files)) data files:")
    foreach(println, files)

    cm = cgrad(:Dark2_6, length(files))
    for (i, file) in enumerate(files)
        data = readdlm(file, '\t', skipstart=1)
        label = init_vals[i]
        if size(data, 2) >= 3
            x_vals = data[:, 2]  
            v_vals = data[:, 3]  
            t_vals = data[:, 1]
            lines!(ax, t_vals, x_vals, linewidth=4, label = label, alpha = 0.7, color = cm[i])
            lines!(ax2, t_vals, v_vals, linewidth=4, label = label, alpha = 0.7, color = cm[i])
        else
            @warn "cos jest nie tak z $(file) "
        end
    end

    Legend(fig[2,1:2], ax, labelsize = 25, orientation = :horizontal, nbanks = 2, framevisible = false)
    if isdisplay
        display(fig)
    else 
        save("projekt_Duffing/graphics/time$filename.pdf", fig)
    end
    return fig
end
## =======================================================================
# === TRAJECTORIES FOR 1 GIVEN COMBINATION OF PARAMS ζ, α, β with γ≠0 ===
## =======================================================================
function sila_wymuszajaca(p::Params, isdisplay)
    init_vals = [
        L"(x_0,\;v_0)=(-1.0, \;0.5)", 
        L"(x_0,\;v_0)=(-1., \;2.0)",
        L"(x_0,\;v_0)=(0.0, \;0.2)", 
        L"(x_0,\;v_0)=(0.0, \;0.5)",
        L"(x_0,\;v_0)=(0.5, \;-1.7)",
        L"(x_0,\;v_0)=(0.5, \;-2.0)", 
        ]
    points = [[-1.0, 0.5], [-1.0, 2.0], [0.0, 0.2], [0.0, 0.5], [0.5, -1.7], [0.5, -2.0]]
    xs = range(-5.0, 5.0, length = 30)
    vs = range(-5., 5., length = 30)

    X = [v for x in xs, v in vs]
    V = [-2.0 * p.zeta * v - p.alpha * x - p.beta * x^3 for x in xs, v in vs]
    mag = sqrt.(X.^2 .+ V.^2) .* 10
    X ./= mag
    V ./= mag
    
    # ============================================================
    fig = Figure(size = (1000, 500))
    ax = Axis(fig[1, 1], xlabel=L"x", ylabel=L"$v$", xlabelsize = 30, limits = ((-4., 4.1), (-5., 5.1)),
    ylabelsize = 30, title=L"\zeta=0.05, \; \alpha=-1.0, \; \beta=5.0, \; \omega=2.0, \;\gamma=2.5", titlesize = 30, xticklabelsize = 20, yticklabelsize = 20)
    arrows!(ax, xs, vs, X, V, arrowsize=8, lengthscale = 2., linecolor=:gray, linewidth=8, alpha = 0.2)
    # ============================================================
    data_dir = "/home/marta/Documents/studia/dynamical_modelling/projekt_Duffing/build/data/"
    # begin_file = "duff_b_pos_aneg_zeta05_gamma25_omega2"
    begin_file = "duff_b_pos_aneg_zeta05_gamma25_omega2"

    files = filter(f -> occursin(begin_file, f) && endswith(f, "00000.txt"), readdir(data_dir; join=true))

    println("Found $(length(files)) data files:")
    foreach(println, files)

    cm = cgrad(:RdBu_6, 6)

    for (i, file) in enumerate(files)
        data = readdlm(file, '\t', skipstart=1)
        label = init_vals[i]
        if size(data, 2) >= 3
            x_vals = data[:, 2]  
            y_vals = data[:, 3]  
            lines!(ax, x_vals, y_vals, linewidth=3, label = label, alpha = 0.7, color = cm[i])
            scatter!(ax, [points[i][1]], [points[i][2]], marker = :circle, markersize = 20, color = cm[i])
        else
            @warn "cos jest nie tak z $(file) "
        end
    end

    Legend(fig[1,2], ax, labelsize = 25, framevisible = false)
    if isdisplay
        display(fig)
    else 
        save("projekt_Duffing/graphics/sila_wymuszajaca_gamma02.pdf", fig)
    end
    return fig
end

## ===========================================================================
# === TIME-DEPENDENCIES FOR 1 GIVEN COMBINATION OF PARAMS ζ, α, β with γ≠0 ===
## ===========================================================================
function sila_wymuszajaca_time(p::Params, isdisplay)
    init_vals = [
        L"(x_0,\;v_0)=(-1.0, \;0.5)", 
        L"(x_0,\;v_0)=(-1., \;2.0)",
        L"(x_0,\;v_0)=(0.0, \;0.2)", 
        L"(x_0,\;v_0)=(0.0, \;0.5)",
        L"(x_0,\;v_0)=(0.5, \;-1.7)",
        L"(x_0,\;v_0)=(0.5, \;-2.0)", 
        ]
    # ============================================================
    fig = Figure(size = (1000, 500))
    ax = Axis(fig[1, 1], xlabel=L"t", ylabel=L"$x$", xlabelsize = 30,# limits = ((-0.1, 10.2), (-1.2, 1.2)),
    ylabelsize = 30, title=L"\zeta=0.05, \; \alpha=-1.0, \; \beta=5.0, \; \omega=2.0, \;\gamma=2.5", titlesize = 30, xticklabelsize = 20, yticklabelsize = 20)
    ax2 = Axis(fig[1, 2], xlabel=L"t", ylabel=L"$v$", xlabelsize = 30,# limits = ((-0.1, 10.2), (-2.5, 2.7)),
    ylabelsize = 30, titlesize = 30, xticklabelsize = 20, yticklabelsize = 20)
    
    # ============================================================
    data_dir = "/home/marta/Documents/studia/dynamical_modelling/projekt_Duffing/build/data/"
    begin_file = "b_pos_aneg_zeta05_gamma25_omega2"
    # begin_file = "ab_pos_zeta05_gamma02_omega10"
    
    files = filter(f -> occursin(begin_file, f) && endswith(f, ".txt"), readdir(data_dir; join=true))
    
    println("Found $(length(files)) data files:")
    foreach(println, files)
    
    cm = cgrad(:Dark2_6, length(files))
    for (i, file) in enumerate(files)
        data = readdlm(file, '\t', skipstart=1)
        label = init_vals[i]
        if size(data, 2) >= 3
            x_vals = data[:, 2]  
            v_vals = data[:, 3]  
            t_vals = data[:, 1]
            lines!(ax, t_vals, x_vals, linewidth=4, label = label, alpha = 0.7, color = cm[i])
            lines!(ax2, t_vals, v_vals, linewidth=4, label = label, alpha = 0.7, color = cm[i])
        else
            @warn "cos jest nie tak z $(file) "
        end
    end
    
    Legend(fig[2,1:2], ax, labelsize = 25, orientation = :horizontal, nbanks = 2, framevisible = false)
    if isdisplay
        display(fig)
    else 
        save("projekt_Duffing/graphics/time-gamma022.pdf", fig)
    end
    return fig
end

## =======================================================================
# === POINCARE MAP WITH PARAMS ζ, α, β with γ≠0 ==========================
## =======================================================================
function plot_poincare(isdisplay)
    filename = "projekt_Duffing/build/data/poincare_analysis_poincare.txt"
    data = readdlm(filename)
    x = data[:, 2]
    v = data[:, 3]
    fig = Figure(size = (800, 500))
    ax = Axis(fig[1, 1], xlabel=L"x", ylabel=L"$v$", xlabelsize = 30,
    ylabelsize = 30, title=L"\zeta=0.15, \; \alpha=-1.0, \; \beta=1.0, \; \gamma=0.36, \;\omega=1.2", titlesize = 30, xticklabelsize = 20, yticklabelsize = 20)
    scatter!(ax, x, v, markersize = 8)
    if isdisplay
        display(fig)
    else 
        save("projekt_Duffing/graphics/poincare_period1.pdf", fig)
    end
    return fig
end
## =======================================================================
# === LAPUNOV EXPONENT ζ, α, β with respect to γ =========================
## =======================================================================
function get_lapunow_exponen(isdisplay)
    data = readdlm("projekt_Duffing/build/data/lyapunov_lyapunov_gamma.txt", comments=true)
    
    x = data[3:end, 1]
    v = data[3:end, 2]
    fig = Figure(size = (1000, 500))
    ax = Axis(fig[1, 1], xlabel=L"\gamma", ylabel=L"$\lambda$", xlabelsize = 30,
    ylabelsize = 30, title=L"\text{wykładnik Lapunowa}", titlesize = 30, xticklabelsize = 20, yticklabelsize = 20)
    lines!(ax, 
    x, v
    )
    hlines!(ax, 0.0, 0.0, color = :gray, linestyle = :dash)
    if isdisplay
        display(fig)
    else 
        save("projekt_Duffing/graphics/lapunow_gamma.pdf", fig)
    end
end

## =======================================================================
# === BIFURCATION DIAGRAM OF PARAMS ζ, α, β with RESPECT TO γ ===
## =======================================================================
function plot_bifur(isdisplay)
    amp_file = "projekt_Duffing/build/data/bifurcation_gamma_scan_gamma.txt"
    
    data_a = readdlm(amp_file, '\t', comments=true)
    gam = data_a[:,1]
    amp = data_a[:,2]
    
    fig2 = Figure(size=(1000,500))
    ax2 = Axis(fig2[1,1], xlabel=L"$\gamma$", ylabel=L"x",
    title=L"Bifurkacja: amplituda do $\gamma$",
    xlabelsize = 30, ylabelsize = 30, titlesize = 30, xticklabelsize = 20, yticklabelsize = 20)
    scatter!(ax2, gam, amp, markersize=1, alpha=0.7)
    if isdisplay
        display(fig2)
    else 
        save("projekt_Duffing/graphics/biffur_gamma.png", fig2)
    end
end
## =======================================================================
# ========================= POINCARE MAP GIF =============================
## =======================================================================
function poincare_gif()
    data = readdlm("projekt_Duffing/build/data/bifurcation_gamma_scan_gamma.txt", skipstart = true)
    γ = data[:, 1]
    x = data[:, 2]
    v = data[:, 3]
    title_obs = Observable("Mapa Poincare")
    fig = Figure(size = (1000, 500))
    ax = Axis(fig[1,1], xlabel = L"x", ylabel = L"v", title = title_obs, xlabelsize = 30, ylabelsize = 30, titlesize = 30, xticklabelsize = 20, yticklabelsize = 20)
    scatter_fig = CairoMakie.scatter!(ax, x[begin:201], v[begin:201], markersize = 4)
    record()
end
## =======================================================================
# === ENERGY ANALYSIS WITH RESPECT TO TIME AND X: ζ, α, β with γ≠0 =======
## =======================================================================

function run_energy_analysis(isdisplay)
    file2 = readdlm("projekt_Duffing/build/data/energy_energy.txt")
    file1 = readdlm("projekt_Duffing/build/data/energy_energy_.txt")
    file = readdlm("projekt_Duffing/build/data/energy_energy_no_zeta.txt")
    fig = Figure(resolution=(1000,800))
    ax2 = Axis(fig[1,1], xlabel=L"$t$", ylabel=L"\text{energia}",
    title=L"\text{Energia w czasie}",
    xlabelsize = 30, ylabelsize = 30, titlesize = 30, xticklabelsize = 20, yticklabelsize = 20)
    ax1 = Axis(fig[2,1], xlabel=L"$x$", ylabel=L"\text{energia}",
    title=L"\text{Energia w przestrzeni}",
    xlabelsize = 30, ylabelsize = 30, titlesize = 30, xticklabelsize = 20, yticklabelsize = 20)
    
    x = file[2:end, 2]
    t = file[2:end, 1]
    e_kinetic = file[2:end, 4]
    e_potential = file[2:end, 5]
    e_input = file[2:end, 6]	
    e_total = file[2:end, 7]
    lines!(ax2, t, e_kinetic, alpha=0.7, label = L"\text{kinetyczna}", linewidth = 3, color = :seagreen)
    lines!(ax2, t, e_potential, alpha=0.7, label = L"\text{potencjalna}", linewidth = 3, color = :royalblue1)
    # lines!(ax2, t, e_input, alpha=0.7, label = L"\text{praca siły wymuszającej}")
    lines!(ax2, t, e_total, alpha=0.7, label = L"\text{całkowita}", linewidth = 3, color = :purple)
    lines!(ax1, x, e_kinetic, alpha=0.7, #label = L"\text{kinetyczna}", 
    linewidth = 3, color = :seagreen)
    lines!(ax1, x, e_potential, alpha=0.7, #label = L"\text{potencjalna}",  
    linewidth = 3, color = :royalblue1)
    # lines!(ax2, t, e_input, alpha=0.7, label = L"\text{praca siły wymuszającej}")
    lines!(ax1, x, e_total, alpha=0.7,# label = L"\text{całkowita}",
     linewidth = 3, color = :purple)

    Legend(fig[3, 1], ax2, orientation = :horizontal, framevisible = false, labelsize = 27)
    
    if isdisplay
        display(fig)
    else 
        save("projekt_Duffing/graphics/energy_zeta0_gamma2.pdf", fig)
    end
end
function run_ode_trajs(display)
    for i in eachindex(p)
        get_trajectories(p[i], filenames[i], titles[i], display)
    end
end
function run_time_ode(display)
    for i in eachindex(p)
        get_time_dependencies(p[i], filenames[i], titles[i], display)
    end
end
function run_example_force_trajs(display)
    p = Params(0.05, -1.0, 0.25, 2.5, 2.0)
    sila_wymuszajaca(p, display)
end
function run_example_time_force(display)
    p = Params(0.05, -1.0, 0.25, 2.5, 2.0)
    sila_wymuszajaca_time(p, display)
end

## =======================================================================
# === POINCARE MAP GIF WITH CHANGING GAMMA ============================
## =======================================================================
function poincare_gamma_gif()
    filename = "projekt_Duffing/build/data/bifurcation_gamma_scan_gamma.txt"
    data = readdlm(filename, '\t', Float64)  
    
    gamma_vals = data[:, 1]
    x_vals = data[:, 2]
    v_vals = data[:, 3]
    
    valid_indices = .!(isnan.(gamma_vals) .| isnan.(x_vals) .| isnan.(v_vals))
    gamma_vals = gamma_vals[valid_indices]
    x_vals = x_vals[valid_indices]
    v_vals = v_vals[valid_indices]
    
    unique_gammas = unique(gamma_vals)
    points_per_gamma = 1000  
    n_gammas = length(unique_gammas)
    
    fig = Figure(size=(800, 600))
    title_obs = Observable("Mapa Poincaré: γ = $(round(unique_gammas[1], digits=3))")
    
    ax = Axis(fig[1,1], 
             xlabel = L"x", ylabel = L"v", 
             title = title_obs,
             xlabelsize = 30, ylabelsize = 30, titlesize = 25,
             xticklabelsize = 20, yticklabelsize = 20,
             limits = ((-1.2, 1.6), (-0.7, 1.2))
            )  
    
    x_obs = Observable(Float64[])
    v_obs = Observable(Float64[])
    
    scatter_plot = scatter!(ax, x_obs, v_obs, 
                           markersize = 3, 
                           color = :blue, 
                           alpha = 0.7)

    gif_filename = "projekt_Duffing/graphics/poincare_gamma_evolution.gif"
    
    record(fig, gif_filename, 1:n_gammas, framerate = 10) do i
        start_idx = (i-1) * points_per_gamma +1
        end_idx = i * points_per_gamma
        
        current_gamma = unique_gammas[i]
        current_x = x_vals[start_idx:end_idx]
        current_v = v_vals[start_idx:end_idx]
        
        title_obs[] = "Mapa Poincaré: γ = $(round(current_gamma, digits=3))"
        
        x_obs[] = current_x
        v_obs[] = current_v
        
    end
    
    return fig
end
## ============== QUANTUM ANALYSIS ======================
function run_tdse()
    fig = Figure(size = (1000, 500))
    title_obs = Observable("Ewolucja czasowa")
    ax1_color = :blue
    ax = Axis(fig[1,1], xlabel = L"x \; (\text{a.u.})", ylabel = L"|\psi(x,\; t)|^2", 
    xlabelsize = 30, ylabelsize = 30, 
    title = title_obs, xticklabelsize = 18, titlesize = 20,
    yticklabelsize = 18,
    leftspinecolor = ax1_color, 
    yticklabelcolor = ax1_color, ylabelcolor = ax1_color, ytickcolor = ax1_color)
    ax2_color = :red
    ax2 = Axis(fig[1,1], ylabel = L"U(x,\; t)\; \text{(a.u.)}", ylabelsize = 30, titlesize = 20,
    yticklabelsize = 18, yaxisposition = :right, rightspinecolor = ax2_color, 
    yticklabelcolor = ax2_color, ylabelcolor = ax2_color, ytickcolor = ax2_color)
    hidexdecorations!(ax2)
    xlims!(ax, -3, 3)
    ylims!(ax, -0.1, 1.5)  
    ylims!(ax2, -0.8, 2.0)
    xlims!(ax2, -3.0, 3.0)
    
    x = Observable(Float64[])
    ψ = Observable(Float64[])
    
    p = Params(0.0, -1.0, 1.0, 0.3, 1.0)
    V(x_val, t) = 0.5 * p.alpha * x_val^2 + 0.25 * p.beta * x_val^4 - x_val * p.gamma * cos(p.omega * t)
    V_obs = Observable(Float64[])
    
    scatter_plot = CairoMakie.lines!(ax, x, ψ, linewidth=3, color=:blue)
    potential_plot = CairoMakie.lines!(ax2, x, V_obs, linewidth=3, color=:red, linestyle=:dash)
    
    # ============================================================
    data_dir = "projekt_Duffing/build/data/quantum"
    begin_file = "duffing_gauss2_step_"

    files = filter(f -> occursin(begin_file, f) && endswith(f, ".txt"), readdir(data_dir; join=true))
    
    sort!(files, by = f -> parse(Int, match(r"step_(\d+)", f).captures[1]))

    println("Found $(length(files)) data files")
    
    obs_file = joinpath(data_dir, "duffing_gauss2_observables.txt")
    obs_data = readdlm(obs_file, '\t', skipstart=1)
    times = obs_data[:, 2]  

    
    filename = "projekt_Duffing/graphics/tdse_evolution_gauss2.gif"
    record(fig, filename, 1:length(files), framerate = 50) do i
        time = times[i]
        title_obs[] = "Ewolucja czasowa t=$(round(time, digits=3))"
        
        data = readdlm(files[i], '\t', skipstart=1)
        x[] = data[:, 1]  
        ψ[] = data[:, 2]
        V_obs[] = [V(x_val, time) for x_val in data[:, 1]]

        if i % 10 == 0
            println("Frame $i: t=$time, max(|ψ|²)=$(maximum(ψ[]))")
        end
    end
    
    println("Animation saved to $filename")
end
function get_expectation_values(isdisplay::Bool)
    data = readdlm("projekt_Duffing/build/data/quantum/duffing_gauss1_observables.txt", skipstart = 1)
    t = data[:, 2]
    x = data[:, 3]
    p = data[:, 4]
    E = data[:, 5]
    fig = Figure(size = (1000, 500))
    ax1_color = :royalblue1
    ax = Axis(fig[1,1], xlabel = L"t\; (\text{a.u.})", ylabel = L"\langle x\rangle\; (\text{a.u.})", 
    xlabelsize = 40, ylabelsize = 40, xticklabelsize = 30, yticklabelsize = 30,
    leftspinecolor = ax1_color, yaxisposition = :left, 
    yticklabelcolor = ax1_color, ylabelcolor = ax1_color, ytickcolor = ax1_color)
    ax2_color = :crimson
    ax2 = Axis(fig[1,1], ylabel = L"\langle p \rangle\; \text{(a.u.)}", ylabelsize = 40, titlesize = 30,
    yticklabelsize = 30, yaxisposition = :right, rightspinecolor = ax2_color, 
    yticklabelcolor = ax2_color, ylabelcolor = ax2_color, ytickcolor = ax2_color)
    hidexdecorations!(ax2)
    # ylims!(ax2, 8, 10)
    
    lines!(ax, t, x, color = ax1_color, linewidth = 4)
    lines!(ax2, t, p, color = ax2_color, linewidth = 4)
    # lines!(ax2, t, E, color = :purple, linewidth = 4)
    if isdisplay
        display(fig)
    else 
        save("projekt_Duffing/graphics/QM_xp_exp_val2.pdf", fig)
    end
end

# =========== energy =========
function get_exp_energy(is_display::Bool)
    data = readdlm("projekt_Duffing/build/data/quantum/duffing_gauss1_observables.txt", skipstart = 1)
    data2 = readdlm("projekt_Duffing/build/data/quantum/duffing_gauss2_observables.txt", skipstart = 1)
    t = data[:, 2]

    E = data[:, 5]
    E2 = data2[:, 5]
    fig = Figure(size = (1000, 500))
    ax1_color = :black
    ax = Axis(fig[1,1], xlabel = L"t\; (\text{a.u.})", ylabel = L"\langle E\rangle\; (\text{a.u.})", 
    xlabelsize = 35, ylabelsize = 35, xticklabelsize = 25, yticklabelsize = 25,
    leftspinecolor = ax1_color, yaxisposition = :left, 
    yticklabelcolor = ax1_color, ylabelcolor = ax1_color, ytickcolor = ax1_color)
    
    lines!(ax, t, E, color = :maroon3, linewidth = 4, label = L"x_0=-1\;\mathrm{a.u.}")
    lines!(ax, t, E2, color = :seagreen3, linewidth = 4, label = L"x_0=1\;\mathrm{a.u.}")
    axislegend(ax, labelsize = 30, framevisible = false, position = :rb)
    if is_display
        display(fig)
    else 
        save("projekt_Duffing/graphics/QM_E_exp_val2.pdf", fig)
    end
end