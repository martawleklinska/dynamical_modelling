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
function get_trajectories(p::Params, filename, title)
    init_vals = [
        L"(x_0,\;v_0)=(-1.0, \;0.5)", 
        L"(x_0,\;v_0)=(-1., \;2.0)",
        L"(x_0,\;v_0)=(0.0, \;0.2)", 
        L"(x_0,\;v_0)=(0.0, \;0.5)",
        L"(x_0,\;v_0)=(0.5, \;-1.7)",
        L"(x_0,\;v_0)=(0.5, \;-2.0)", 
        ]
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
        else
            @warn "cos jest nie tak z $(file) "
        end
    end

    scatter!(ax, 0.0, 0.0, color = :black, markersize = 20, marker = :star5)
    if p.alpha * p.beta < 0
        scatter!(ax, sqrt(-p.alpha/p.beta), 0.0, color = :blue, markersize = 20, marker = :star5)
        scatter!(ax, -sqrt(-p.alpha/p.beta), 0.0, color = :red, markersize = 20, marker = :star5)
    end
    Legend(fig[1,2], ax, labelsize = 25)
    display(fig)
    # save("projekt_Duffing/graphics/$filename.pdf", fig)
    return fig
end

## ============================================================================
# === TIME-DEPENDENCIES FOR 9 GIVEN COMBINATIONS OF PARAMS ζ, α, β with γ=0 ===
## =======================================================================
function get_time_dependencies(p::Params, filename, title)
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
    display(fig)
    # save("projekt_Duffing/graphics/time$filename.pdf", fig)
    return fig
end
## =======================================================================
# === TRAJECTORIES FOR 1 GIVEN COMBINATION OF PARAMS ζ, α, β with γ≠0 ===
## =======================================================================
function sila_wymuszajaca(p::Params)
    init_vals = [
        L"(x_0,\;v_0)=(-1.0, \;0.5)", 
        L"(x_0,\;v_0)=(-1., \;2.0)",
        L"(x_0,\;v_0)=(0.0, \;0.2)", 
        L"(x_0,\;v_0)=(0.0, \;0.5)",
        L"(x_0,\;v_0)=(0.5, \;-1.7)",
        L"(x_0,\;v_0)=(0.5, \;-2.0)", 
        ]
    xs = range(-2.2, 2.2, length = 20)
    vs = range(-3., 3., length = 20)

    X = [v for x in xs, v in vs]
    V = [-2.0 * p.zeta * v - p.alpha * x - p.beta * x^3 for x in xs, v in vs]
    mag = sqrt.(X.^2 .+ V.^2) .* 15
    X ./= mag
    V ./= mag
    
    # ============================================================
    fig = Figure(size = (1000, 500))
    ax = Axis(fig[1, 1], xlabel=L"x", ylabel=L"$v$", xlabelsize = 30, limits = ((-2., 1.9), (-3., 2.8)),
    ylabelsize = 30, title=L"\zeta=0.05, \; \alpha=1.0, \; \beta=5.0, \; \omega=1.0, \;\gamma=0.2", titlesize = 30, xticklabelsize = 20, yticklabelsize = 20)
    arrows!(ax, xs, vs, X, V, arrowsize=8, lengthscale = 1.5, linecolor=:gray, linewidth=8, alpha = 0.3)
    # ============================================================
    data_dir = "/home/marta/Documents/studia/dynamical_modelling/projekt_Duffing/build/data/"
    # begin_file = "duff_b_pos_aneg_zeta05_gamma25_omega2"
    begin_file = "duff_ab_pos_zeta05_gamma02_omega10"

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
            lines!(ax, x_vals, y_vals, linewidth=4, label = label, alpha = 0.7, color = cm[i])
        else
            @warn "cos jest nie tak z $(file) "
        end
    end

    Legend(fig[1,2], ax, labelsize = 25)
    display(fig)
    # save("projekt_Duffing/graphics/sila_wymuszajaca_gamma02.pdf", fig)
    return fig
end

## ===========================================================================
# === TIME-DEPENDENCIES FOR 1 GIVEN COMBINATION OF PARAMS ζ, α, β with γ≠0 ===
## ===========================================================================
function sila_wymuszajaca_time(p::Params)
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
    ylabelsize = 30, title=L"\zeta=0.05, \; \alpha=1.0, \; \beta=1.0, \; \omega=1.0, \;\gamma=0.2", titlesize = 30, xticklabelsize = 20, yticklabelsize = 20)
    ax2 = Axis(fig[1, 2], xlabel=L"t", ylabel=L"$v$", xlabelsize = 30,# limits = ((-0.1, 10.2), (-2.5, 2.7)),
    ylabelsize = 30, titlesize = 30, xticklabelsize = 20, yticklabelsize = 20)
    
    # ============================================================
    data_dir = "/home/marta/Documents/studia/dynamical_modelling/projekt_Duffing/build/data/"
    # begin_file = "b_pos_aneg_zeta05_gamma25_omega2"
    begin_file = "ab_pos_zeta05_gamma02_omega10"
    
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
    display(fig)
    # save("projekt_Duffing/graphics/time-gamma022.pdf", fig)
    return fig
end

## =======================================================================
# === POINCARE MAP WITH PARAMS ζ, α, β with γ≠0 ==========================
## =======================================================================
function plot_poincare(filename::String)
    data = readdlm(filename, comments=true)
    
    x = data[:, 2]
    v = data[:, 3]
    fig = Figure(size = (800, 500))
    ax = Axis(fig[1, 1], xlabel=L"x", ylabel=L"$v$", xlabelsize = 30,# limits = ((-0.1, 10.2), (-1.2, 1.2)),
    ylabelsize = 30, title=L"\zeta=0.1, \; \alpha=-1.0, \; \beta=0.25, \; \omega=2.0, \;\gamma=2.5", titlesize = 30, xticklabelsize = 20, yticklabelsize = 20)
    scatter!(ax, 
    x, v, markersize = 4
    )
    display(fig)
    # save("projekt_Duffing/graphics/poincare_chaos1.pdf", fig)
    return fig
end

## =======================================================================
# === LAPUNOV EXPONENT ζ, α, β with respect to γ =========================
## =======================================================================
function get_lapunow_exponen()
    data = readdlm("projekt_Duffing/build/data/lyapunov_lyapunov_gamma.txt", comments=true)
    
    x = data[3:end, 1]
    v = data[3:end, 2]
    fig = Figure(size = (1000, 500))
    ax = Axis(fig[1, 1], xlabel=L"\gamma", ylabel=L"$\lambda$", xlabelsize = 30,# limits = ((-0.1, 10.2), (-1.2, 1.2)),
    ylabelsize = 30, title=L"\text{wykładnik Lapunowa}", titlesize = 30, xticklabelsize = 20, yticklabelsize = 20)
    lines!(ax, 
    x, v
    )
    hlines!(ax, 0.0, 0.0, color = :gray, linestyle = :dash)
    display(fig)
    # save("projekt_Duffing/graphics/lapunow_gamma.pdf", fig)
end

## =======================================================================
# === BIFURCATION DIAGRAM OF PARAMS ζ, α, β with RESPECT TO γ ===
## =======================================================================
function plot_bifur()
    amp_file = "projekt_Duffing/build/data/bifurcation_gamma_scan_gamma.txt"
    
    data_a = readdlm(amp_file, '\t', comments=true)
    gam = data_a[:,1]
    amp = data_a[:,2]
    
    fig2 = Figure(resolution=(1000,500))
    ax2 = Axis(fig2[1,1], xlabel=L"$\gamma$", ylabel=L"x",
    title=L"Bifurkacja: amplituda do $\gamma$",
    xlabelsize = 30, ylabelsize = 30, titlesize = 30, xticklabelsize = 20, yticklabelsize = 20)
    scatter!(ax2, gam, amp, markersize=1, alpha=0.7)
    display(fig2)
    # save("projekt_Duffing/graphics/biffur_gamma.png", fig2)
end
## =======================================================================
# === ENERGY ANALYSIS WITH RESPECT TO TIME AND X: ζ, α, β with γ≠0 =======
## =======================================================================

function run_energy_analysis()
    file = readdlm("projekt_Duffing/build/data/energy_energy.txt")
    t = file[2:end, 2]
    x = file[2:end, 1]
    e_kinetic = file[2:end, 4]
    e_potential = file[2:end, 5]
    e_input = file[2:end, 6]	
    e_total = file[2:end, 7]
    fig = Figure(resolution=(1000,800))
    ax2 = Axis(fig[1,1], xlabel=L"$t$", ylabel=L"\text{energia}",
    title=L"\text{Energia w czasie}",
    xlabelsize = 30, ylabelsize = 30, titlesize = 30, xticklabelsize = 20, yticklabelsize = 20)
    ax1 = Axis(fig[2,1], xlabel=L"$t$", ylabel=L"\text{energia}",
    title=L"\text{Energia w przestrzeni}",
    xlabelsize = 30, ylabelsize = 30, titlesize = 30, xticklabelsize = 20, yticklabelsize = 20)
    
    lines!(ax2, t, e_kinetic, alpha=0.7, label = L"\text{kinetyczna}", linewidth = 3)
    lines!(ax2, t, e_potential, alpha=0.7, label = L"\text{potencjalna}", linewidth = 3)
    # lines!(ax2, t, e_input, alpha=0.7, label = L"\text{praca siły wymuszającej}")
    lines!(ax2, t, e_total, alpha=0.7, label = L"\text{całkowita}", linewidth = 3, color = :purple)
    lines!(ax1, x, e_kinetic, alpha=0.7, label = L"\text{kinetyczna}", linewidth = 3)
    lines!(ax1, x, e_potential, alpha=0.7, label = L"\text{potencjalna}", linewidth = 3)
    # lines!(ax2, t, e_input, alpha=0.7, label = L"\text{praca siły wymuszającej}")
    lines!(ax1, x, e_total, alpha=0.7, label = L"\text{całkowita}", linewidth = 3, color = :purple)
    
    Legend(fig[3, 1], ax2, orientation = :horizontal, framevisible = false, labelsize = 27)
    display(fig)
    # save("projekt_Duffing/graphics/energy_analysis.pdf", fig2)
end
function run_ode_trajs()
    for i in eachindex(p)
        get_trajectories(p[i], filenames[i], titles[i])
    end
end
function run_time_ode()
    for i in eachindex(p)
        get_time_dependencies(p[i], filenames[i], titles[i])
    end
end
function run_example_force_trajs()
    p = Params(0.05, 1.0, 5.0, 0.2, 1.0)
    sila_wymuszajaca(p)
end
function run_example_time_force()
    sila_wymuszajaca_time(p)
end
## ============== QUANTUM ANALYSIS ======================
function run_tdse()
    fig = Figure(size = (1000, 500))
    title_obs = Observable("Ewolucja czasowa")
    ax1_color = :blue
    ax = Axis(fig[1,1], xlabel = L"x", ylabel = L"|\psi(x,\; t)|^2", 
    xlabelsize = 30, ylabelsize = 30, 
    title = title_obs, xticklabelsize = 18, titlesize = 20,
    yticklabelsize = 18,
    leftspinecolor = ax1_color, 
    yticklabelcolor = ax1_color, ylabelcolor = ax1_color, ytickcolor = ax1_color)
    ax2_color = :red
    ax2 = Axis(fig[1,1], ylabel = L"U(x,\; t) (a.u.)", ylabelsize = 30, titlesize = 20,
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
    begin_file = "duffing_step_"

    files = filter(f -> occursin(begin_file, f) && endswith(f, ".txt"), readdir(data_dir; join=true))
    
    sort!(files, by = f -> parse(Int, match(r"step_(\d+)", f).captures[1]))

    println("Found $(length(files)) data files")
    
    obs_file = joinpath(data_dir, "duffing_observables.txt")
    obs_data = readdlm(obs_file, '\t', skipstart=1)
    times = obs_data[:, 2]  

    
    filename = "projekt_Duffing/graphics/tdse_evolution.gif"
    record(fig, filename, 1:length(files), framerate = 10) do i
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
