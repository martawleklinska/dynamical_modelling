using CairoMakie
using DelimitedFiles

## ====================== DUFFING OSCILLATOR PROJEKT =========================

# test params
struct Params
    zeta::Float64
    alpha::Float64
    beta::Float64
    gamma::Float64
    omega::Float64
end

function get_trajectories(p::Params, filename, title)
    # p = Params(-3., -1.0, -5.0, 0.0, 0.0)
    # L"$(\zeta, \; \alpha, \; \beta) = (-3., \; -1.0, \; -5.0)$"

    init_vals = [L"(x_0,\;v_0)=(0.5, \;-2.0)", L"(x_0,\;v_0)=(-1., \;2.0)",
     L"(x_0,\;v_0)=(-1.0, \;0.5)", L"(x_0,\;v_0)=(0.5, \;-1.7)",
     L"(x_0,\;v_0)=(0.0, \;0.1)", L"(x_0,\;v_0)=(0.0, \;-0.5)"]
    xs = range(-1.2, 1.2, length = 20)
    vs = range(-2., 2., length = 20)

    X = [v for x in xs, v in vs]
    V = [-2.0 * p.zeta * v + p.alpha * x - p.beta * x^3 for x in xs, v in vs]
    mag = sqrt.(X.^2 .+ V.^2) .* 20
    X ./= mag
    V ./= mag
    
    # ============================================================
    fig = Figure(size = (1000, 500))
    ax = Axis(fig[1, 1], xlabel=L"x", ylabel=L"$v$", xlabelsize = 30, limits = ((-1.2, 1.2), (-2.2, 2.)),
    ylabelsize = 30, title=title, titlesize = 30, xticklabelsize = 20, yticklabelsize = 20)
    arrows!(ax, xs, vs, X, V, arrowsize=8, lengthscale = 1.5, linecolor=:gray, linewidth=8, alpha = 0.6)
    # ============================================================
    data_dir = "/home/marta/Documents/studia/dynamical_modelling/projekt/build/data/"
    begin_file = "duff_$filename"

    files = filter(f -> occursin(begin_file, f) && endswith(f, ".txt"), readdir(data_dir; join=true))

    println("Found $(length(files)) data files:")
    foreach(println, files)

    # cm = cgrad(:Blues, 6)
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
    if p.alpha * p.beta > 0
        scatter!(ax, sqrt(p.alpha/p.beta), 0.0, color = :blue, markersize = 20, marker = :star5)
        scatter!(ax, -sqrt(p.alpha/p.beta), 0.0, color = :red, markersize = 20, marker = :star5)
    end
    Legend(fig[1,2], ax, labelsize = 25)
    # display(fig)
    save("projekt/graphics/$filename.pdf", fig)
    return fig
end

p = [Params(0.1, 1.0, 5.0, 0.0, 0.0),
    Params(0.1, -1.0, 5.0, 0.0, 0.0),
    Params(.1, 1.0, -5.0, 0.0, 0.0),
    Params(.1, -1.0, -5.0, 0.0, 0.0),
    # Params(3., -1.0, -5.0, 0.0, 0.0),
    # Params(3., 1.0, 5.0, 0.0, 0.0),
    Params(-0.1, 1.0, 5.0, 0.0, 0.0),
    Params(-0.1, -1.0, -5.0, 0.0, 0.0),
    # Params(-3., -1.0, -5.0, 0.0, 0.0)
    ]
filenames = ["ab_pos_zeta_small_",
            "b_pos_a_neg_zeta_small_",
            "a_pos_b_neg_zeta_small_",
            "ab_neg_zeta_small_",
            # "ab_neg_zeta_big_",
            # "ab_pos_zeta_big_",
            "ab_pos_zeta_neg_",
            "ab_neg_zeta_neg_x",
            "ab_neg_zeta_neg_big_"
            ]
titles = [L"$(\zeta, \; \alpha, \; \beta) = (0.1,\; 1.0,\; 5.0)$",
          L"$(\zeta, \; \alpha, \; \beta) = (0.1,\; -1.0,\; 5.0)$",
          L"$(\zeta, \; \alpha, \; \beta) = (0.1, \;1.0, \;-5.0)$",
          L"$(\zeta, \; \alpha, \; \beta) = (0.1,\; -1.0,\; -5.0)$",
        #   L"$(\zeta, \; \alpha, \; \beta) = (3., \;-1.0, \;-5.0)$",
        #   L"$(\zeta, \; \alpha, \; \beta) = (3., \;1.0, \;5.0)$",
          L"$(\zeta, \; \alpha, \; \beta) = (-0.1, \;1.0, \;5.0)$",
          L"$(\zeta, \; \alpha, \; \beta) = (-0.1, \;-1.0, \;-5.0)$",
          L"$(\zeta, \; \alpha, \; \beta) = (-3., \;-1.0, \;-5.0)$"
]
# for i in eachindex(p)
#     get_trajectories(p[i], filenames[i], titles[i])
# end

# time dependencies
function get_time_dependencies(p::Params, filename, title)
    init_vals = [L"(x_0,\;v_0)=(0.5, \;-2.0)", L"(x_0,\;v_0)=(-1., \;2.0)",
     L"(x_0,\;v_0)=(-1.0, \;0.5)", L"(x_0,\;v_0)=(0.5, \;-1.7)",
     L"(x_0,\;v_0)=(0.0, \;0.1)", L"(x_0,\;v_0)=(0.0, \;-0.5)"]

    # ============================================================
    fig = Figure(size = (1000, 500))
    ax = Axis(fig[1, 1], xlabel=L"t", ylabel=L"$x$", xlabelsize = 30,# limits = ((-0.1, 10.2), (-1.2, 1.2)),
    ylabelsize = 30, title=title, titlesize = 30, xticklabelsize = 20, yticklabelsize = 20)
    ax2 = Axis(fig[1, 2], xlabel=L"t", ylabel=L"$v$", xlabelsize = 30,# limits = ((-0.1, 10.2), (-2.5, 2.7)),
    ylabelsize = 30, title=title, titlesize = 30, xticklabelsize = 20, yticklabelsize = 20)
    
    # ============================================================
    data_dir = "/home/marta/Documents/studia/dynamical_modelling/projekt/build/data/"
    begin_file = "duff_$filename"

    files = filter(f -> occursin(begin_file, f) && endswith(f, ".txt"), readdir(data_dir; join=true))

    println("Found $(length(files)) data files:")
    foreach(println, files)

    # cm = cgrad(:RdYlBu_11, 6)
    # cm = cgrad(:RdBu_6, 6)
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
    # display(fig)
    save("projekt/graphics/time$filename.pdf", fig)
    return fig
end

# for i in eachindex(p)
#     get_time_dependencies(p[i], filenames[i], titles[i])
# end

##
using CairoMakie, DelimitedFiles, FilePathsBase

function plot_bifur(out_prefix)
    poincare_file = "/home/marta/Documents/studia/dynamical_modelling/projekt/build/data/$(out_prefix)_bifur_poincare.txt"
    amp_file = "/home/marta/Documents/studia/dynamical_modelling/projekt/build/data/$(out_prefix)_bifur_amp.txt"

    # load poincare
    data_p = readdlm(poincare_file, '\t', comments=true)
    gammas_p = data_p[:,1]
    xs_p = data_p[:,3]

    fig = Figure(resolution=(1000,600))
    ax = Axis(fig[1,1], xlabel="gamma", ylabel="x (Poincare)",
              title="Poincare map (stroboscopic samples)")
    scatter!(ax, gammas_p, xs_p, markersize=5.5, color=:black)
    display(fig)

    # # now amplitude
    # data_a = readdlm(amp_file, '\t', comments=true)
    # gam = data_a[:,1]
    # amp = data_a[:,2]

    # fig2 = Figure(resolution=(1000,400))
    # ax2 = Axis(fig2[1,1], xlabel="gamma", ylabel="avg period max x",
    #            title="Bifurcation: amplitude vs gamma")
    # lines!(ax2, gam, amp, linewidth=2)
    # scatter!(ax2, gam, amp, markersize=3)
    # display(fig2)
end

plot_bifur("duffing")  # dostosuj prefiks do tego co zapisałeś
