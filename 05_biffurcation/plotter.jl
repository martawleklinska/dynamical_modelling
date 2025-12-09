using CairoMakie
using DelimitedFiles

function get_wykres_xm()
# ------------------- G(x,m) = m - |x| -------------------
ms = range(-1, 2, length=300)
    x_star_pos = [m > 0 ? m : NaN for m in ms]
    x_star_neg = [m > 0 ? -m : NaN for m in ms]

    fig = Figure(size=(900,400))

    ax1 = Axis(fig[1,1], xlabel=L"m", ylabel=L"x^*", xlabelsize = 30, ylabelsize = 30, title=L"x^*(m) \text{ dla } \dot{x}=m-|x|", titlesize = 30, xticklabelsize = 20, yticklabelsize = 20)
    lines!(ax1, ms, x_star_pos, linewidth=3, label="stabilny")
    lines!(ax1, ms, x_star_neg, linewidth=3, linestyle=:dash, label="niestabilny")
    axislegend(ax1, position=:rc)
    vlines!(ax1, [0], color=:gray, linestyle=:dash)
    xlims!(ax1, (-1,2)); ylims!(ax1, (-2,2))

    save("lab05/graphics/ex1a_xm.pdf", fig)
    display(fig)
end

get_wykres_xm()

## ex1_numerical 
function get_numerical_solutions(m)
    fig = Figure(size = (600, 400))
    ax = Axis(fig[1,1], xlabel = L"t", ylabel = L"$x$", title = L"ex1: $m = 0.0$", 
    xlabelsize = 30, ylabelsize = 30, titlesize = 30, xticklabelsize = 20, 
    yticklabelsize = 20, limits = ((0,2), (-1,1)))
    data_dir = "/home/marta/Documents/studia/dynamical_modelling/lab05/build/data/"
    begin_file = r"ex1_m0.0"

    files = filter(f -> occursin(begin_file, f) && endswith(f, ".txt"), readdir(data_dir; join=true))

    println("Found $(length(files)) data files:")
    foreach(println, files)

    for (i, file) in enumerate(files)
        data = readdlm(file, '\t', skipstart=1)
        # label = init_vals[i]
        if size(data, 2) >= 2
            t_vals = data[:, 1]
            x_vals = data[:, 2] 
            lines!(ax, t_vals, x_vals, linewidth=2,  alpha = 0.7)
        else
            @warn "cos jest nie tak z $(file) "
        end
    end
    lines!(ax, [0,5], [m, m], linewidth = 3, linestyle = :dash)
    lines!(ax, [0,5], [-m, -m], linewidth = 3, linestyle = :dash)
    # display(fig);
    save("lab05/graphics/ex1c_m$m.pdf", fig)
end
get_numerical_solutions(0.0)
## ex2_stabilność
function get_wykres_xp()
        # # ------------------- H(x,p) = (x-1)(x^2+2x-p) -------------------
    fig = Figure(size=(900,400))
    ps_stab = range(-4, 3, length=300)
    x1 = ones(length(ps_stab))
    ps_niestab = range(3, 6, length=300)
    x11 = ones(length(ps_niestab))

    ps_other= range(-1, 6, length = 100)
    ps_other1 = range(-1, 3, length = 100)
    ps_other2 = range(3, 6, length = 100)
    x2 = [-1 + sqrt(1+p) for p in ps_other]
    x31 = [-1 - sqrt(1+p) for p in ps_other1]
    x32 = [-1 - sqrt(1+p) for p in ps_other2]
    # x2[p -> 1+p < 0] .= NaN
    # x3[p -> 1+p < 0] .= NaN

    ax2 = Axis(fig[1,1], xlabel=L"p", titlesize = 30, ylabel=L"x^*", title=L"x^*(p) \text{ dla } \dot{x}=(x-1)(x^2+2x-p)",
    xlabelsize = 30, ylabelsize = 30, )
    lines!(ax2, ps_stab, x1, color=:orange, linewidth=3, label=L"x^*_1=1", linestyle = :dash)
    lines!(ax2, ps_niestab, x11, color=:orange, linewidth=3, label=L"x^*_1=1")
    lines!(ax2, ps_other, x2, color=:green, linewidth=3, linestyle = :dash, label=L"x^*_2=-1+\sqrt{1+p}")
    lines!(ax2, ps_other1, x31, color=:red, linewidth=3, label=L"x^*_3=-1-\sqrt{1+p}")
    lines!(ax2, ps_other2, x32, color=:red, linewidth=3, linestyle=:dash, label=L"x^*_3=-1-\sqrt{1+p}")
    
    vlines!(ax2, [-1], color=:gray, linestyle=:dash)
    axislegend(ax2, position=:lb, labelsize = 20)
    # display(fig)
    save("lab05/graphics/ex2a_stab.pdf", fig)
end

get_wykres_xp()

## ## ex1_numerical 
function get_numerical_solutions_ex2(p)
    fig = Figure(size = (800, 400))
    ax = Axis(fig[1,1], xlabel = L"t", ylabel = L"$x$", title = L"ex2: $p = 5$", 
    xlabelsize = 25, ylabelsize = 25, titlesize = 25, xticklabelsize = 20, 
    yticklabelsize = 20, limits = ((0,2), (-4,3)))
    ax2 = Axis(fig[1,2], xlabel = L"t", ylabel = L"$x$", title = L"ex2: $p = 5$", 
    xlabelsize = 25, ylabelsize = 25, titlesize = 25, xticklabelsize = 20, 
    yticklabelsize = 20, limits = ((0,1), (0.5,2.)))
    
    data_dir = "/home/marta/Documents/studia/dynamical_modelling/lab05/build/data/"
    begin_file = r"ex2"

    files = filter(f -> occursin(begin_file, f) && endswith(f, ".txt"), readdir(data_dir; join=true))

    println("Found $(length(files)) data files:")
    foreach(println, files)

    for (i, file) in enumerate(files)
        data = readdlm(file, '\t', skipstart=1)
        # label = init_vals[i]
        if size(data, 2) >= 2
            t_vals = data[:, 1]
            x_vals = data[:, 2] 
            lines!(ax, t_vals, x_vals, linewidth=2,  alpha = 0.7)
            lines!(ax2, t_vals, x_vals, linewidth = 2, alpha = 0.7)
        else
            @warn "cos jest nie tak z $(file) "
        end
    end
    # lines!(ax, [0,5], [m, p], linewidth = 3, linestyle = :dash)
    # lines!(ax, [0,5], [-m, -m], linewidth = 3, linestyle = :dash)
    # display(fig);
    save("lab05/graphics/ex2b_p$p.pdf", fig)
end
get_numerical_solutions_ex2(5)