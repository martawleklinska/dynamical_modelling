using CairoMakie
using DelimitedFiles
using FilePathsBase

## ============================================================
α = 3.0
β = 1.
γ = 1.
δ = 2.

function get_arrowed_plot_lotka_volterry()
    
    xs = range(-3, 5, length=20)
    ys = range(-2, 4., length=20)
    
    U = [α*x - β*x*y for x in xs, y in ys]
    V = [-γ*y + δ*x*y for x in xs, y in ys]
    
    mag = sqrt.(U.^2 .+ V.^2).*3
    U ./= mag
    V ./= mag
    
    # ============================================================
    fig = Figure(resolution=(1000, 650))
    ax = Axis(fig[1, 1], xlabel=L"$x$", ylabel=L"$y$", xlabelsize = 30,
    ylabelsize = 30, titlesize = 30, xticklabelsize = 25, yticklabelsize = 25,
    title="α=$α, β=$β, γ=$γ, δ=$δ", limits = ((-3, 5.), (-2, 4.)))

    arrows!(ax, xs, ys, U, V, arrowsize=8, linecolor=:gray, linewidth=1, alpha = 0.6)
    # x = 0 
    lines!(ax, zeros(length(ys)), ys, color = :blue, linewidth = 8., alpha = 0.3, label = L"x=0")
    #  y = 0 
    lines!(ax, xs, zeros(length(xs)), color = :red, alpha = 0.3, linewidth = 8., label = L"y=0")

    #  y = α/β
    lines!(ax, xs, fill(α/β, length(xs)), color = :blue, alpha = 0.3, linewidth = 8., label = L"y=\alpha/\beta")
    # x = γ/δ  
    lines!(ax, fill(γ/δ, length(ys)), ys, color = :red, alpha = 0.3, linewidth = 8., label = L"x=\gamma/\delta")

    scatter!(ax, [γ/δ], [α/β], color = :black, markersize = 18, marker = :star5)#, label = L"\text{punkt staly} ($\gamma/\delta, \alpha/\beta) ")
    scatter!(ax, [0.0], [0.0], color = :black, markersize = 18, marker = :star5)#, label = L"\text{punkt staly} ($0.0, 0.0$)")
    # ============================================================
    data_dir = "/home/marta/Documents/studia/dynamical_modelling/lab02/build/data/"
    files = filter(f -> occursin(r"y0", f) && endswith(f, ".txt"), readdir(data_dir; join=true))

    println("Found $(length(files)) data files:")
    foreach(println, files)

    for file in files
        data = readdlm(file, '\t', skipstart=1)
        if size(data, 2) >= 3
            x_vals = data[:, 2]  
            y_vals = data[:, 3]  
            lines!(ax, x_vals, y_vals, linewidth=2, label=basename(file))
        else
            @warn "cos jest nie tak z $(file) "
        end
    end

    Legend(fig[1,2], ax, labelsize = 20)

    fig
    # save("lotka_volterra.pdf", fig)
end

# get_arrowed_plot_lotka_volterry()


## zadanie \dot{x}=kx(1-x)
function ex1(k)
    fig = Figure(resolution = (1000, 600))
    ax = Axis(fig[1, 1], xlabel=L"$t$", ylabel=L"$x$", xlabelsize = 30,
    ylabelsize = 30, titlesize = 30, xticklabelsize = 25, yticklabelsize = 25,
    title=L"\dot{x} = kx(1-x), \; k=-1", limits = ((0.0, 2.0), (-2., 2.)))
    data_dir = "/home/marta/Documents/studia/dynamical_modelling/lab02/build/data/"
    files = filter(f -> occursin(r"init_k-1", f) && endswith(f, ".txt"), readdir(data_dir; join=true))

    println("Found $(length(files)) data files:")
    foreach(println, files)

    for file in files
        data = readdlm(file, '\t', skipstart=1)
        if size(data, 1) >= 2
            x_vals = data[:, 1]  
            y_vals = data[:, 2]  
            lines!(ax, x_vals, y_vals, linewidth=2, label=basename(file))
        else
            @warn "cos jest nie tak z $(file) "
        end
    end
    Legend(fig[1,2],ax, position = :rb)
    save("ex1_k-1.pdf", fig)
end
# ex1(-1)

## zadanie 1d
function get_analytical_func(x0, t)
    k = 1 
    vals = Float64[]
    if x0<0.5
        for (i, tim) in enumerate(t)
            push!(vals, x0 * exp(k * tim))
        end
    else
        for (i, tim) in enumerate(t)
            push!(vals, (1+(x0-1) * exp(k * tim)))
        end
    end
        return vals
end

function ex1d()
    fig = Figure(resolution = (1000, 600))
    ax = Axis(fig[1, 1], xlabel=L"$t$", ylabel=L"$x$", xlabelsize = 30,
    ylabelsize = 30, titlesize = 30, xticklabelsize = 25, yticklabelsize = 25,
    title=L"\dot{x} = kx(1-x), \; k=1", limits = ((0.0, 2.0), (-2., 2.)))
    data_dir = "/home/marta/Documents/studia/dynamical_modelling/lab02/build/data/"
    files = filter(f -> occursin(r"init_k1d", f) && endswith(f, ".txt"), readdir(data_dir; join=true))

    println("Found $(length(files)) data files:")
    foreach(println, files)
    init_vals = [-0.1, 0.1, 0.9, 1.1]
    for (i, file) in enumerate(files)
        data = readdlm(file, '\t', skipstart=1)
        if size(data, 1) >= 2
            label = init_vals[i]
            x_vals = data[:, 1]  
            y_vals = data[:, 2]  
            y_anal = get_analytical_func(label, x_vals)
            lines!(ax, x_vals, y_vals, linewidth=2, label=basename(file))
            lines!(ax, x_vals, y_anal, linewidth = 2, linestyle = :dash, label = "x0=$label")
        else
            @warn "cos jest nie tak z $(file) "
        end
    end
    Legend(fig[1,2],ax, position = :rb)
    save("lab02/ex1d_k1_anal.pdf", fig)
    # fig
end

ex1d()

## ex1e
function get_potential(k)
    xs = range(-1, 2, length=40)
    V = [k/3 * x^3 - k/2 * x^2 for x in xs]

    fig = Figure(resolution = (1000, 600))
    ax = Axis(fig[1, 1], xlabel=L"$x$", ylabel=L"$V(x)$", xlabelsize = 30,
    ylabelsize = 30, titlesize = 30, xticklabelsize = 25, yticklabelsize = 25,
    title=L"\text{ćwiczenie 1e}, $k=1$")
    lines!(ax, xs, V)
    save("lab02/ex1e_potentialk1.pdf", fig)
    # fig
end
get_potential(1)


## ex2

function ex2()
    init_vals = [
                L"(x_0, y_0)=(-2.0, 0.5)",
                L"(x_0, y_0)=(-3.5, 1.0)",
                L"(x_0, y_0)=(0.0, -0.5)",
                L"(x_0, y_0)=(0.1, 0.1)", 
                L"(x_0, y_0)=(0.5, 0.2)",
                L"(x_0, y_0)=(2.0, 1.5)",
                ]
    xs = range(-30, 30, length=40)
    ys = range(-30., 30., length=40)
    
    U = [x * (y - 1) for x in xs, y in ys]
    V = [3 * x - 2 * y + x ^ 2- 2 * y^2 for x in xs, y in ys]
    
    mag = sqrt.(U.^2 .+ V.^2).*3
    U ./= mag
    V ./= mag
    
    # ============================================================
    fig = Figure(resolution=(1000, 650))
    ax = Axis(fig[1, 1], xlabel=L"$x$", ylabel=L"$y$", xlabelsize = 30,
    ylabelsize = 30, titlesize = 30, xticklabelsize = 25, yticklabelsize = 25,
    title=L"\text{ćwiczenie 2}", limits = ((-10., 10.), (-10., 10.)))

    arrows!(ax, xs, ys, U, V, arrowsize=10,lengthscale = 1.5, linecolor=:gray, linewidth=8, alpha = 0.6)
    # x = 0 
    # lines!(ax, zeros(length(ys)), ys, color = :blue, linewidth = 8., alpha = 0.3, label = L"x=0")
    # #  y = 0 
    # lines!(ax, xs, zeros(length(xs)), color = :red, alpha = 0.3, linewidth = 8., label = L"y=0")

    # #  y = α/β
    # lines!(ax, xs, fill(1, length(xs)), color = :blue, alpha = 0.3, linewidth = 8., label = L"y=\alpha/\beta")
    # # x = γ/δ  
    # lines!(ax, fill(-3, length(ys)), ys, color = :red, alpha = 0.3, linewidth = 8., label = L"x=\gamma/\delta")

    
    # ============================================================
    # data_dir = "/home/marta/Documents/studia/dynamical_modelling/lab02/build/data/"
    # files = filter(f -> occursin(r"sys2_traj", f) && endswith(f, ".txt"), readdir(data_dir; join=true))

    # println("Found $(length(files)) data files:")
    # foreach(println, files)

    # for (i, file) in enumerate(files)
    #     data = readdlm(file, '\t', skipstart=1)
    #     label = init_vals[i]
    #     if size(data, 2) >= 3
    #         x_vals = data[:, 2]  
    #         y_vals = data[:, 3]  
    #         lines!(ax, x_vals, y_vals, linewidth=4, label=label, alpha = 0.7)
    #     else
    #         @warn "cos jest nie tak z $(file) "
    #     end
    # end
    Δ(x) = 1 + 2x^2 + 6x
    y_plus(x)  = (-1 + sqrt(Δ(x)))/2
    y_minus(x) = (-1 - sqrt(Δ(x)))/2

    ys_plus = [Δ(x) >= 0 ? y_plus(x) : NaN for x in xs]
    ys_minus= [Δ(x) >= 0 ? y_minus(x) : NaN for x in xs]

    lines!(ax, xs, ys_plus, linewidth = 10, alpha = 0.3, color = :red, label=L"\text{1) }\dot{y} = 0")
    lines!(ax, xs, ys_minus, linewidth = 10, alpha = 0.3, color = :orange, label=L"\text{2) }\dot{y} = 0")
    vlines!(ax, [0], linewidth = 10, alpha = 0.3, color = :blue, label=L"$x=0$")
    hlines!(ax, [1], linewidth = 10, alpha = 0.3, color = :blue, label=L"$y=1$")

    scatter!(ax, [0.0], [0.0], color = :black, markersize = 25, marker = :star5)#, label = L"\text{punkt staly} ($\gamma/\delta, \alpha/\beta) ")
    scatter!(ax, [-4.0], [1.0], color = :black, markersize = 25, marker = :star5)#, label = L"\text{punkt staly} ($0.0, 0.0$)")
    scatter!(ax, [1.0], [1.0], color = :black, markersize = 25, marker = :star5)#, label = L"\text{punkt staly} ($0.0, 0.0$)")
    scatter!(ax, [0.0], [-1.0], color = :black, markersize = 25, marker = :star5)#, label = L"\text{punkt staly} ($0.0, 0.0$)")
    Legend(fig[1,2], ax, labelsize = 24)

    # fig
    save("lab02/ex2_izokliny.pdf", fig)
end
ex2()

## ex2c
function get_xt_yt()
    fig = Figure(resolution=(1000, 650))
    ax = Axis(fig[1, 1], xlabel=L"$t$", ylabel=L"$y$", xlabelsize = 30,
    ylabelsize = 30, titlesize = 30, xticklabelsize = 28, yticklabelsize = 28,
    title=L"\text{ćwiczenie 2c}", limits = ((-0.03, 2.5), (-5., 5.)))

    init_vals = [
                L"(x_0, y_0)=(-2.0, 0.5)",
                L"(x_0, y_0)=(-3.5, 1.0)",
                L"(x_0, y_0)=(0.0, -0.5)",
                L"(x_0, y_0)=(0.1, 0.1)", 
                L"(x_0, y_0)=(0.5, 0.2)",
                L"(x_0, y_0)=(2.0, 1.5)",
                ]
    xs = range(-30, 30, length=40)
    ys = range(-30., 30., length=40)
    data_dir = "/home/marta/Documents/studia/dynamical_modelling/lab02/build/data/"
    files = filter(f -> occursin(r"sys2_traj", f) && endswith(f, ".txt"), readdir(data_dir; join=true))

    println("Found $(length(files)) data files:")
    foreach(println, files)

    for (i, file) in enumerate(files)
        data = readdlm(file, '\t', skipstart=1)
        label = init_vals[i]
        if size(data, 2) >= 3
            t_vals = data[:, 1]
            x_vals = data[:, 2]  
            y_vals = data[:, 3]  
            lines!(ax, t_vals, y_vals, linewidth=4, label=label, alpha = 0.7)
        else
            @warn "cos jest nie tak z $(file) "
        end
    end
    axislegend(ax, position = :rt, labelsize = 25)
    # fig
    save("lab02/ex2c_yt.pdf", fig)
end
get_xt_yt()