using CairoMakie
using DelimitedFiles

function get_phase_graph()
    xs = range(-3, 3., length=30)
    ys = range(-3, 3., length=30)
    
    
    X = [-2*y-x*((x^2+y^2)^2-4*(x^2+y^2)+1) for x in xs, y in ys]
    Y = [2*x-y*((x^2+y^2)^2-4*(x^2+y^2)+1) for x in xs, y in ys]
 
    
    mag = sqrt.(X.^2 .+ Y.^2).*30
    X ./= mag
    Y ./= mag
    
    # ============================================================
    fig = Figure()
    ax = Axis(fig[1, 1], xlabel=L"$\dot{x}$", ylabel=L"$\dot{y}$", xlabelsize = 25, limits = ((-1,1), (-1,1)),
    ylabelsize = 25, title=L"\text{ćwiczenie 1b}", titlesize = 25, xticklabelsize = 20, yticklabelsize = 20)
    arrows!(ax, xs, ys, X, Y, arrowsize=8,lengthscale = 1.5, linecolor=:gray, linewidth=8, alpha = 0.6)
    
    # display(fig)
    save("lab04/graphics/ex1b_male.pdf", fig)
end
get_phase_graph()

## trajectories
function get_trajs()
    # 
    init_vals = [
                L"x_0=0.2",
                L"x_0=0.4",
                L"x_0=0.5175",
                L"x_0=0.5180", 
                L"x_0=1.0",
                L"x_0=1.8,",
                L"x_0=2.5"
                ]
    xs = range(-3, 3., length=30)
    ys = range(-3, 3., length=30)
    
    
    X = [-2*y-x*((x^2+y^2)^2-4*(x^2+y^2)+1) for x in xs, y in ys]
    Y = [2*x-y*((x^2+y^2)^2-4*(x^2+y^2)+1) for x in xs, y in ys]
 
    
    mag = sqrt.(X.^2 .+ Y.^2).*10
    X ./= mag
    Y ./= mag

    # ============================================================
    fig = Figure(size = (1000, 600))
    ax = Axis(fig[1, 1], xlabel=L"$x$", ylabel=L"$y$", xlabelsize = 30,
    ylabelsize = 30, titlesize = 30, xticklabelsize = 25, yticklabelsize = 25,
    title=L"\text{ćwiczenie 1c}")

    arrows!(ax, xs, ys, X, Y, arrowsize=10,lengthscale = 1.5, linecolor=:gray, linewidth=8, alpha = 0.6)
    
    # ============================================================
    data_dir = "/home/marta/Documents/studia/dynamical_modelling/lab04/build/data/"
    begin_file = r"sys3_traj_x0_"

    files = filter(f -> occursin(begin_file, f) && endswith(f, ".txt"), readdir(data_dir; join=true))

    println("Found $(length(files)) data files:")
    foreach(println, files)

    for (i, file) in enumerate(files)
        data = readdlm(file, '\t', skipstart=1)
        label = init_vals[i]
        if size(data, 2) >= 3
            x_vals = data[:, 2]  
            y_vals = data[:, 3]  
            lines!(ax, x_vals, y_vals, linewidth=4, label = label, alpha = 0.7)
        else
            @warn "cos jest nie tak z $(file) "
        end
    end
    Legend(fig[1,2], ax, labelsize = 20)
    # display(fig)
    save("lab04/graphics/ex1c.pdf", fig)
end
get_trajs()



## zależność czasowa
function get_xt_yt()
    fig = Figure(size = (1500, 600))

    ax = Axis(fig[1, 1], xlabel=L"$t$", ylabel=L"$x$", xlabelsize = 30,
    ylabelsize = 30, titlesize = 30, xticklabelsize = 25, yticklabelsize = 25,
    title=L"\text{ćwiczenie 1c}")

    init_vals = [
                L"x_0=0.2",
                L"x_0=0.4",
                L"x_0=0.5175",
                L"x_0=0.5180", 
                L"x_0=1.0",
                L"x_0=1.8,",
                L"x_0=2.5"
                ]

    data_dir = "/home/marta/Documents/studia/dynamical_modelling/lab04/build/data/"
    begin_file = r"sys3_traj_x0_"

    files = filter(f -> occursin(begin_file, f) && endswith(f, ".txt"), readdir(data_dir; join=true))

    println("Found $(length(files)) data files:")
    foreach(println, files)

    for (i, file) in enumerate(files)
        data = readdlm(file, '\t', skipstart=1)
        label = init_vals[i]
        if size(data, 2) >= 3
            t_vals = data[:, 1]
            x_vals = data[:, 2] # theta 
            y_vals = data[:, 3]  # omega
            lines!(ax, t_vals, x_vals, linewidth=4, label=label, alpha = 0.7)
        else
            @warn "cos jest nie tak z $(file) "
        end
    end
    ax = Axis(fig[1,2], xlabel=L"$t$", ylabel=L"$y$", xlabelsize = 30,
    ylabelsize = 30, titlesize = 30, xticklabelsize = 28, yticklabelsize = 28,)
    for (i, file) in enumerate(files)
        data = readdlm(file, '\t', skipstart=1)
        label = init_vals[i]
        if size(data, 2) >= 3
            t_vals = data[:, 1]
            x_vals = data[:, 2] # theta 
            y_vals = data[:, 3]  # omega
            lines!(ax, t_vals, y_vals, linewidth=4, label=label, alpha = 0.7)
        else
            @warn "cos jest nie tak z $(file) "
        end
    end
    Legend(fig[1,3], ax, position = :rt, L"\text{model 3}", labelsize = 30, titlesize = 30)
    
    # filenames = ["ex1_5_time_dep_model1", "ex1_5_time_dep_model2", "ex1_5_time_dep_model3"]
    # filename = filenames[model]
    # display(fig)
    save("lab04/graphics/ex1c_time.pdf", fig)
end
get_xt_yt()

## radius eqn
function get_radius_theta()
    fig = Figure(size = (1500, 600))
    # titlelayout = GridLayout(fig[0, 1:2], tellwidth = false)
    # Label(titlelayout[1, 1], "ćwiczenie 1e", textsize = 20)
    # rowgap!(titlelayout, 0)
    ax = Axis(fig[1, 1], xlabel=L"$t$", ylabel=L"$r$", xlabelsize = 30,
    ylabelsize = 30, titlesize = 30, xticklabelsize = 25, yticklabelsize = 25,)

    init_vals = [
                L"r_0=0.2",
                L"r_0=0.4",
                L"r_0=0.5175",
                L"r_0=0.5180", 
                L"r_0=1.0",
                L"r_0=1.8,",
                L"r_0=2.5"
                ]

    data_dir = "/home/marta/Documents/studia/dynamical_modelling/lab04/build/data/"
    begin_file = r"rad__x0_"

    files = filter(f -> occursin(begin_file, f) && endswith(f, ".txt"), readdir(data_dir; join=true))

    println("Found $(length(files)) data files:")
    foreach(println, files)

    for (i, file) in enumerate(files)
        data = readdlm(file, '\t', skipstart=1)
        label = init_vals[i]
        if size(data, 2) >= 2
            t_vals = data[:, 1] # time
            x_vals = data[:, 2] # radius
            lines!(ax, t_vals, x_vals, linewidth=4, label=label, alpha = 0.7)
        else
            @warn "cos jest nie tak z $(file) "
        end
    end
    ax2 = Axis(fig[1,2], xlabel=L"$t$", ylabel=L"$\vartheta$", xlabelsize = 30,
    ylabelsize = 30, titlesize = 30, xticklabelsize = 28, yticklabelsize = 28)
    ts = LinRange(0, 10, 500)
    lines!(ax2, ts, 2 * ts)
    Legend(fig[2,1:2], orientation = :horizontal, ax, position = :rt, L"\text{warunki początkowe}", labelsize = 30, titlesize = 30)
    # display(fig)
    save("lab04/graphics/ex1e_time.pdf", fig)
end
get_radius_theta()

## biegunowe
function polar()
    fig = Figure(size = (1500, 600))

    ax = Axis(fig[1, 1], xlabel=L"$r$", ylabel=L"$\vartheta$", xlabelsize = 30,
    ylabelsize = 30, titlesize = 30, xticklabelsize = 25, yticklabelsize = 25,)

    init_vals = [
                L"r_0=0.2",
                L"r_0=0.4",
                L"r_0=0.5175",
                L"r_0=0.5180", 
                L"r_0=1.0",
                L"r_0=1.8,",
                L"r_0=2.5"
                ]

    data_dir = "/home/marta/Documents/studia/dynamical_modelling/lab04/build/data/"
    begin_file = r"rad__x0_"

    files = filter(f -> occursin(begin_file, f) && endswith(f, ".txt"), readdir(data_dir; join=true))

    println("Found $(length(files)) data files:")
    foreach(println, files)

    for (i, file) in enumerate(files)
        data = readdlm(file, '\t', skipstart=1)
        label = init_vals[i]
        if size(data, 2) >= 2
            t_vals = data[:, 1] # time
            x_vals = data[:, 2] # radius
            lines!(ax, x_vals, 2*t_vals, linewidth=4, label=label, alpha = 0.7)
        else
            @warn "cos jest nie tak z $(file) "
        end
    end

    Legend(fig[2,1], orientation = :horizontal, ax, 
    position = :rt, L"\text{warunki początkowe}", 
    labelsize = 30, titlesize = 30, framevisible = false)
    # display(fig)
    save("lab04/graphics/ex1e_polar.pdf", fig)
end
polar()

## ex2 
function get_phase_graph_ex2()
    xs = range(-3, 3., length=30)
    ys = range(-3, 3., length=30)
    
    
    X = [u for u in xs, y in ys]
    Y = [2 * (1 - 1/3 * u) * u - x for u in xs, x in ys]
 
    
    mag = sqrt.(X.^2 .+ Y.^2).*30
    X ./= mag
    Y ./= mag
    
    # ============================================================
    fig = Figure()
    ax = Axis(fig[1, 1], xlabel=L"$\dot{u}$", ylabel=L"$\dot{x}$", xlabelsize = 25, limits = ((-1,1), (-1,1)),
    ylabelsize = 25, title=L"\text{ćwiczenie 2:portret fazowy}", titlesize = 25, xticklabelsize = 20, yticklabelsize = 20)
    arrows!(ax, xs, ys, X, Y, arrowsize=8,lengthscale = 1.5, linecolor=:gray, linewidth=8, alpha = 0.6)
    
    # display(fig)
    save("lab04/graphics/ex2_phase.pdf", fig)
end
get_phase_graph_ex2()

##
function ex2_trajs()
    xs = range(-5, 5, length=40)
    us = range(-5, 5, length=40)

    X = [u for x in xs, u in us]                            # dx/dt = u
    U= [- (2/3) * u^3 + 2 * u - x for x in xs, u in us]   # du/dt = -2/3 u^3 + 2u - x

    mag = sqrt.(X.^2 .+ U.^2)
    mag[mag .== 0] .= 1
    X ./= mag
    U ./= mag

    fig = Figure(size = (1000, 600))
    ax = Axis(fig[1, 1],
        xlabel=L"$x$", ylabel=L"$u$",xlabelsize = 25, 
        ylabelsize = 25, title=L"\text{ćwiczenie 2: rozwiązania}", titlesize = 25, xticklabelsize = 20, yticklabelsize = 20,
        limits = ((-5,5), (-5,5)))

    arrows!(ax, xs, us, X, U, arrowsize=10, lengthscale=0.5,
        linecolor=:gray, linewidth=4, alpha=0.6)

    # ============================================================
    data_dir = "/home/marta/Documents/studia/dynamical_modelling/lab04/build/data/"
    begin_file = r"ex2"

    files = filter(f -> occursin(begin_file, f) && endswith(f, ".txt"), readdir(data_dir; join=true))

    println("Found $(length(files)) data files:")
    foreach(println, files)

    for (i, file) in enumerate(files)
        data = readdlm(file, '\t', skipstart=1)
        if size(data, 2) >= 3
            x_vals = data[:, 2]  
            u_vals = data[:, 3]  
            lines!(ax, x_vals, u_vals, linewidth=4, alpha = 0.7)
        else
            @warn "cos jest nie tak z $(file) "
        end
    end

    # display(fig)
    save("lab04/graphics/ex2_trajs.pdf", fig)
end

ex2_trajs()

##
function get_xtut()
    fig = Figure(size = (1500, 600))

    ax = Axis(fig[1, 1], xlabel=L"$t$", ylabel=L"$x$", xlabelsize = 30,
    ylabelsize = 30, titlesize = 30, xticklabelsize = 25, yticklabelsize = 25)

    init_vals = [
                L"(x_0,\; u_0)=(0.0,\; 4.0)",
                L"(x_0,\; u_0)=(0.0,\; -4.0)", 
                L"(x_0,\; u_0)=(-0.7,\; -0.5)",
                L"(x_0,\; u_0)=(-1.5,\; -3)",
                L"(x_0,\; u_0)=(1.,\; -3.0)", 
                L"(x_0,\; u_0)=(0.5,\; -3.)"
                ]

    data_dir = "/home/marta/Documents/studia/dynamical_modelling/lab04/build/data/"
    begin_file = r"ex2"

    files = filter(f -> occursin(begin_file, f) && endswith(f, ".txt"), readdir(data_dir; join=true))

    println("Found $(length(files)) data files:")
    foreach(println, files)

    for (i, file) in enumerate(files)
        data = readdlm(file, '\t', skipstart=1)
        label = init_vals[i]
        if size(data, 2) >= 3
            t_vals = data[:, 1]
            x_vals = data[:, 2] # theta 
            y_vals = data[:, 3]  # omega
            lines!(ax, t_vals, x_vals, linewidth=4, label = label, alpha = 0.7)
        else
            @warn "cos jest nie tak z $(file) "
        end
    end
    ax = Axis(fig[1,2], xlabel=L"$t$", ylabel=L"$u$", xlabelsize = 30,
    ylabelsize = 30, titlesize = 30, xticklabelsize = 28, yticklabelsize = 28,)
    for (i, file) in enumerate(files)
        data = readdlm(file, '\t', skipstart=1)
        label = init_vals[i]
        if size(data, 2) >= 3
            t_vals = data[:, 1]
            x_vals = data[:, 2] # theta 
            y_vals = data[:, 3]  # omega
            lines!(ax, t_vals, y_vals, linewidth=4, label = label, alpha = 0.7)
        else
            @warn "cos jest nie tak z $(file) "
        end
    end
    Legend(fig[2,1:2], ax, orientation=:horizontal, framevisible=false, L"\text{warunki początkowe}", labelsize = 20, titlesize = 30)
    
    # display(fig)
    save("lab04/graphics/ex2_time.pdf", fig)
end
get_xtut()