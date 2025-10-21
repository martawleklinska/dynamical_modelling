using CairoMakie
using DelimitedFiles
using FilePathsBase

## ============================================================

function portret_fazowy(;model::Int = 1)
    xs = range(-10, 10., length=20)
    ys = range(-5, 5., length=20)
    
    if model == 1
        X = [y for x in xs, y in ys]
        Y = [-x for x in xs, y in ys]
    elseif model == 2
        X = [y for x in xs, y in ys]
        Y = [-x+1/6*x^3 for x in xs, y in ys]
    else
        X = [y for x in xs, y in ys]
        Y = [-sin(x) for x in xs, y in ys]
    end
    
    mag = sqrt.(X.^2 .+ Y.^2).*3
    X ./= mag
    Y ./= mag
    
    title_labels = [L"\text{model 1}", L"\text{model 2}", L"\text{model 3}"]
    
    # ============================================================
    fig = Figure()
    ax = Axis(fig[1, 1], xlabel=L"$\vartheta$ \text{[rad]}", ylabel=L"$\dot{\vartheta}$ \text{[rad/s]}", xlabelsize = 25,
    ylabelsize = 25, title=title_labels[model], titlesize = 25, xticklabelsize = 20, yticklabelsize = 20)
    arrows!(ax, xs, ys, X, Y, arrowsize=3, linecolor=:gray, linewidth=1, alpha = 0.6)
    
    # punkty stałe
    if model == 1
        scatter!(ax, 0.0, 0.0, color = :blue, markersize = 25, marker = :star5)
    elseif model == 2
        scatter!(ax, 0.0, 0.0, color = :blue, markersize = 25, marker = :star5)
        scatter!(ax, -sqrt(6), 0.0, color = :red, markersize = 25, marker = :star5)
        scatter!(ax, sqrt(6), 0.0, color = :red, markersize = 25, marker = :star5)
    else
        scatter!(ax, 0.0, 0.0, color = :blue, markersize = 25, marker = :star5)
        scatter!(ax, π, 0.0, color = :red, markersize = 25, marker = :star5)
        scatter!(ax, -π, 0.0, color = :red, markersize = 25, marker = :star5)
        scatter!(ax, 2*π, 0.0, color = :blue, markersize = 25, marker = :star5)
        scatter!(ax, -2*π, 0.0, color = :blue, markersize = 25, marker = :star5)
        scatter!(ax, -3*π, 0.0, color = :red, markersize = 25, marker = :star5)
        scatter!(ax, 3*π, 0.0, color = :red, markersize = 25, marker = :star5)
    end

    
    filename_labels = ["ex1_first_approx", "ex1_second_approx", "ex1_sine_approx"]
    # display(fig)
    filename_label = filename_labels[model]
    save("lab03/graphics/$filename_label.pdf", fig)
end

portret_fazowy(model = 1)

## ex2

function ex2(;model::Int = 1)
    init_vals = [
                L"(\vartheta_0, \omega_0)=(0.0, 0.5)",
                L"(\vartheta_0, \omega_0)=(0.0, 1.0)",
                L"(\vartheta_0, \omega_0)=(0.0, 1.5)",
                L"(\vartheta_0, \omega_0)=(0.0, 2.0)", 
                L"(\vartheta_0, \omega_0)=(0.0, 2.5)",
                L"(\vartheta_0, \omega_0)=(0.0, 3.0)",
                ]
    xs = range(-10, 10., length=20)
    ys = range(-5, 5., length=20)

    if model == 1
        X = [y for x in xs, y in ys]
        Y = [-x for x in xs, y in ys]
    elseif model == 2
        X = [y for x in xs, y in ys]
        Y = [-x+1/6*x^3 for x in xs, y in ys]
    else
        X = [y for x in xs, y in ys]
        Y = [-sin(x) for x in xs, y in ys]
    end
    
    mag = sqrt.(X.^2 .+ Y.^2).*3
    X ./= mag
    Y ./= mag
    if model == 1
        titlenames = L"\text{model 1}"
    elseif model == 2
        titlenames = L"\text{model 2}"  
    else 
        titlenames = L"\text{model 3}"
    end
    # ============================================================
    fig = Figure(size = (1000, 600))
    ax = Axis(fig[1, 1], xlabel=L"$\vartheta$ \text{[rad]}", ylabel=L"$\dot{\vartheta}$ \text{[rad/s]}", xlabelsize = 30,
    ylabelsize = 30, titlesize = 30, xticklabelsize = 25, yticklabelsize = 25,
    title=titlenames, limits = ((-10., 10.), (-5., 5.)))

    arrows!(ax, xs, ys, X, Y, arrowsize=10,lengthscale = 1.5, linecolor=:gray, linewidth=8, alpha = 0.6)
    
    # ============================================================
    data_dir = "/home/marta/Documents/studia/dynamical_modelling/lab03/build/data/"
    if model == 1
        begin_file = r"pend_model1"
    elseif model == 2
        begin_file = r"pend_model2"
    else
        begin_file = r"pend_model3"
    end
    files = filter(f -> occursin(begin_file, f) && endswith(f, ".txt"), readdir(data_dir; join=true))

    println("Found $(length(files)) data files:")
    foreach(println, files)

    for (i, file) in enumerate(files)
        data = readdlm(file, '\t', skipstart=1)
        label = init_vals[i]
        if size(data, 2) >= 3
            x_vals = data[:, 2]  
            y_vals = data[:, 3]  
            lines!(ax, x_vals, y_vals, linewidth=4, label=label, alpha = 0.7)
        else
            @warn "cos jest nie tak z $(file) "
        end
    end
    Legend(fig[1,2], ax, labelsize = 30)

    # display(fig)
    filenames = ["ex1_rozw_1-approx", "ex1_rozw_2-approx", "ex1_rozw_sin"]
    filename = filenames[model]
    save("lab03/graphics/$filename.pdf", fig)
end
ex2(model = 1)

## zależność czasowa
function get_xt_yt(;model::Int = 1)
    fig = Figure(size = (1500, 600))
    if model == 1
        lims = ((-0.03, 10.), (-3.2, 3.2))
        lims2 = lims
    elseif model == 2
        lims = ((-0.03, 6), (-5, 20.))
        lims2 = lims
    else
        lims = ((-0.03, 10.), (-3., 26.))
        lims2 = ((-0.03, 10.), (-1.7, 3.1))
    end
    ax = Axis(fig[1, 1], xlabel=L"$t$ [s]", ylabel=L"$\vartheta$ [rad]", xlabelsize = 30,
    ylabelsize = 30, titlesize = 30, xticklabelsize = 28, yticklabelsize = 28,
    title=L"$\vartheta(t)$", limits = lims)
    init_vals = [
                L"(\vartheta_0, \omega_0)=(0.0, 0.5)",
                L"(\vartheta_0, \omega_0)=(0.0, 1.0)",
                L"(\vartheta_0, \omega_0)=(0.0, 1.5)",
                L"(\vartheta_0, \omega_0)=(0.0, 2.0)", 
                L"(\vartheta_0, \omega_0)=(0.0, 2.5)",
                L"(\vartheta_0, \omega_0)=(0.0, 3.0)",
                ]

    data_dir = "/home/marta/Documents/studia/dynamical_modelling/lab03/build/data/"
    if model == 1
        begin_file = r"pend_model1"
    elseif model == 2
        begin_file = r"pend_model2"
    else
        begin_file = r"pend_model3"
    end
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
    ax = Axis(fig[1,2], xlabel=L"$t$ \text{[s]}", ylabel=L"$\omega$ \text{[rad/s]}", xlabelsize = 30,
    ylabelsize = 30, titlesize = 30, xticklabelsize = 28, yticklabelsize = 28,
    title=L"$\omega(t)$", limits = lims2)
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
    
    filenames = ["ex1_5_time_dep_model1", "ex1_5_time_dep_model2", "ex1_5_time_dep_model3"]
    filename = filenames[model]
    # display(fig)
    save("lab03/graphics/$filename.pdf", fig)
end
get_xt_yt(model =3)

## energia w czasie
function get_energy_plot(;model::Int = 1)
    fig = Figure(size = (1000, 600))
    lims = ((-0.03, 10.), (-0.05, 500))
    if model == 1
        titlename = L"\text{model 1}"
    elseif model == 2
        titlename = L"\text{model 2}"
    else 
        titlename = L"\text{model 3}"
    end

    ax = Axis(fig[1, 1], xlabel=L"$t$ \text{[s]}", ylabel=L"$E(t)$ \text{[J]}", xlabelsize = 30,
    ylabelsize = 30, titlesize = 30, xticklabelsize = 28, yticklabelsize = 28,
    title=titlename, limits = lims)
    init_vals = [
                L"(\vartheta_0, \omega_0)=(0.0, 0.5)",
                L"(\vartheta_0, \omega_0)=(0.0, 1.0)",
                L"(\vartheta_0, \omega_0)=(0.0, 1.5)",
                L"(\vartheta_0, \omega_0)=(0.0, 2.0)", 
                L"(\vartheta_0, \omega_0)=(0.0, 2.5)",
                L"(\vartheta_0, \omega_0)=(0.0, 3.0)",
                ]


    data_dir = "/home/marta/Documents/studia/dynamical_modelling/lab03/build/data/"
    if model == 1
        begin_file = r"pend_model1"
    elseif model == 2
        begin_file = r"pend_model2"
    else
        begin_file = r"pend_model3"
    end
    files = filter(f -> occursin(begin_file, f) && endswith(f, ".txt"), readdir(data_dir; join=true))

    println("Found $(length(files)) data files:")
    foreach(println, files)

    for (i, file) in enumerate(files)
        data = readdlm(file, '\t', skipstart=1)
        label = init_vals[i]
        if size(data, 2) >= 4
            t_vals = data[:, 1]
            e_vals = data[:, 4] # energy 
            lines!(ax, t_vals, e_vals, linewidth=4, label = label, alpha = 0.7)
        else
            @warn "cos jest nie tak z $(file) "
        end
    end

    # Legend(fig[1,2], ax,  position = :rb, labelsize = 30)
    # display(fig)
    filenames = ["energy_plot_model1", "energy_plot_model2", "energy_plot_model3"]
    filename = filenames[model]
    save("lab03/graphics/$filename.pdf", fig)
end
get_energy_plot(model = 1)

## energia w czasie
function get_energy_plot_omega05(;model::Int = 1)
    init_vals = [
                L"\text{model 1}",
                L"\text{model 2}",
                L"\text{model 3}",
                ]
    fig = Figure(size = (1000, 600))
    ax = Axis(fig[1, 1], xlabel=L"$t$ \text{[s]}", ylabel=L"$|E(t)-E(0)|$ [\text{J}]", xlabelsize = 30,
    ylabelsize = 30, titlesize = 30, xticklabelsize = 28, yticklabelsize = 28,
    title=init_vals[model])#, limits = ((-0.03, 10.), (11.98, 12.02)))


    data_dir = "/home/marta/Documents/studia/dynamical_modelling/lab03/build/data/"
    if model == 1
        begin_file = r"pend_model1"
    elseif model == 2
        begin_file = r"pend_model2"
    else
        begin_file = r"pend_model3"
    end
    files = filter(f -> occursin(begin_file, f) && endswith(f, "w0_0.500000.txt"), readdir(data_dir; join=true))

    println("Found $(length(files)) data files:")
    foreach(println, files)

    for (i, file) in enumerate(files)
        data = readdlm(file, '\t', skipstart=1)
        label = init_vals[1]
        if size(data, 2) >= 3
            t_vals = data[:, 1]
            e_vals = abs.(data[:, 4] .-data[1,4])# energy 
            lines!(ax, t_vals, e_vals, linewidth=4, label = L"$(\vartheta_0, \omega_0)=(0.0, 0.5)$", alpha = 0.7)
        else
            @warn "cos jest nie tak z $(file) "
        end
    end

    axislegend(ax, position = :rb, labelsize = 30)
    # display(fig)
    filenames = ["energy_diff", "energy_diff-model2", "energy_diff-model3"]
    filename = filenames[model]
    save("lab03/graphics/$filename.pdf", fig)
    # display(fig)
end
get_energy_plot_omega05(model = 1)