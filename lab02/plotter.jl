using DelimitedFiles
using CairoMakie

data_y05 = readdlm("/home/marta/Documents/studia/dynamical_modelling/lab02/build/data/lv_rk4_y0_0.500000.txt", '\t', skipstart = 1)

function plot_lotka_volterra(data)
    t = data[:, 1]
    x = data[:, 2]
    y = data[:, 3]
    
    fig = Figure(resolution = (800, 600))
    
    ax1 = Axis(fig[1, 1], xlabel = "x (prey)", ylabel = "y (predator)", title = "Phase Portrait")
    lines!(ax1, x, y, linewidth = 2, color = :blue)
    scatter!(ax1, [x[1]], [y[1]], color = :green, markersize = 10, label = "Start")
    scatter!(ax1, [x[end]], [y[end]], color = :red, markersize = 10, label = "End")
    axislegend(ax1)
    
    ax2 = Axis(fig[1, 2], xlabel = "Time", ylabel = "Population", title = "Time Evolution")
    lines!(ax2, t, x, label = "Prey (x)", linewidth = 2, color = :blue)
    lines!(ax2, t, y, label = "Predator (y)", linewidth = 2, color = :red)
    axislegend(ax2)
    
    display(fig)
    return fig
end

function plot_vector_field(data)
    t = data[:, 1]
    x = data[:, 2]
    y = data[:, 3]
    
    step = max(1, length(x) ÷ 20)  
    indices = 1:step:length(x)
    
    x_sub = x[indices]
    y_sub = y[indices]
    
    us = [x[i+1]-x[i] for i in indices if i < length(x)]
    vs = [y[i+1]-y[i] for i in indices if i < length(y)]
    
    n = min(length(x_sub), length(us))
    x_sub = x_sub[1:n]
    y_sub = y_sub[1:n]
    us = us[1:n]
    vs = vs[1:n]
    
    fig = Figure(resolution = (600, 600))
    ax = Axis(fig[1, 1], xlabel = "x ", ylabel = "y", title = "Vector Field")
    
    # lines!(ax, x, y, linewidth = 1, color = :lightblue, alpha = 0.7)
    
    arrows2d!(ax, x_sub, y_sub, us, vs)
    
    display(fig)
    return fig
end

# plot_lotka_volterra(data_y05)
plot_vector_field(data_y05)