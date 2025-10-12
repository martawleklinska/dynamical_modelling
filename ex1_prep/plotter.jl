using DelimitedFiles
using CairoMakie

data = readdlm("ex1_prep/solution.csv", '\t', skipstart=1)

function plot_data(data)
    t = data[:, 1]          
    y_euler = data[:, 2]    
    y_rk4 = data[:, 3]      
    y_exact = data[:, 4]    
    
    fig = Figure(resolution = (800, 600))
    ax = Axis(fig[1,1], xlabel = "Time t", ylabel = "y(t)")
    
    lines!(ax, t, y_euler, label = "Euler", linewidth = 2)
    lines!(ax, t, y_rk4, label = "RK4", linewidth = 2)
    lines!(ax, t, y_exact, label = "Exact", linewidth = 2, linestyle = :dash)
    
    axislegend(ax, position = :rt)
    
    display(fig)
    
    return fig
end

plot_data(data)

## oscillaotr

data = readdlm("ex1_prep/vector_field.csv", ',', skipstart=1)
function get_arrows()
    x, y, dx, dy = data[1, :], data[2, :], data[3, :], data[4, :]
    fig = Figure(resolution = (800, 600))

    ax = Axis(fig[1,1], xlabel = "x", ylabel = "y")
    arrows2d!(ax, x, y, dx, dy)
    display(fig);
end
get_arrows()


## logistic_map

logistic_data = readdlm("ex1_prep/log_map.txt")
function get_log_map()
    fig = Figure(resolution = (800, 600))
    ax = Axis(fig[1,1], xlabel = "Time step", ylabel = "x value")
    
    n_timesteps, n_initial_values = size(logistic_data)
    time_steps = 0:(n_timesteps-1)
    
    for i in 1:n_initial_values
        initial_val = logistic_data[1, i]  
        evolution = logistic_data[:, i]   
        lines!(ax, time_steps, evolution, label = "x₀ = $(initial_val)", linewidth = 2)
    end
    
    axislegend(ax, position = :rt)
    display(fig)
    return fig
end
get_log_map()

## next log_map
logistic_data = readdlm("ex1_prep/log_map_r.txt")
function get_log_map_r()
    r = logistic_data[:, 1]
    x = logistic_data[:, 2]

    fig = Figure(resolution = (800, 600))
    ax = Axis(fig[1,1], xlabel = "r", ylabel = "last 1000 x")
    
    scatter!(r, x)
    
    # axislegend(ax, position = :rt)
    display(fig)
    return fig
end
get_log_map_r()