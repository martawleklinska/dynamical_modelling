using DelimitedFiles
using CairoMakie

data_0 = readdlm("lab01/build/data/init_cond_0.000000.txt", '\t', skipstart=1)
data_05 = readdlm("lab01/build/data/init_cond_0.500000.txt", '\t', skipstart=1)
data_10 = readdlm("lab01/build/data/init_cond_1.000000.txt", '\t', skipstart=1)
data_20 = readdlm("lab01/build/data/init_cond_2.000000.txt", '\t', skipstart=1)

data8_0 = readdlm("lab01/build/data/init_cond8_0.000000.txt", '\t', skipstart=1)
data8_05 = readdlm("lab01/build/data/init_cond8_0.500000.txt", '\t', skipstart=1)
data8_10 = readdlm("lab01/build/data/init_cond8_1.000000.txt", '\t', skipstart=1)
data8_20 = readdlm("lab01/build/data/init_cond8_2.000000.txt", '\t', skipstart=1)

err_0 = readdlm("lab01/build/data/error_0.000000.txt", '\t')
err_05 = readdlm("lab01/build/data/error_0.500000.txt", '\t')
err_10 = readdlm("lab01/build/data/error_1.000000.txt", '\t')
err_20 = readdlm("lab01/build/data/error_2.000000.txt", '\t')


errel_0 = readdlm("lab01/build/data/errorrel_0.000000.txt", '\t')
errel_05 = readdlm("lab01/build/data/errorrel_0.500000.txt", '\t')
errel_10 = readdlm("lab01/build/data/errorrel_1.000000.txt", '\t')
errel_20 = readdlm("lab01/build/data/errorrel_2.000000.txt", '\t')


err8_0 = readdlm("lab01/build/data/error8_0.000000.txt", '\t')
err8_05 = readdlm("lab01/build/data/error8_0.500000.txt", '\t')
err8_10 = readdlm("lab01/build/data/error8_1.000000.txt", '\t')
err8_20 = readdlm("lab01/build/data/error8_2.000000.txt", '\t')


errel8_0 = readdlm("lab01/build/data/errorrel8_0.000000.txt", '\t')
errel8_05 = readdlm("lab01/build/data/errorrel8_0.500000.txt", '\t')
errel8_10 = readdlm("lab01/build/data/errorrel8_1.000000.txt", '\t')
errel8_20 = readdlm("lab01/build/data/errorrel8_2.000000.txt", '\t')

function plot_data()
    x = data_0[:, 1]   
    x_20 = data_20[:, 1]       
    init0 = data_0[:, 2]    
    init05 = data_05[:, 2]      
    init10 = data_10[:, 2]    
    init20 = data_20[:, 2]

    x8 = data8_0[:, 1]   
    x_208 = data8_20[:, 1]       
    init08 = data8_0[:, 2]    
    init058 = data8_05[:, 2]      
    init108 = data8_10[:, 2]    
    init208 = data8_20[:, 2]

    fig = Figure(resolution = (800, 600))
    ax = Axis(fig[1,1], xlabel = L"$x$", ylabel = L"$y(t)$", limits = ((0.0, 1.0), (-2., 5.)),
                xlabelsize = 30, ylabelsize = 30, title = L"$\frac{dy}{dx}=y^2(x)$", titlesize = 30,
                xticklabelsize = 20, yticklabelsize = 20)
    
    lines!(ax, x, init0, label = "init cond=0.0", linewidth = 2)
    lines!(ax, x, init05, label = "init cond=0.5", linewidth = 2)
    lines!(ax, x, init10, label = "init cond = 1.0", linewidth = 2)
    lines!(ax, x_20, init20, label = "init cond = 2.0", linewidth = 2)

        
    lines!(ax, x8, init08, label = "init cond8=0.0", linewidth = 2)
    lines!(ax, x8, init058, label = "init cond8=0.5", linewidth = 2)
    lines!(ax, x8, init108, label = "init cond 8= 1.0", linewidth = 2)
    lines!(ax, x_208, init208, label = "init cond8 = 2.0", linewidth = 2)
    

    axislegend(ax, position = :lt, labelsize = 20)
    
    display(fig)
    
    return fig
end

function plot_error()
    x = err_05[:, 1]   
    x_20 = err_20[:, 1]    
    error0 = err_0[:, 2]    
    error05 = err_05[:, 2]      
    error10 = err_10[:, 2]    
    error20 = err_20[:, 2]

    x8 = err8_0[:, 1]   
    x_208 = err8_20[:, 1]    
    error08 = err8_0[:, 2]    
    error058 = err8_05[:, 2]      
    error108 = err8_10[:, 2]    
    error208 = err8_20[:, 2]

    fig = Figure(resolution = (800, 600))
    ax = Axis(fig[1,1], xlabel = L"$x$", ylabel = L"\text{error}",
                xlabelsize = 30, ylabelsize = 30, title = L"$\frac{dy}{dx}=y^2(x)$", titlesize = 30,
                yscale = log, xticklabelsize = 20, yticklabelsize = 20)
    
    lines!(ax, x, error0, label = "init cond=0.0", linewidth = 2)
    lines!(ax, x, error05, label = "init cond=0.5", linewidth = 2)
    lines!(ax, x, error10, label = "init cond = 1.0", linewidth = 2)
    lines!(ax, x_20, error20, label = "init cond = 2.0", linewidth = 2)

    
    lines!(ax, x8, error08, label = "init cond pd=0.0", linewidth = 2)
    lines!(ax, x8, error058, label = "init cond pd=0.5", linewidth = 2)
    lines!(ax, x8, error108, label = "init cond pd= 1.0", linewidth = 2)
    lines!(ax, x_208, error208, label = "init cond pd= 2.0", linewidth = 2)
    

    axislegend(ax, position = :rb, labelsize = 20)
    
    display(fig)

end

function plot_rel_error()
    x = errel_0[:, 1]   
    x_20 = errel_20[:, 1]    
    error0 = errel_0[:, 2]    
    error05 = errel_05[:, 2]      
    error10 = errel_10[:, 2]    
    error20 = errel_20[:, 2]

    x8 = errel8_0[:, 1]   
    x_208 = errel8_20[:, 1]    
    error08 = errel8_0[:, 2]    
    error058 = errel8_05[:, 2]      
    error108 = errel8_10[:, 2]    
    error208 = errel8_20[:, 2]


    fig = Figure(resolution = (800, 600))
    ax = Axis(fig[1,1], xlabel = L"$x$", ylabel = L"\text{error}",
                xlabelsize = 30, ylabelsize = 30, title = L"$\frac{dy}{dx}=y^2(x)$", titlesize = 30,
                yscale = log)
    
    lines!(ax, x, error0, label = "init cond=0.0", linewidth = 2)
    lines!(ax, x, error05, label = "init cond=0.5", linewidth = 2)
    lines!(ax, x, error10, label = "init cond = 1.0", linewidth = 2)
    lines!(ax, x_20, error20, label = "init cond = 2.0", linewidth = 2)
  
    lines!(ax, x8, error08, label = "init cond pd=0.0", linewidth = 2)
    lines!(ax, x8, error058, label = "init cond pd=0.5", linewidth = 2)
    lines!(ax, x8, error108, label = "init cond pd= 1.0", linewidth = 2)
    lines!(ax, x_208, error208, label = "init cond pd= 2.0", linewidth = 2)
    

    axislegend(ax, position = :rb, labelsize = 20)
    
    display(fig)
end

# plot_data()
plot_error()
# plot_rel_error()
