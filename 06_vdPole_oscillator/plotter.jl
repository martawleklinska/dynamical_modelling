# using CairoMakie, LinearAlgebra

# function lambda_vs_m(ε; mmin=0.0, mmax=3.0, n=800)
#     ms = range(mmin, mmax, length=n)
#     re = zeros(n); im = zeros(n)
#     re2 = zeros(n); im2 = zeros(n)
#     for (i,m) in enumerate(ms)
#         # eigenvals from lambda^2 + ε(m^2-1)λ + 1 = 0
#         a = 1.0
#         b = ε*(m^2 - 1.0)
#         c = 1.0
#         disc = b^2 - 4*a*c
#         if disc >= 0
#             λ1 = (-b + sqrt(disc))/(2a)
#             λ2 = (-b - sqrt(disc))/(2a)
#             re[i] = real(λ1); im[i] = imag(λ1)
#             re2[i] = real(λ2); im2[i] = imag(λ2)
#         else
#             # complex conjugate pair
#             realpart = -b/(2a)
#             imagpart = sqrt(-disc)/(2a)
#             re[i] = realpart; im[i] = imagpart
#             re2[i] = realpart; im2[i] = -imagpart
#         end
#     end
#     return ms, re, im, re2, im2
# end

# for ε in (0.5)
#     ms, r1, im1, r2, im2 = lambda_vs_m(ε; mmin=0.0, mmax=3.0, n=800)
#     fig = Figure(resolution=(800,500))
#     ax = Axis(fig[1,1], xlabel=L"$m$", ylabel = L"\lambda", title="ε=$(ε)", xlabelsize = 25, ylabelsize = 25, titlesize = 25,
#             xticklabelsize = 20, yticklabelsize = 20)
#     lines!(ax, ms, r1, label=L"\mathfrak{Re} (\lambda_0)")
#     lines!(ax, ms, r2, linestyle=:dash, label=L"\mathfrak{Re} (\lambda_1)")
#     lines!(ax, ms, im1, color=:red, label=L"\mathfrak{Im} (\lambda_0)")
#     lines!(ax, ms, im2, color=:red, linestyle=:dash, label=L"\mathfrak{Im} (\lambda_1)")
#     # hlines!(ax, [0.0], color=:black, linestyle=:dot)
#     axislegend(ax, labelsize = 20)
#     # display(fig)
#     save("lab06/graphics/ex1_eps0.5_eigen.pdf", fig)
# end


# ## get_trajectories
# """
#         {0.1, 0.2},   // x0=0.1, y0=0.2
#         {-0.5, 1.0},  // x0=-0.5, y0=1.0
#         {2.0, -1.5},  // x0=2.0, y0=-1.5
#         {0.0, 0.0}    // x0=0.0, y0=0.0
# """
# function get_trajectories()
#     init_vals = [L"(x_0,\;v_0)=(0.1, \;0.2)", L"(x_0,\;v_0)=(-0.5, \;1.0)",
#      L"(x_0,\;v_0)=(2.0, \;-1.5)", L"(x_0,\;v_0)=(0.0, \;0.0)"]
#     xs = range(-1., 2.5, length = 20)
#     vs = range(-2.5, 2.5, length = 20)

#     X = [v for x in xs, v in vs]
#     V = [- 0.5 * (x^2-1) * v-x + 1.0 for x in xs, v in vs]
#     mag = sqrt.(X.^2 .+ V.^2) .* 10
#     X ./= mag
#     V ./= mag
    
#     # ============================================================
#     fig = Figure(size = (1000, 500))
#     ax = Axis(fig[1, 1], xlabel=L"x", ylabel=L"$v$", xlabelsize = 30, limits = ((-1.,2.5), (-1.5, 2.5)),
#     ylabelsize = 30, title=L"\epsilon=0.5, \; m = 1.0", titlesize = 30, xticklabelsize = 20, yticklabelsize = 20)
#     arrows!(ax, xs, vs, X, V, arrowsize=8, lengthscale = 1.5, linecolor=:gray, linewidth=8, alpha = 0.6)
#     # ============================================================
#     data_dir = "/home/marta/Documents/studia/dynamical_modelling/lab06/build/data/"
#     begin_file = "eps05m1_"

#     files = filter(f -> occursin(begin_file, f) && endswith(f, ".txt"), readdir(data_dir; join=true))

#     println("Found $(length(files)) data files:")
#     foreach(println, files)

#     # cm = cgrad(:Blues, 6)
#     cm = cgrad(:RdBu_6, 6)

#     for (i, file) in enumerate(files)
#         data = readdlm(file, '\t', skipstart=1)
#         label = init_vals[i]
#         if size(data, 2) >= 3
#             x_vals = data[:, 2]  
#             y_vals = data[:, 3]  
#             lines!(ax, x_vals, y_vals, linewidth=4, label = label, alpha = 0.7, color = cm[i])
#         else
#             @warn "cos jest nie tak z $(file) "
#         end
#     end

#     scatter!(ax, 1., 0.0, color = :black, markersize = 20, marker = :star5)
#     # if p.alpha * p.beta > 0
#     #     scatter!(ax, sqrt(p.alpha/p.beta), 0.0, color = :blue, markersize = 20, marker = :star5)
#     #     scatter!(ax, -sqrt(p.alpha/p.beta), 0.0, color = :red, markersize = 20, marker = :star5)
#     # end
#     Legend(fig[2,1], ax, labelsize = 25, framevisible = false, orientation = :horizontal, nbanks = 2)
#     display(fig)
#     # save("lab06/graphics/eps05m1.pdf", fig)
#     return fig
# end
# get_trajectories()

# ## get time dependencies
# function get_time_dependencies(filename)
#     init_vals = [L"(x_0,\;v_0)=(0.1, \;0.2)", L"(x_0,\;v_0)=(-0.5, \;1.0)",
#      L"(x_0,\;v_0)=(2.0, \;-1.5)", L"(x_0,\;v_0)=(0.0, \;0.0)"]
    
#      fig = Figure(size = (1000, 500))
#     ax = Axis(fig[1, 1], xlabel=L"t", ylabel=L"$x$", xlabelsize = 30, limits = ((-0.1, 20.2), (-1.2, 3.2)),
#     ylabelsize = 30, title="eps=$filename", titlesize = 30, xticklabelsize = 20, yticklabelsize = 20)

#     # ============================================================
#     data_dir = "/home/marta/Documents/studia/dynamical_modelling/lab06/build/data/"
#     begin_file = filename

#     files = filter(f -> occursin(begin_file, f) && endswith(f, ".txt"), readdir(data_dir; join=true))

#     println("Found $(length(files)) data files:")
#     foreach(println, files)

#     # cm = cgrad(:RdYlBu_11, 6)
#     # cm = cgrad(:RdBu_6, 6)
#     cm = cgrad(:Dark2_6, length(files))
#     for (i, file) in enumerate(files)
#         data = readdlm(file, '\t', skipstart=1)
#         label = init_vals[i]
#         if size(data, 2) >= 3
#             x_vals = data[:, 2]  
#             v_vals = data[:, 3]  
#             t_vals = data[:, 1]
#             lines!(ax, t_vals, x_vals, linewidth=4, label = label, alpha = 0.7, color = cm[i])
#             # lines!(ax2, t_vals, v_vals, linewidth=4, label = label, alpha = 0.7, color = cm[i])
#         else
#             @warn "cos jest nie tak z $(file) "
#         end
#     end

#     Legend(fig[2,1], ax, labelsize = 25, orientation = :horizontal, nbanks = 2, framevisible = false)
#     display(fig)
#     # save("lab06/graphics/time$filename.pdf", fig)
#     return fig

# end
# get_time_dependencies("eps05m1_")


# Automated plotting utilities for eps and m
using CairoMakie, LinearAlgebra, DelimitedFiles

# ------------------------------------------------------------
# Helper: create directory string in correct format
# ------------------------------------------------------------
function eps_m_dir(ε, m)
    if ε == 0.0
        return "eps0"  # no m subfolders
    else
        # eps05..., eps3_...
        eps_str = ε == 0.5 ? "eps05_" : ε == 3.0 ? "eps3_" : error("unsupported ε")
        return string(eps_str, "m", m)
    end
end

# ------------------------------------------------------------
# Helper: TeX title formatter
# ------------------------------------------------------------
function tex_title(ε, m)
    if ε == 0.0
        return L"\varepsilon = 0.0"
    elseif ε == 3.0 && m == 0.5
        return L"\varepsilon = 3.0,\; m = 0.5"
    elseif ε == 3.0 && m == 0.8
        return L"\varepsilon = 3.0,\; m = 0.8"
    elseif ε == 3.0 && m == 1.0
        return L"\varepsilon = 3.0,\; m = 0.8"
    elseif ε == 3.0 && m == 1.1
        return L"\varepsilon = 3.0,\; m = 1.1"
    elseif ε == 3.0 && m == 2.0
        return L"\varepsilon = 3.0,\; m = 2.0"
    elseif ε == 3.0 && m == 2.5
        return L"\varepsilon = 3.0,\; m = 2.5"
        # ==== eps 0.5
    elseif ε == 0.5 && m == 0.5
        return L"\varepsilon = 0.5,\; m = 0.5"
    elseif ε == 0.5 && m == 1.0
        return L"\varepsilon = 0.5,\; m = 1.0"
    elseif ε == 0.5 && m == 1.5
        return L"\varepsilon = 0.5,\; m = 1.5"
    elseif ε == 0.5 && m == 2.0
        return L"\varepsilon = 0.5,\; m = 2.0"
    elseif ε == 0.5 && m == 2.5
        return L"\varepsilon = 0.5,\; m = 2.5"
    else
        return L"\varepsilon = ε,\; m = m"
    end
end

# ------------------------------------------------------------
# Function: compute eigenvalues as before
# ------------------------------------------------------------
function lambda_vs_m(ε; mmin=0.0, mmax=3.0, n=800)
    ms = range(mmin, mmax, length=n)
    re = zeros(n); im = zeros(n)
    re2 = zeros(n); im2 = zeros(n)
    for (i,m) in enumerate(ms)
        a = 1.0
        b = ε*(m^2 - 1.0)
        c = 1.0
        disc = b^2 - 4*a*c
        if disc >= 0
            λ1 = (-b + sqrt(disc))/(2a)
            λ2 = (-b - sqrt(disc))/(2a)
            re[i] = real(λ1); im[i] = imag(λ1)
            re2[i] = real(λ2); im2[i] = imag(λ2)
        else
            realpart = -b/(2a)
            imagpart = sqrt(-disc)/(2a)
            re[i] = realpart; im[i] = imagpart
            re2[i] = realpart; im2[i] = -imagpart
        end
    end
    return ms, re, im, re2, im2
end

# ------------------------------------------------------------
# Automated eigenvalue plot generator
# ------------------------------------------------------------
function plot_eigenvalues(ε, m)
    dir = eps_m_dir(ε, m)
    title_tex = tex_title(ε, m)

    ms, r1, im1, r2, im2 = lambda_vs_m(ε)

    fig = Figure(resolution=(800,500))
    ax = Axis(fig[1,1], xlabel=L"m", ylabel=L"\lambda", title=title_tex,
        xlabelsize=25, ylabelsize=25, titlesize=25,
        xticklabelsize = 20, yticklabelsize = 20, limits = ((0.0, 3.0), (-7, 3.0)))

    lines!(ax, ms, r1, label=L"\Re(\lambda_0)")
    lines!(ax, ms, r2, linestyle=:dash, label=L"\Re(\lambda_1)")
    lines!(ax, ms, im1, color=:red, label=L"\Im(\lambda_0)")
    lines!(ax, ms, im2, color=:red, linestyle=:dash, label=L"\Im(\lambda_1)")

    axislegend(ax, labelsize=20, position = :rc)
    save("lab06/graphics/eigen_ε$(ε)_m$(m).pdf", fig)
    # display(fig)
    return fig
end
# plot_eigenvalues(3, 0.5)

# ------------------------------------------------------------
# Trajectory plot automation
# ------------------------------------------------------------
function get_trajectories(ε, m)
    dir = eps_m_dir(ε, m)
    prefix = dir 
    title_tex = tex_title(ε, m)

    init_vals = [
                L"(x_0, v_0)=(-0.5,1.0)",
                L"(x_0, v_0)=(0.0,0.0)",
                L"(x_0, v_0)=(0.1,0.2)", 
                 L"(x_0, v_0)=(2.0,-1.5)", 
                 ]

    data_dir = "/home/marta/Documents/studia/dynamical_modelling/lab06/build/data/"
    files = filter(f -> occursin(prefix, f) && endswith(f, ".txt"), readdir(data_dir; join=true))
    
    xs = range(-2., 2., length=20)
    vs = range(-4., 5.5, length=20)

    X = [v for x in xs, v in vs]
    V = [-ε*(x^2-1)*v - x + m for x in xs, v in vs]
    mag = sqrt.(X.^2 .+ V.^2) .* 7
    X ./= mag; V ./= mag

    fig = Figure(size=(800,500))
    ax = Axis(fig[1,1], xlabel=L"x", ylabel=L"v", title=title_tex,
              xlabelsize=30, ylabelsize=30, titlesize=30, xticklabelsize=20, yticklabelsize=20)

    arrows!(ax, xs, vs, X, V, arrowsize=8, lengthscale=1.2,
            linecolor=:gray, linewidth=8, alpha=0.6)

    cm = cgrad(:rainbow, length(files))

    for (i, file) in enumerate(files)
        data = readdlm(file, '\t', skipstart=1)
        label = init_vals[i]
        if size(data,2) >= 3
            x_vals = data[:,2]; y_vals = data[:,3]
            lines!(ax, x_vals, y_vals, linewidth = 2, label = label)
        end
    end

    Legend(fig[2,1], ax, labelsize=25, framevisible=false, orientation=:horizontal, nbanks=2)
    save("lab06/graphics/trajectories_ε$(ε)_m$(m).png", fig)
    # display(fig)
    return fig
end
get_trajectories(3.0, 1.0)

# ------------------------------------------------------------
# Time dependencies plot
# ------------------------------------------------------------
function get_time_dependencies(ε, m)
    dir = eps_m_dir(ε, m)
    prefix = dir 
    title_tex = tex_title(ε, m)

    init_vals = [
                L"(x_0, v_0)=(-0.5,1.0)",
                L"(x_0, v_0)=(0.0,0.0)",
                L"(x_0, v_0)=(0.1,0.2)", 
                 L"(x_0, v_0)=(2.0,-1.5)", 
                 ]
    fig = Figure(size=(1000,500))
    ax = Axis(fig[1,1], xlabel=L"t", ylabel=L"x", title=title_tex,
              xlabelsize=30, ylabelsize=30, titlesize=30,
              limits=((-0.1,5.2),(-1.8,3.)), xticklabelsize=20, yticklabelsize=20)

    data_dir = "/home/marta/Documents/studia/dynamical_modelling/lab06/build/data/"
    files = filter(f -> occursin(prefix, f) && endswith(f, ".txt"), readdir(data_dir; join=true))

    cm = cgrad(:Dark2_6, length(files))

    for (i,file) in enumerate(files)
        data = readdlm(file, '\t', skipstart=1)
        if size(data,2) >= 3
            t_vals = data[:,1]; x_vals=data[:,2]
            lines!(ax, t_vals, x_vals, linewidth=4, label=init_vals[i], alpha=0.7, color=cm[i])
        end
    end

    Legend(fig[2,1], ax, labelsize=25, orientation=:horizontal, nbanks=2, framevisible=false)
    save("lab06/graphics/time_ε$(ε)_m$(m).pdf", fig)
    # display(fig)
    return fig
end

# get_time_dependencies(0.5, 0.5)