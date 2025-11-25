using CairoMakie

f(x, a) = a * x - x^3
fprime(x, a) = a - 3x^2

## ex2: x_n(x)
x0s = [-0.8, -0.3, 0.2, 0.7]
function iterate_map(a, x0, N)
    xs = Vector{Float64}(undef, N+1)
    xs[1] = x0
    for n in 1:N
        xs[n+1] = f(xs[n], a)
    end
    return xs
end

function mapn(a::Float64)
    N = 300
    n = 0:1:N
    fig = Figure(size = (800, 600))
    ax = Axis(fig[1,1], xlabel = L"$n$", ylabel = L"$x_n$", xlabelsize = 25, ylabelsize = 25, title = "a=$a", titlesize = 25)
    for x0 in x0s
        xs = iterate_map(a, x0, N)
        CairoMakie.scatter!(ax, n, xs, label="x0=$(x0)")
    end
    Legend(fig[1,2], ax)
    save("lab07/graphics/ex2$a.pdf", fig)
    return fig
end

as = [0.9, 1.1, 2.1, 2.46, 2.62, 2.83]
for a in as
    mapn(a)
end

## ex4: bifurcation
function bifurcation_diagram(a_min, a_max; na=500, x0s=collect(range(-1.5,1.5,length=40)),
                             transient=500, sample=200)
    fig = Figure(size = (1000, 600))
    ax = Axis(fig[1,1], xlabel = L"$a$", ylabel = L"$x$", xlabelsize = 25, ylabelsize = 25, 
              title = "Diagram bifurkacji", titlesize = 25)
    as = range(a_min, a_max, length=na)
    pts_a = Float64[]
    pts_x = Float64[]

    for (i, a) in enumerate(as)
        for x0 in x0s
            x = x0
            for i in 1:transient
                x = f(x,a)
            end
            for i in 1:sample
                x = f(x,a)
                push!(pts_a, a)
                push!(pts_x, x)
            end
        end
    end
    CairoMakie.scatter!(ax, pts_a, pts_x, markersize=1, alpha=0.7)
    
    save("lab07/graphics/bifurcation_diagram.png", fig)
    display(fig)
    return fig
end
bifurcation_diagram(0, 3)

## ex5: wykładnik lapunowa
function lapunow(a, x0; N = 300)
    x = f(x0, a)  
    sumln = 0.0
    for i in 1:N
        d = abs(fprime(x, a))
        d = max(d, 1e-12)  
        sumln += log(d)
        x = f(x, a)
    end
    
    return sumln / N
end

function get_lambdaa(a_min::Float64, a_max::Float64; x0=0.2 ,if_x0::Bool=true)
    x0s = [0.1, 0.5, 1.0, 1.5]
    Ns = [10, 50, 100, 300, 600]
    na = 400
    as = range(a_min, a_max, length = na)
    λ = zeros(na)  
    λ1 = zeros(na)  
    λ2 = zeros(na)  
    λ3 = zeros(na)  
    fig = Figure(size = (800, 500))
    ax = Axis(fig[1,1], xlabel = L"$a$", ylabel = L"$\lambda$", 
    title = L"\text{Wykładnik Lapunowa przy różnych } N", xlabelsize = 25, ylabelsize = 25, 
    titlesize = 25, limits = ((1,3), (-5, 1)))
    
    for (i, a) in enumerate(as)
        if if_x0
            λ[i] = lapunow(a, x0s[1])
            λ1[i] = lapunow(a, x0s[2])
            λ2[i] = lapunow(a, x0s[3])
            λ3[i] = lapunow(a, x0s[4])
        else
            λ[i] = lapunow(a, x0; N =Ns[1])
            λ1[i] = lapunow(a, x0; N = Ns[2])
            λ2[i] = lapunow(a, x0; N =Ns[3])
            λ3[i] = lapunow(a, x0; N = Ns[4])
        end
    end
    if if_x0
        x1=x0s[1]
        x2=x0s[2]
        x3=x0s[3]
        x4=x0s[4]
        lines!(ax, as, λ, label = "x0=$x1")
        lines!(ax, as, λ1, label = "x0=$x2")
        lines!(ax, as, λ2, label = "x0=$x3")
        lines!(ax, as, λ3, label = "x0=$x4")
    else
        N1 = Ns[1]
        N2 = Ns[2]
        N3 = Ns[3]
        N4 = Ns[4]
        N5 = Ns[5]
        lines!(ax, as, λ, label = "N=$N1")
        lines!(ax, as, λ1, label = "N=$N2")
        lines!(ax, as, λ2, label = "N=$N3")
        lines!(ax, as, λ3, label = "N=$N4")
    end
    hlines!(ax, [0], color=:red, linestyle=:dash, alpha=0.7) 
    Legend(fig[1,2], ax)
    save("lab07/graphics/ex7_lapunow$if_x0.png", fig)
    # display(fig)
    return fig
end
get_lambdaa(0., 3., if_x0=false)

## ex6
for a in as
    println("a=$a: lambda: ", lapunow(a, 1.5))
end

