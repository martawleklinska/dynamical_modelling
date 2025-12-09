using CairoMakie
using Printf
using DelimitedFiles

struct MapParams 
    K::Float64
    b::Float64
end

function step(θ::Float64, L::Float64, p::MapParams)
    if p.b == 0.0
        Lnew = L + p.K * sin(θ)
        θ_new = θ + Lnew
    else 
        exponen = exp(-p.b)
        Lnew = (L + p.K * sin(θ)) * exponen
        θ_new = θ + (L + p.K * sin(θ)) * (1- exponen)/p.b
    end
    return (θ_new, Lnew)
end

function simulate_traj(θ0::Float64, L0::Float64, p::MapParams; N=1000)
    θ = θ0; L = L0
    θ_arr = zeros(Float64, N+1); L_arr = zeros(Float64, N+1)
    θ_arr[1] = θ; L_arr[1] = L
    for n in 1:N
        θ, L = step(θ, L, p)
        θ_arr[n+1] = θ
        L_arr[n+1] = L
    end
    return (θ_arr, L_arr)
end

function poincare_points(θ0::Float64, L0::Float64, p::MapParams; N=1000, discard=200)
    θ = θ0; L = L0
    pts = Float64[]
    ptsL = Float64[]
    for n in 1:N
        θ, L = step(θ, L, p)
        if n > discard
            push!(pts, mod(θ, 2π))
            push!(ptsL, L)
        end
    end
    return (pts, ptsL)
end


function run_ex1()
    Ns = 1000
    discard = 50
    Ks = [0.5, 0.971635, 5.0]
    bs = [0.0, 0.001, 0.01]
    
    initials = [
        (0.0, 0.1*2π),
        (0.2*2π, 0.0),
        (0.2*2π, 0.65*2π),
    ]
    a_grid = 0.0:0.02:0.98

    mkpath("lab08/data/") 

    for b in bs
        for K in Ks
            p = MapParams(K, b)
            prefix = @sprintf("lab08/data/K_%0.6f_b_%0.6f", K, b)

            fig = Figure(resolution=(800, 500))
            ax = Axis(fig[1,1], title="L_n - K=$(K), b=$(b)", xlabel = L"$n$", ylabel = L"L_n", ylabelsize = 20, xlabelsize = 20, titlesize = 20, xticklabelsize = 15, yticklabelsize = 15)
            ax2 = Axis(fig[1,2], title="θ_n - K=$(K), b=$(b)", xlabel = L"$n$", ylabel = L"$\vartheta$", ylabelsize = 20, xlabelsize = 20, titlesize = 20, yticklabelsize = 15)
            for (i, (θ0, L0)) in enumerate(initials)
                θarr, Larr = simulate_traj(θ0, L0, p, N=Ns)
                n = 0:Ns
                CairoMakie.scatter!(ax2, n, mod.(θarr, 2π), markersize = 3,  label=L"\theta_n=%$(round(θ0, digits=2))")
                CairoMakie.scatter!(ax, n, Larr, label="L_n=$(round(L0, digits=2))", markersize = 3)
            end
            axislegend(ax, labelsize = 20)
            axislegend(ax2, labelsize = 20)
            save(prefix*"_timeseries.png", fig)

            allθ = Float64[]
            allL = Float64[]
            for a in a_grid
                θ0 = 0.5*2π
                L0 = a*2π
                ptsθ, ptsL = poincare_points(θ0, L0, p, N=Ns, discard=discard)
                append!(allθ, ptsθ)
                append!(allL, ptsL)
            end
            fig = Figure(size = (800, 500))
            ax = Axis(fig[1,1],  xlabel=L"\vartheta \text{ mod} 2\pi", ylabel=L"L", title="Poincaré K=$(K), b=$(b)",
                    ylabelsize = 20, xlabelsize = 20, titlesize = 20, xticklabelsize = 15, yticklabelsize = 15)
            CairoMakie.scatter!(ax, allθ, allL, markersize=2)
            save(prefix*"_poincare.png", fig)
            
            writedlm(prefix*"_poincare.csv", hcat(allθ, allL), ',')
            println("Saved results for K=$(K) b=$(b)")
        end
    end
end


# run_ex1()

## ex4: do animation
function generate_K_gifs()
    Ks = 0.1:0.1:5.0
    bs = [0.0, 0.001, 0.01]
    Ns = 1000
    discard = 50
    
    a_grid = 0.0:0.02:0.98
    θ0 = 0.5*2π

    mkpath("lab08/gifs/")

    for b in bs
        filename = "lab08/gifs/poincare_b_$(round(b, digits=3)).gif"
        
        fig = Figure(size = (800, 500))
        title_obs = Observable("Poincaré Map - K=$(round(Ks[1], digits=1)), b=$(round(b, digits=3))")
        ax = Axis(fig[1,1], 
                 xlabel = L"\vartheta \text{ mod } 2\pi", 
                 ylabel = L"L", 
                 title = title_obs,
                 ylabelsize = 22, 
                 xlabelsize = 22, 
                 titlesize = 22, 
                 xticklabelsize = 18, 
                 yticklabelsize = 18,)
        
        allθ_obs = Observable(Float64[])
        allL_obs = Observable(Float64[])
        
        scatter_plot = CairoMakie.scatter!(ax, allθ_obs, allL_obs, markersize=2, color=:blue)
        
        xlims!(ax, 0, 2π)
        ylims!(ax, -2π+2.5, 2π+2.5)
        
        println("Generating for b = $b")
        
        record(fig, filename, 1:length(Ks), framerate = 9) do i 
            K = Ks[i]
            p = MapParams(K, b)
            
            empty!(allθ_obs.val)
            empty!(allL_obs.val)
            
            for a in a_grid
                L0 = a*2π
                ptsθ, ptsL = poincare_points(θ0, L0, p, N=Ns, discard=discard)
                append!(allθ_obs.val, ptsθ)
                append!(allL_obs.val, ptsL)
            end
            
            title_obs[] = "Poincaré Map - K=$(round(K, digits=1)), b=$(round(b, digits=3))"
            
            notify(allθ_obs)
            notify(allL_obs)
        end
        
        println("Saved GIF: $filename")
    end
    
    println("done")
end

generate_K_gifs()
