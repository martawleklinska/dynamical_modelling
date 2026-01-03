using Plots, DelimitedFiles, Printf
using StatsBase
gr()  

# duffing params
alpha = -1e-2
beta  =  1e-4
gamma =  0.5
omega =  1e-3
m     =  1.0
dt = 0.001


function create_wigner_animation()
    output_paths = [
        "build/output/",  
        "dynamical_modelling/projekt_Duffing/moyal_solver/build/output/",               
        "projekt_Duffing/moyal_solver/build/output/"             
    ]
    
    output_dir = nothing
    for path in output_paths
        if isdir(path)
            output_dir = path
            println("Found output directory at: $path")
            break
        end
    end
    
    if output_dir === nothing
        error("Output directory not found! Make sure you're in the moyal_solver/ directory and have run the C++ program first.")
    end
    
    wigner_files = filter(f -> startswith(f, "wigner_") && endswith(f, ".dat"), 
                         readdir(output_dir))
    sort!(wigner_files)
    
    if isempty(wigner_files)
        error("No Wigner function files found in $output_dir")
    end
    
    println("Found $(length(wigner_files)) Wigner function files")
    
    first_file = joinpath(output_dir, wigner_files[1])
    data = readdlm(first_file)
    x_coords = data[:, 1]
    p_coords = data[:, 2]
    
    x_unique = sort(unique(x_coords))
    p_unique = sort(unique(p_coords))
    nx = length(x_unique)
    np = length(p_unique)
    
    println("Grid size: $nx × $np")
    println("x range: $(minimum(x_unique)) to $(maximum(x_unique))")
    println("p range: $(minimum(p_unique)) to $(maximum(p_unique))")

    w_min, w_max = 0.3, -0.1
    for (i, filename) in enumerate(wigner_files)
        if i % 10 == 1  # Check every 10th file 
            data = readdlm(joinpath(output_dir, filename))
            w_vals = data[:, 3]
            w_min = min(w_min, minimum(w_vals))
            w_max = max(w_max, maximum(w_vals))
        end
    end
    println("Wigner function range: $w_min to $w_max")
    
    println("Creating animation...")
    anim = Animation()
    
    animation_files = wigner_files[1:4:end]
    
    for (i, filename) in enumerate(animation_files)
        
        step_str = match(r"wigner_(\d+)\.dat", filename).captures[1]
        step = parse(Int, step_str)
        time_val = step * 1.0  
        time = step * dt
        
        data = readdlm(joinpath(output_dir, filename))
        wigner_real = data[:, 3]
        
        W = reshape(wigner_real, np, nx)'
        W = sign.(W) .* abs.(W).^(1/4)
        
        p = Plots.heatmap(
            x_unique, p_unique, W',
            xlabel = "Położenie x",
            ylabel = "Pęd p",
            title  = @sprintf("Funkcja Wignera (t = %.1f)", time_val),
            c = :RdBu,
            clims  = (w_min, w_max),
                # xlims = (-4.5, 3),
                # ylims = (-6, 6),
            size   = (800, 600),
            dpi    = 100,
            aspect_ratio = :auto,
            )
            
        frame(anim, p)
        
        if i % 5 == 0
            println("frame $i/$(length(animation_files))")
        end
    end
    
    gif_filename = "projekt_Duffing/moyal_solver/graphics/wigner_evolution.gif"
    println("Saving animation to $gif_filename...")
    gif(anim, gif_filename, fps=5)
    
    println("Animation saved as $gif_filename")
    
    println("Creating key snapshots...")
    animation_files = wigner_files[1:20:end]

    for (i, filename) in enumerate(animation_files)
        
        step_str = match(r"wigner_(\d+)\.dat", filename).captures[1]
        step = parse(Int, step_str)
        time_val = step * 1.0  
        time = step * dt
        
        data = readdlm(joinpath(output_dir, filename))
        
        wigner_real = data[:, 3]
        W = reshape(wigner_real, np, nx)'
        W = sign.(W) .* abs.(W).^(1/4)
        time_val = step * 10.0
        
        p = Plots.heatmap(x_unique, p_unique, W',
                   xlabel="Położenie x", 
                   ylabel="Pęd p",
                   title=@sprintf("Funkcja Wignera t = %.1f", time_val),
                   color=:RdBu,
                   clims=(w_min, w_max),
                #    xlims = (-20, 40),
                #    ylims = (-0.7, 0.8),
                   size=(800, 600),
                   dpi=150)
        png_filename = @sprintf("projekt_Duffing/moyal_solver/graphics/wigner_snapshot_t%.0f.png", time_val)
        savefig(p, png_filename)
        println("✓ Saved $png_filename")
    end
    
    println("\n Visualization complete!")
    println("Files created:")
    println("  - wigner_snapshot_t*.png (key snapshots)")
end

create_wigner_animation()
##
function hamiltonian() 
    x_unique = range(-0, 30, 100)
    p_unique = range(-0.5, 0.5, 100)
    Vx = 0.5*alpha .* x_unique.^2 .+
        0.25*beta  .* x_unique.^4 .-
        gamma .* x_unique .* cos(omega * 0.0001)
    H = @. (p_unique'^2)/(2m) + Vx
    Emin = minimum(H)
    Emax = 0.0001 
    levels = range(Emin, Emax, length=50)  
    p = contour(
    x_unique, p_unique, H',
    levels = levels,
    linewidth = 1.0,
    linecolor = :black,
    alpha = 0.4
    )
return p
end
hamiltonian()