using Plots, DelimitedFiles, Printf
using StatsBase
gr()  

# duffing params
alpha = -1e-2
beta  =  1e-4
gamma =  0.5
omega =  1e-3
m     =  1.0
dt = 0.0001


function create_wigner_animation()
    output_paths = [
        "structured_moyal_solver/build/output/",  
        "output/",               
        "../output/"             
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
        error("Output directory not found! Make sure you're in the structured_moyal_solver/ directory and have run the C++ program first.")
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
    
    println("Determining global color scale...")
    w_min, w_max = 0.3, -0.1
    for (i, filename) in enumerate(wigner_files)
        if i % 10 == 1  # Check every 10th file for efficiency
            data = readdlm(joinpath(output_dir, filename))
            w_vals = data[:, 3]
            w_min = min(w_min, minimum(w_vals))
            w_max = max(w_max, maximum(w_vals))
        end
    end
    println("Wigner function range: $w_min to $w_max")
    
    println("Creating animation...")
    anim = Animation()
    
    animation_files = wigner_files[1:1:50]
    
    for (i, filename) in enumerate(animation_files)
        
        step_str = match(r"wigner_(\d+)\.dat", filename).captures[1]
        step = parse(Int, step_str)
        time_val = step * 10.0  
        time = step * dt
        
        
        data = readdlm(joinpath(output_dir, filename))
        wigner_real = data[:, 3]
        
        W = reshape(wigner_real, np, nx)'
        
        p = heatmap(
            x_unique, p_unique, W',
            xlabel = "Position x",
            ylabel = "Momentum p",
            title  = @sprintf("Wigner Function (t = %.1f)", time_val),
            c = :RdBu,
            clims  = (w_min, w_max),
            xlims = (-20, 40),
            ylims = (-0.7, 0.8),
            size   = (800, 600),
            dpi    = 100,
            aspect_ratio = :auto
            )
            
            x_uni = range(-0, 30, 100)
            p_uni = range(-0.5, 0.5, 100)
            Vx = 0.5*alpha .* x_uni.^2 .+
                0.25*beta  .* x_uni.^4 .-
                gamma .* x_uni .* cos(omega * time)
            H = @. (p_uni'^2)/(2m) + Vx
            Emin = minimum(H)
            Emax = 0.0001 
            levels = range(Emin, Emax, length=50)  
            contour!(
            x_uni, p_uni, H',
            levels = levels,
            linewidth = 1.0,
            linecolor = :black,
            alpha = 0.4
            )
        
        frame(anim, p)
        
        if i % 5 == 0
            println("Processed frame $i/$(length(animation_files))")
        end
    end
    
    gif_filename = "wigner_evolution.gif"
    println("Saving animation to $gif_filename...")
    gif(anim, gif_filename, fps=5)
    
    println("✓ Animation saved as $gif_filename")
    
    println("Creating key snapshots...")
    key_steps = [0, 000001, 00002, 00003, 00004]
    for step in key_steps
        filename = @sprintf("wigner_%08d.dat", step)
        filepath = joinpath(output_dir, filename)
        
        if isfile(filepath)
            data = readdlm(filepath)
            wigner_real = data[:, 3]
            W = reshape(wigner_real, np, nx)'
            time_val = step * 10.0
            
            p = heatmap(x_unique, p_unique, W',
                       xlabel="position x", 
                       ylabel="momentum p",
                       title=@sprintf("Wigner Function at t = %.1f", time_val),
                       color=:RdBu,
                       clims=(w_min, w_max),
                    #    xlims = (-20, 40),
                    #    ylims = (-0.7, 0.8),
                       size=(800, 600),
                       dpi=150)
            x_unique = range(-0, 30, 100)
            p_unique = range(-0.5, 0.5, 100)
            Vx = 0.5*alpha .* x_unique.^2 .+
                0.25*beta  .* x_unique.^4 .-
                gamma .* x_unique .* cos(omega * 0.0001)
            H = @. (p_unique'^2)/(2m) + Vx
            Emin = minimum(H)
            Emax = 0.0001 
            levels = range(Emin, Emax, length=50)  
            contour!(
            x_unique, p_unique, H',
            levels = levels,
            linewidth = 1.0,
            linecolor = :black,
            alpha = 0.4
            )
            png_filename = @sprintf("wigner_snapshot_t%.0f.png", time_val)
            savefig(p, png_filename)
            println("✓ Saved $png_filename")
        end
    end
    
    println("\n Visualization complete!")
    println("Files created:")
    println("  - $gif_filename (animated)")
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