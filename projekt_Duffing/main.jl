include("plotter.jl")

function run_all_tasks(;display::Bool = true)
    # run_ode_trajs(display)             # plot 9 considered cases of combinations of params with γ=0
    # run_time_ode(display)              # plot 9 considered cases od time dependencies with γ=0
    # run_example_force_trajs(display)   # plot one double well trajectories with γ≠0
    # run_example_time_force(display)    # plot one double well time-dependencies with γ≠0
    # plot_poincare(display)             # plot poincare map 
    # poincare_gamma_gif()               # SAVE a Poincare map gif in γ=[0.25, 0.5]
    # get_lapunow_exponen(display)       # plot Lapunov exponent with respect to γ
    # plot_bifur(display)                # plot bifurcation diagram with respect to γ
    run_energy_analysis(display)       # plot energy with respect to x and t
    # run_tdse()                         # SAVE gif of |ψ(x,t)|^2(x)
    # get_expectation_values(display)    # get <x>(t), <p>(t)
end

run_all_tasks(display = false)