include("plotter.jl")

function run_all_tasks()
    # run_ode_trajs()             # plot 9 considered cases of combinations of params with γ=0
    # run_time_ode()              # plot 9 considered cases od time dependencies with γ=0
    # run_example_force_trajs()   # plot one double well trajectories with γ≠0
    # run_example_time_force()    # plot one double well time-dependencies with γ≠0
    # plot_poincare("duffing")    # plot poincare map 
    # get_lapunow_exponen()       # plot Lapunov exponent with respect to γ
    # plot_bifur("duffing")       # plot bifurcation diagram with respect to γ
    run_energy_analysis()       # plot energy with respect to x and t
    # run_tdse()                  # SAVE gif of |ψ(x,t)|^2(x)
end

run_all_tasks()