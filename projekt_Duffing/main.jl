include("plotter.jl")

function run_all_tasks()
    run_ode_trajs()
    run_time_ode()
    run_example_force_trajs()
    run_example_time_force()
    plot_poincare("duffing")
    get_lapunow_exponen()
    plot_bifur("duffing")
    run_energy_analysis()
end

run_all_tasks()