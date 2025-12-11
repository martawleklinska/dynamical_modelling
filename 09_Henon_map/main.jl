include("tasks.jl")

as = [0.50, 1.10, 1.25, 1.40]
b = 0.3
N = 2000
xs, ys = henon_iter(0.0, 0.0, 1.4, b, 50000)

function run_all_tasks()
    # run_task1()
    # run_task2()
    # run_task3()
    # run_task4()
    # run_task5()
    # run_task6()
    # run_ex7()
    run_ex8()
end

run_all_tasks()
