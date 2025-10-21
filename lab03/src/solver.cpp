#include "solver.hpp"
#include <iostream>
#include <fstream>
#include <filesystem>
#include <cmath>


Solver_pendulum::Solver_pendulum(int model_type) : model(model_type) {}

int Solver_pendulum::func(double t, const double y[], double f[], void *params) {
    (void)t;
    Solver_pendulum *self = static_cast<Solver_pendulum*>(params);
    double theta = y[0];
    double omega = y[1];

    f[0] = omega; // dθ/dt = ω

    switch (self->model) {
        case 1: f[1] = -theta; break;                         // pierwsze przyblizenie
        case 2: f[1] = -(theta - pow(theta,3)/6.0); break;    // drugie przyblizenie
        case 3: f[1] = -sin(theta); break;                    // pełne
    }

    return GSL_SUCCESS;
}

void Solver_pendulum::jacobian(double theta, double omega, gsl_matrix *J, int model_type) {
    double dfdtheta = 0.0;
    double dfdomega = 1.0;
    double dgdtheta, dgdomega = 0.0;

    if (model_type == 1)
        dgdtheta = -1.0;
    else if (model_type == 2)
        dgdtheta = -(1 - theta*theta/2.0);
    else
        dgdtheta = -cos(theta);

    gsl_matrix_set(J, 0, 0, dfdtheta);
    gsl_matrix_set(J, 0, 1, dfdomega);
    gsl_matrix_set(J, 1, 0, dgdtheta);
    gsl_matrix_set(J, 1, 1, dgdomega);
}

void Solver_pendulum::solve(double t0, double t1, int n_steps) {
    std::filesystem::create_directories("data");

    gsl_odeiv2_system sys = {func, nullptr, 2, this};
    gsl_odeiv2_driver *driver = gsl_odeiv2_driver_alloc_y_new(
        &sys, gsl_odeiv2_step_rk8pd, (t1 - t0) / n_steps, abs_eps, rel_eps);

    std::vector<double> omega0_vals = {0.5, 1.0, 1.5, 2.0, 2.5, 3.0};
    double theta0 = 0.0;

    for (double omega0 : omega0_vals) {
        std::string fname = "data/pend_model" + std::to_string(model) + "_w0_" + std::to_string(omega0) + ".txt";
        std::ofstream ofs(fname);
        ofs << "t\ttheta\tomega\tE\n";

        double y[2] = {theta0, omega0};
        double t = t0;
        double m = 1.0, g = 9.8, l = 9.8;

        for (int i = 0; i <= n_steps; ++i) {
            double ti = t0 + i * (t1 - t0) / n_steps;
            int status = gsl_odeiv2_driver_apply(driver, &t, ti, y);
            if (status != GSL_SUCCESS) break;

            // energia
            double kinetic = 0.5 * m * l * l * y[1]*y[1];
            double potential = m * g * l * (1 - cos(y[0]));
            double E = kinetic + potential;

            ofs << ti << "\t" << y[0] << "\t" << y[1] << "\t" << E << "\n";
        }
        ofs.close();
    }

    gsl_odeiv2_driver_free(driver);
}

void Solver_pendulum::analyze_fixed_points() {
    std::vector<double> fixed_points = {0.0, M_PI}; // sinθ = 0
    gsl_matrix *J = gsl_matrix_alloc(2,2);

    for (double th : fixed_points) {
        jacobian(th, 0.0, J, model);

        double a = gsl_matrix_get(J,0,0);
        double b = gsl_matrix_get(J,0,1);
        double c = gsl_matrix_get(J,1,0);
        double d = gsl_matrix_get(J,1,1);

        double tr = a + d;
        double det = a*d - b*c;
        double disc = tr*tr - 4*det;

        std::cout << "Punkt staly: theta=" << th << " omega=0\n";
        std::cout << "  Tr=" << tr << "  Det=" << det << "  Δ=" << disc << "\n";

        if (det < 0) std::cout << "  -> Siodło\n";
        else if (tr == 0 && det > 0) std::cout << " -> Centrum\n";
        else if (tr < 0 && det > 0) std::cout << " -> Stabilny\n";
        else if (tr > 0 && det > 0) std::cout << " -> Niestabilny\n";
        std::cout << std::endl;
    }

    gsl_matrix_free(J);
}
