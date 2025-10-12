#include<iostream>
#include<fstream>
#include<cmath>
#include<vector>

class LogisticMap{
    public:
        double func(double x, double r); 
        std::vector<std::vector<double> > evolution;
        // std::vector<std::vector<double> > evo_r;
        std::vector<double> init_values = {0.2, 0.4, 0.5, 0.9, 0.999};
        // std::vector<double> r_vals;
        LogisticMap() : evolution(50, std::vector<double>(5)) {};  
        void do_evolution();
        void save_log_map_ex1();
        void save_log_map_ex2();
};

double LogisticMap::func(double x, double r){
    return r * x * (1-x);
}
void LogisticMap::do_evolution(){
    for (int j = 0; j < init_values.size(); j++){
        evolution[0][j] = init_values[j];
    }
    
    double r = 2.;  
    for (int i = 1; i < evolution.size(); i++) {        
        for (int j = 0; j < init_values.size(); j++){    
            evolution[i][j] = func(evolution[i-1][j], r);
        }
    }
}

void LogisticMap::save_log_map_ex1(){
    std::ofstream output("log_map.txt");
    do_evolution();  

    for (int i = 0; i < evolution.size(); i++){
        for (int j = 0; j < init_values.size(); j++){
            output << evolution[i][j] << "\t";
        }
    output << "\n";
    }
    output.close();
}

void LogisticMap::save_log_map_ex2() {
    std::vector<double> r_vals(300); // 1 to 4 step 0.01 → (4-1)/0.01 = 300
    for (int i = 0; i < 300; i++) {
        r_vals[i] = 1.0 + i * 0.01;
    }

    const int steps = 10000;
    std::vector<std::vector<double>> evolution(r_vals.size(), std::vector<double>(steps, 0.0));

    double x0 = 0.5;
    for (size_t i = 0; i < r_vals.size(); ++i)
        evolution[i][0] = x0;

    for (int n = 1; n < steps; ++n) {
        for (size_t i = 0; i < r_vals.size(); ++i) {
            double r = r_vals[i];
            evolution[i][n] = func(evolution[i][n - 1], r);
        }
    }

    std::ofstream output("log_map_r.txt");
    for (size_t i = 0; i < r_vals.size(); ++i) {
        for (int n = steps - 1000; n < steps; ++n) {
            output << r_vals[i] << "\t" << evolution[i][n] << "\n";
        }
    }
    output.close();
}


int main(){
    LogisticMap lm;
    // lm.save_log_map_ex1();

    lm.save_log_map_ex2();
    return 0;
}
