#include <iostream>
#include <fstream>
#include <cmath>
#include <map>
#include <thread>
#include <mutex>
#include <limits>
#include <iomanip>

double function(double x, double y);
std::map<std::string, double> read_configuration(std::string path);
double volume(double x, double y, double dx, double dy);
double calculate_volume(double x_from, double x_to, double y_from, double y_to, double absolute_error, double relative_error);

const double a = 20;
const double b = 0.2;
const double c = 2*M_PI;

int main(int argc, char *argv[]) {

    if (argc < 2) {
        std::cerr << "No input file\nThe input should be like:\nprogram_name.exe config_file_name_and_path" << std::endl;
        return 1;
    }

    std::map<std::string, double> configuration;

    configuration = read_configuration(argv[1]);

    double res = 0;

    



//    ___________________________________

//    for (auto &x: configuration) {
//        std::cout << x.first << " " << x.second << std::endl;
//    }
//    std::cout << "Hello, World!" << std::endl;

    res = calculate_volume(
            configuration["x_from"],
            configuration["x_to"],
            configuration["y_from"],
            configuration["y_to"],
            configuration["absolute_error"],
            configuration["relative_error"]);

    std::cout << std::setprecision (std::numeric_limits<double>::digits10 + 1) << res << std::endl;

//    ___________________________________


    return 0;
}

double function(double x, double y) {
    return  -1 * a * exp(-1 * b * sqrt(0.5 * (x*x + y*y) ) ) -
            exp( 0.5 * ( cos(c*x) + cos(c*y) ) ) +
            a +
            exp(1);
}

std::map<std::string, double> read_configuration(std::string path) {
    std::ifstream file(path);

    std::map<std::string, double> configuration;

    if (!file) {
        std::cerr << "Error opening input file" << std::endl;
        return configuration;
    }

    std::string category_name;
    double category_value;
    while ( file >> category_name) {
        if (!(file >> category_value)) {
            std::cerr << "Error reading configuration file\nNo value for " << category_name << std::endl;
            return configuration;
        }
        configuration[category_name] = category_value;
    }
}

double volume(double x, double y, double dx, double dy) {
    return dx*dy*function(x, y);
}

double calculate_volume(double x_from, double x_to, double y_from, double y_to, double absolute_error, double relative_error) {

    double previous_res = std::numeric_limits<double>::infinity();
    double res = 0;
    int n = std::pow(2, 4);
    while (fabs(previous_res - res) > absolute_error or fabs(2 * (previous_res - res) / (previous_res + res)) > relative_error) {
        previous_res = res;
        res = 0;
        double dx = (x_to - x_from)/n;
        double dy = (y_to - y_from)/n;
        for (size_t i_x = 0; i_x < n; ++i_x) {
            for (size_t i_y = 0; i_y < n; ++i_y)  {
                res += volume(i_x*dx + x_from, i_y*dy + y_from, dx, dy);
            }
        }
        n *= 2;
    }
    std::cout << n << std::endl;

    return res;

}