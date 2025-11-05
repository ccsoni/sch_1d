#pragma once
#include <cstdio>
#include <cstdlib>
#include <cstdint>
#include <cassert>
#include <complex>
#include <iostream>
#include <iomanip>
#include <map>
#include <string>
#include <vector>

template <typename T>
inline T SQR(T x)
{
  return x*x;
}

template <typename T>
inline T CUBE(T x)
{
  return x*x*x;
}

template <typename T>
inline T QUAD(T x)
{
  return x*x*x*x;
}

using namespace std::complex_literals;

class run_param {
public:
  int32_t nmesh_x, nmesh_v;
  int32_t nstep;

  std::string model_name;

  double tnow, tend, dt_output;
  double dtime;

  double rho;  // dt/(dx)^2
  double hbar;

  // phase space resolution where sigma_x*sigma_v = hbar/2
  double sigma_x, sigma_v;

  double xmax, xmin;
  double delta_x;

  double vmax, vmin;
  double delta_v;

  double mass;
  double Kene, Wene;

  int output_indx;

  int noutput;
  std::vector<double> output_timing;

  double next_output_timing()
  {
    return output_timing[output_indx];
  }

  void print_help(void)
  {
    std::cout <<
      "Mandatory options:\n -T t_end\n -t dt_output\n -m model_name\n" <<
      std::endl;
    std::cout <<
      "Additional options:\n -d hbar\n -r dt/(dx^2)\n -N nmesh\n" <<
      std::endl;
  }

  void init_run(int argc, char **argv)
  {
    bool mandatory_options_missed = false;

    //parse option
    std::map<std::string, std::string> options;
    for(int i=1;i<argc;i++) {
      std::string key=argv[i];
      if(i+1<argc && argv[i+1][0]!='-') {
	options[key] = argv[i+1];
      }else{
	options[key] = "true";
      }
    }

    if(options.count("--help")) {
      print_help();
      std::exit(EXIT_SUCCESS);
    }

    if(options.count("-T")) {
      this->tend = std::stof(options["-T"]);
      printf("# tend = %14.6e\n", this->tend);
    }else{
      mandatory_options_missed = true;
    }

    if(options.count("-t")) {
      this->dt_output = std::stof(options["-t"]);
      printf("# dt_output = %14.6e\n", this->dt_output);
    }else{
      mandatory_options_missed = true;
    }

    if(options.count("-m")) {
      this->model_name = options["-m"];
      printf("# model_name : %s\n", this->model_name.c_str());
    }else{
      mandatory_options_missed = true;
    }

    if(options.count("-N")) {
      this->nmesh_x = std::stoi(options["-N"]);
    }else{
      this->nmesh_x = 128;
    }
    printf("# nmesh_x = %d\n", this->nmesh_x);

    if(options.count("-d")) {
      this->hbar = std::stof(options["-d"]);
    }else{
      this->hbar = 0.1;
    }
    printf("# hbar = %14.6e\n", this->hbar);

    if(options.count("-r")) {
      this->rho = std::stof(options["-r"]);
    }else{
      this->rho = 1.0;
    }
    printf("# rho = %14.6e\n", this->rho);

    if(mandatory_options_missed) {
      print_help();
      exit(EXIT_FAILURE);
    }

    // setting up the output timing
    double output_time = 0.0;
    do {
      output_time += this->dt_output;
      output_timing.push_back(output_time);
    }while(output_time + this->dt_output < this->tend);
    output_timing.push_back(this->tend);

    this->output_indx = 0;

    // list the output timing
    std::cout << "# output timing" << std::endl;
    noutput = 0;
    for(double time : output_timing) {
      std::cout << std::scientific << std::setprecision(6) << "# " << time << std::endl;
      noutput++;
    }
    printf("# noutput = %d\n", noutput);

  }
};

constexpr double ROOT_2PI=1.772453851;
constexpr double QUAD_ROOT_2PI=1.583233487;
constexpr double HALF_M_PI=(M_PI*0.5);

constexpr int32_t MODEL_NAME_LENGTH=128;

using complexf = std::complex<float>;
using complexd = std::complex<double>;
