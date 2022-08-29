#include "molecularDynamics.h"

// This file contains the implementation of the 
// MD Class member functions

// Molecular Dynamics Constructor
MolecularDynamics::MolecularDynamics(){
  binsize = 0.05;
  t_for_equil = 20000;
  t_for_meas = 1000000; // 1000000
  h = 0.01;
  phi = 0.6;
  n_for_equil = t_for_equil / h;
  n_of_meas = t_for_meas;
  n_bet_meas = 1 / h; // only measure v every unit of time
  v2 = 0;
};

// Molecular Dynamics Destructor
MolecularDynamics::~MolecularDynamics(){
  binsize = 0;
  t_for_equil = 0;
  t_for_meas = 0;
  h = 0;
  phi = 0;
  n_for_equil = 0;
  n_of_meas = 0;
  n_bet_meas = 0; // only measure v every unit of time
  v2 = 0;
};

// Molecular Dynamics - Write to File
void MolecularDynamics::write_to_file() {
  ofstream fileOut;
  fileOut.open("md_analysis.csv", ios::app);
  if(fileOut) {
    fileOut << "Temperature," << v2 << endl;
    fileOut << "binsize," << "histogram," << "expected gaussian" << endl;
    for(ihist = 0; ihist < 2*NH + 1; ++ihist) {
      u = pow((ihist - NH)*binsize, 2);
      // Only print out non-zero entries
      if(hist[ihist] != 0 || fabs(u) < 9 * v2 ) {
        hist[ihist] /= binsize*n_of_meas*N; // normalize the histogram 
        // Wirte to file histogram and expected Gaussian
        fileOut << binsize*(ihist-NH) << "," <<  hist[ihist] << "," << exp(-0.5*u/v2) / (sqrt(2*3.14159*v2)) << endl;
      }
    };
  }
};

// Molecular Dynamics - Bin Velocity
void MolecularDynamics::bin() {
  for(i = 0; i < N; ++i){
    v2 += v[i]*v[i]; // collect data for <v^2>
    iv = NH + round(v[i]/ binsize); // add to histogram
    hist[iv] +=1;
  };
};

// Molecular Dynamics - Simulation Step
void MolecularDynamics::simulation_step() {
  for(istep = 0; istep < n_bet_meas; ++istep) {
    position_verlet(N, h, y, v);
  };
};

// Molecular Dynamics - Measure
void MolecularDynamics::measure() {
  // Steps for measurment, loop over meas
  for(i_of_meas = 0; i_of_meas < n_of_meas; ++i_of_meas) {
    // steps between measurement
    simulation_step();
    bin();
  };

  v2/= n_of_meas * N; // average v^2 is the temperature
  cout << "v P(v) P(v)(Boltz)" << endl;
  write_to_file();
};

// Molecular Dynamics - Equilibrate
void MolecularDynamics::equilibrate() {
  for(istep = 0; istep < n_for_equil; ++istep) {
    position_verlet(N, h, y, v);
  };
};

// Molecular Dynamics - Initialize
void MolecularDynamics::initialize(){
  // Intialize Historgram
  for(int i = 0; i < 2*NH + 1; ++i) {
    hist[i] = 0; // Initialize histogram
  };
  // Intitialize positions and velocity
  for(int i = 0; i < N; ++i){
    v[i] =  sin(TWOPI*i/N + phi); // initialize velocity
    y[i] = 0; // Initialize y's (positions)
  };
};


// Position Verlet Algorithm - Integrations
void MolecularDynamics::position_verlet(int N, double h, double x[], double v[]) {
  // half-step in x
  for(int i = 0; i < N; ++i) {
    x[i] += 0.5*h*v[i];
  };
  // Full-step in vv; v[0] and V[N - 1] treated separately 
  // becuase of periodic bc's
  v[0] += h*(force(x[0] - x[1]) - force(x[N-1] - x[0]));
  for(int i = 1; i < N-1; ++i) {
    if(isnan(v[i])){
      return;
    }
    v[i] += h*(force(x[i] - x[i+1]) - force(x[i-1] - x[i]));
  };

  v[N-1] += h * (force(x[N-1] - x[0]) - force(x[N-2] - x[N-1]));
  // half-step in x
  for(int i = 0; i < N; ++i) {
    x[i] += 0.5*h*v[i];
  };

};

// Computes Forces
double MolecularDynamics::force(double y) {
  return -y - y*y - y*y*y;
};
