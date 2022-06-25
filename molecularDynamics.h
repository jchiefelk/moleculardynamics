#include <iostream>
#include <cctype>
#include <cstring>
#include <cmath>
using namespace std;
/***
 * 1D Molecular Dynamics Code
 *****/

const double TWOPI = 6.283185307179586;
const int N = 21;
const int NH = 301;

class MolecularDynamics {
  public:
    // Constructors
    MolecularDynamics();
    void position_verlet(int N, double h, double x[], double v[]);
    double force(double y);
    void initialize();
    void equilibrate();
    void measure();


    // Destructors
    ~MolecularDynamics();

  private:
    double y[N];
    double v[N];
    double h;
    double phi;
    double t_for_equil;
    double t_for_meas;
    double sumv;
    double binsize;
    double v2;
    double u;
    int i;
    int time;
    int n_for_equil;
    int n_of_meas;
    int n_bet_meas;
    int istep;
    int i_of_meas;
    int ihist;
    int iv;
    double hist[2*NH + 1];
};
