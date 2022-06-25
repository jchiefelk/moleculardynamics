#include "molecularDynamics.h"

int main() {
  MolecularDynamics md;
  md.initialize();
  md.equilibrate();
  md.measure();
	
  return 0;
};
