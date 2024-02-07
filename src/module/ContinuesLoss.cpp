#include "crpropa/module/ContinuesLoss.h"

#include <string>

using namespace crpropa;

PiZeroLoss::PiZeroLoss(ref_ptr<Density> dens) : Module(), density(dens) {
	setDescription("PiZeroLoss");
}

void PiZeroLoss::process(Candidate* cand) const {
	// candidate information 
    double step = cand -> getNextStep();
	Vector3d pos = cand -> current.getPosition();
    double E = cand -> current.getEnergy() / GeV; // current energy in GeV

    // current density
	double nHI = density -> getHIDensity(pos);
	double nH2 = density -> getH2Density(pos);
    double nGas = (nHI + nH2) * ccm; // total gas density in cm^-3

    // calculate loss 
	double dEdT = 3.85e-16 * nGas * pow(E, 1.28) * pow(E + 200, -0.2); // energy loss per time (GeV / sec)
    double dE = dEdT * step / c_light; // only relativistic movement, in GeV

    // reduce energy 
    cand -> current.setEnergy((E - dE) * GeV); 
}