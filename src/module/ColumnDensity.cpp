#include "crpropa/module/ColumnDensity.h"

ColumnDensity::ColumnDensity(ref_ptr<Density> dens) : Module(), density(dens) {
    setDescription("ColumnDensity");
}

void ColumnDensity::process(Candidate* cand) const {
    // candidate information 
    double step = cand -> getCurrentStep();
    Vector3d pos = cand -> current.getPosition();

    double nGas = density -> getNucleonDensity(pos); // current density

    // calculate column density from the current step
    double CD = nGas * step; 

    // add to existing column density
    if (cand -> hasProperty(key)) {
        double CD_old = cand -> getProperty(key);
        CD += CD_old;
    }

    // store column density as a property of the candidate
    cand -> setProperty(key, CD);
}

void ColumnDensity::setKey(std::string k) {
    key = k;
}
std::string ColumnDensity::getKey() const {
    return key;
}