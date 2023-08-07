#include "crpropa/diffusionTensor/DiffusionTensor.h"
#include "crpropa/Vector3.h"

using namespace crpropa;

DiffusionPowerLaw::DiffusionPowerLaw(double alpha, double rig, double norm, double epsilon) : 
    alpha(alpha), rig_norm(rig), norm(norm), epsilon(epsilon) { }

Vector3d DiffusionPowerLaw::getTensorDiagonal(double E, int id, double B) const {
    Vector3d diff = Vector3d(1, epsilon, epsilon); 
    diff *= norm * pow(E / rig_norm, alpha); 
    return diff;
}

Vector3d DiffusionCRINGE::getTensorDiagonal(double E) const {
    double gamma = E / mass_proton / c_squared;
    double beta = sqrt(1 - 1 / gamma / gamma); 
    double rigidity = E / GeV; // rigidity in GV for a proton or electron

    double D = norm * beta * pow(rigidity / breaks[0], slopes[0]);
    for(int i = 0; i < 4; i++) {
        D *= pow(1 + pow(rigidity / breaks[i], 1 / softness[i]), softness[i] * (slopes[i+1] - slopes[i]));
    }

    return Vector3d(D);
}