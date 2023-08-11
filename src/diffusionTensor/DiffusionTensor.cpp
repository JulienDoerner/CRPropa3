#include "crpropa/diffusionTensor/DiffusionTensor.h"
#include "crpropa/Vector3.h"

using namespace crpropa;

DiffusionPowerLaw::DiffusionPowerLaw(double alpha, double rig, double norm, double epsilon) : rig_norm(rig), norm(norm) { 
        setAlpha(alpha);
        setEpsilon(epsilon);
    }

Vector3d DiffusionPowerLaw::getTensorDiagonal(double E, int id, double B) const {
    Vector3d diff = Vector3d(1, epsilon, epsilon); 
    diff *= norm * pow(E / rig_norm, alpha); 
    return diff;
}


void DiffusionPowerLaw::setAlpha(double alpha) {
    alpha = alpha;
}
void DiffusionPowerLaw::setEpsilon(double eps) {
    epsilon = eps;
}

double DiffusionPowerLaw::getAlpha() const {
    return alpha;
}

double DiffusionPowerLaw::getEpsilon() const {
    return epsilon;
}

// ------------------------------------------------------------------------

DiffusionBrokenPowerlaw::DiffusionBrokenPowerlaw(double a1, double a2, double Eb, double norm, double eps) {
    setAlpha1(a1); 
    setAlpha2(a2);
    setEnergyBreak(Eb);
    setNorm(norm);
    setEpsilon(eps);
}

void DiffusionBrokenPowerlaw::setAlpha1(double a1) {
    alpha1 = a1;
}

void DiffusionBrokenPowerlaw::setAlpha2(double a2) {
    alpha2 = a2;
}

void DiffusionBrokenPowerlaw::setEnergyBreak(double Eb) {
    E_break = Eb;
}

void DiffusionBrokenPowerlaw::setNorm(double n) {
    norm = n;
}

void DiffusionBrokenPowerlaw::setEpsilon(double e) {
    epsilon = e;
}

Vector3d DiffusionBrokenPowerlaw::getTensorDiagonal(double E, int id, double B) const{
    double r = E / E_break; 
    Vector3d diag(1, epsilon, epsilon);
    if (r < 1) 
        diag *= norm * pow(r, alpha1);
    else
        diag *= norm * pow(r, alpha2);
    
    return diag;
}

double DiffusionBrokenPowerlaw::getAlpha1() const {
    return alpha1;
}

double DiffusionBrokenPowerlaw::getAlpha2() const {
    return alpha2;
}

double DiffusionBrokenPowerlaw::getEnergyBreak() const {
    return E_break;
}

double DiffusionBrokenPowerlaw::getNorm() const {
    return norm;
}

double DiffusionBrokenPowerlaw::getEpsilon() const {
    return epsilon;
}

// ------------------------------------------------------------------------

Vector3d DiffusionCRINGE::getTensorDiagonal(double E, int id, double B) const {
    double gamma = E / mass_proton / c_squared;
    double beta = sqrt(1 - 1 / gamma / gamma); 
    double rigidity = E / GeV; // rigidity in GV for a proton or electron

    double D = norm * beta * pow(rigidity / breaks[0], slopes[0]);
    for(int i = 0; i < 4; i++) {
        D *= pow(1 + pow(rigidity / breaks[i], 1 / softness[i]), softness[i] * (slopes[i+1] - slopes[i]));
    }

    return Vector3d(D);
}