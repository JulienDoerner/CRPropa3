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

DiffusionMultipleBreaks::DiffusionMultipleBreaks(double index, double eps, double norm, double Enorm) : 
    initialIndex(index), eps(eps), norm(norm), Enorm(Enorm) {
        index_s.push_back(index);
        break_s.push_back(Enorm);
        norm_s.push_back(norm);
    }

double DiffusionMultipleBreaks::getDiffusionCoefficient(double E) const {
    // find break iB = 0 is the first reference
    int iB = 0;

    while((iB + 1 < index_s.size())) { // check if next break exists
        if (E < break_s[iB + 1]) 
            break;
        iB++;
    }

    double r = E / break_s[iB];
    double D = norm_s[iB]; 
    double i = index_s[iB];


    return D * pow(r, i);
}

void DiffusionMultipleBreaks::addBreak(double index, double E) {
    if (E < break_s.back())
        throw std::runtime_error("next break must have a larger energy than all breaks before");
    
    norm_s.push_back(getDiffusionCoefficient(E));
    index_s.push_back(index);
    break_s.push_back(E);
}

Vector3d DiffusionMultipleBreaks::getTensorDiagonal(double E, int id, double B) const {
    return Vector3d(1, eps, eps) * getDiffusionCoefficient(E);
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

// ------------------------------------------------------------------------

DiffusionBPLsoft::DiffusionBPLsoft(double norm, double breakEnergy, double index1, double index2, double epsilon) : 
    norm(norm), breakEnergy(breakEnergy), index1(index1), index2(index2), epsilon(epsilon) {}

Vector3d DiffusionBPLsoft::getTensorDiagonal(double E, int id, double B) const {
    double x = E / breakEnergy;
    return norm * pow(x, index1) / (1 + pow(x, index2 - index1)) * Vector3d(1, epsilon, epsilon);
}

void DiffusionBPLsoft::setEpsilon(double eps) {
    epsilon = eps;
}
double DiffusionBPLsoft::getEpsilon() const {
    return epsilon;
}

void DiffusionBPLsoft::setNorm(double n) {
    norm = n;
}
double DiffusionBPLsoft::getNorm() const {
    return norm;
}

void DiffusionBPLsoft::setBreakEnergy(double E) {
    breakEnergy = E;
}
double DiffusionBPLsoft::getBreakEnergy() const {
    return breakEnergy;
}

void DiffusionBPLsoft::setIndex1(double a) {
    index1 = a;
}
double DiffusionBPLsoft::getIndex1() const {
    return index1;
}

void DiffusionBPLsoft::setIndex2(double a) {
    index2 = a;
}
double DiffusionBPLsoft::getIndex2() const {
    return index2;
}

// ------------------------------------------------------------------------