#include "crpropa/magneticField/timedependentField/MagneticFieldTimeDecay.h"

using namespace crpropa;

MagneticFieldTimeDecay::MagneticFieldTimeDecay(ref_ptr<MagneticField> field, double t0, double n) :
        magneticField(field) {
    setDecayTime(t0);
    setDecayExponent(n);
}

Vector3d MagneticFieldTimeDecay::getField(const Vector3d &pos, double z, double t) const {
    Vector3d B(0, 0, 0);
    if (magneticField.valid()) {
        B = magneticField->getField(pos, z, t);
        double decayFactor = pow((1 + t / t0), -n);
        B *= decayFactor; // Apply the decay factor to the magnetic field
    }

    return B;
}

void MagneticFieldTimeDecay::setField(ref_ptr<MagneticField> field) {
    magneticField = field;
}

ref_ptr<MagneticField> MagneticFieldTimeDecay::getField() const {
    return magneticField;
}

void MagneticFieldTimeDecay::setDecayTime(double t0) {
    if (t0 <= 0) {
        throw std::runtime_error("MagneticFieldTimeDecay: Decay time scale (t0) must be greater than 0");
    }
    this->t0 = t0;
}

double MagneticFieldTimeDecay::getDecayTime() const {
    return t0;
}

void MagneticFieldTimeDecay::setDecayExponent(double n) {
    if (n <= 0) {
        throw std::runtime_error("MagneticFieldTimeDecay: Decay exponent (n) must be greater than 0");
    }
    this->n = n;
}

double MagneticFieldTimeDecay::getDecayExponent() const {
    return n;
}

std::string MagneticFieldTimeDecay::getDescription() const {
    std::string s = "Magnetic field with time decay: B(t) = B0 * (1 + t / t0)^(-n)";
    s += "\n  Decay time scale (t0): " + std::to_string(t0);
    s += "\n  Decay exponent (n): " + std::to_string(n);
    if (magneticField.valid()) {
        s += "\n  Base magnetic field: " + magneticField->getDescription();
    } else {
        s += "\n  Base magnetic field: None";
    }
    
    return s;
}