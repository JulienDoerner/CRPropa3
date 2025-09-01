#ifndef CRPROPA_MAGNETICFIELDTIMEDECAY_h
#define CRPROPA_MAGNETICFIELDTIMEDECAY_h

#include "crpropa/magneticField/MagneticField.h"


using namespace crpropa;

/**
 * \addtogroup MagneticField
 * @{
 */


/**
    @class MagneticFieldTimeDecay 
    @brief Magnetic field with a decay in time 

    The magnetic field will be scaled down depending on the time. 
    The magnetic field is given by the formula:
    B(t) = B0 * (1 + t / t0)^(-n)
 */

class MagneticFieldTimeDecay : public MagneticField {
private: 
    ref_ptr<MagneticField> magneticField; 
    double t0; //< decay time scale 
    double n;  //< decay exponent

public: 
    MagneticFieldTimeDecay(ref_ptr<MagneticField> field, double t0, double n = 6./5);

    Vector3d getField(const Vector3d &pos, double z, double t) const override;

    void setField(ref_ptr<MagneticField> field);
    ref_ptr<MagneticField> getField() const;

    void setDecayTime(double t0);
    double getDecayTime() const;

    void setDecayExponent(double n);
    double getDecayExponent() const;

    std::string getDescription() const;
};

/** @}*/

#endif // CRPROPA_MAGNETICFIELDTIMEDECAY_h