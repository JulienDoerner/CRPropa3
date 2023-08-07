#ifndef CRPROPA_DIFFUSIONTENSOR_H
#define CRPROPA_DIFFUSIONTENSOR_H

#include "crpropa/Units.h"
#include "crpropa/Referenced.h"
#include "crpropa/ParticleID.h"
#include "crpropa/Vector3.h"
#include "crpropa/Common.h"

#include <vector>
#include <cmath>

namespace crpropa {

/*
    @class DiffusionTensor
    @brief Abstract base class for different diffusion tensor

    The function getTensorDiagonal returns the values of the diagonal of the diffusion tensor. 
    The first entry is alligned with the magnetic field and the others are perpendicular.
*/
class DiffusionTensor : public Referenced {
public:
    /*
        Returns the diagonal entries of the diffusion tensor. 
        The first element is parallel to the magnetic field line. The second and third entry are the perpendicluar directions.

        @param Energy: Energy of the particle 
        @param id: Id of the particle (to have access to the charge e.g.)
        @param magneticField: magnetic field strength at the position
    */
    virtual ~DiffusionTensor() {};

    virtual Vector3d getTensorDiagonal(double Energy, int id = nucleusId(1,1), double magneticField = 1 * muG) const {
        return Vector3d(0.);
    }

    virtual std::string getDescription() {
		return "Abstract Diffusion Tensor Module\n";
	}
};


class DiffusionPowerLaw : public DiffusionTensor {
private: 
    double alpha = 1./3; // powerlaw index
    double rig_norm = 4 * GeV; // normalisation Energy
    double norm = 6.1e24; // diffusion coefficient normalisation 
    double epsilon = 0.1; // ratio between parallel and perpendicular diffusion

public:
    DiffusionPowerLaw(double alpha, double rig_norm, double norm, double epsilon);

    Vector3d getTensorDiagonal(double E, int id = 0, double magneticField = 0.) const;
};


class DiffusionCRINGE : public DiffusionTensor {
private: 
    const double norm = 5.18e24; // norm of the diffusion tensor in m2/s
    const double slopes[5] = { 0.0116, 0.566, 0.159, 0.453, 1.050 }; // spectral indicies of the energy dependence
    const double breaks[4] = { pow(10,0.711), pow(10, 2.571), pow(10, 4.23), pow(10, 5.89)}; // rigidity [GV], position of the breaks
    const double softness[4] = { 0.0630, 0.75, 0.167, 0.022 }; // softness of the spectral breaks

public: 
    Vector3d getTensorDiagonal(double Energy) const;

    std::string getDescription() {
        return "Diffusion Tensor with the best fit of the CRINGE model.\n";
    }

};

} // namespace crpropa

#endif // CRPROPA_DIFFUSIONTENSOR_H