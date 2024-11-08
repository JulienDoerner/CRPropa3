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
    DiffusionPowerLaw(double alpha = 1./3, double rig_norm = 4 * GeV, double norm = 6.1e24, double epsilon = 0.1);
    Vector3d getTensorDiagonal(double E, int id = 0, double magneticField = 0.) const;

    void setAlpha(double alpha); 
    void setEpsilon(double eps);

    double getAlpha() const;
    double getEpsilon() const;
};

class DiffusionBrokenPowerlaw: public DiffusionTensor {
private: 
    double alpha1, alpha2; // indices before/after break 
    double E_break; // Energy of the break
    double norm; // norm at the break point 
    double epsilon; // anisotropy of the tensor

public: 
    DiffusionBrokenPowerlaw(double alpha1, double alpha2, double Eb, double norm, double eps= 0.1);

    void setAlpha1(double a1);
    void setAlpha2(double a2); 
    void setEnergyBreak(double Eb);
    void setNorm(double norm); 
    void setEpsilon(double eps); 

    Vector3d getTensorDiagonal(double E, int id = 0, double B = 0) const; 

    double getAlpha1() const;
    double getAlpha2() const;
    double getEnergyBreak() const;
    double getNorm() const;
    double getEpsilon() const;
};


class DiffusionMultipleBreaks: public DiffusionTensor {
private: 
    double initialIndex;            //< power law index before the first break
    double norm;                    //< normalisation at a given energy
    double Enorm;                   //< normalisation energy
    double eps;                     //< anisotropy of the diffusion tensor

    std::vector<double> index_s;    //< power law slope after break s
    std::vector<double> break_s;    //< energy of the break
    std::vector<double> norm_s;  //< normalisation at a given break 

public:
    DiffusionMultipleBreaks(double initialIndex, double eps, double norm = 6.1e24 , double Enorm = 4 * GeV);

    double getDiffusionCoefficient(double E) const;
    Vector3d getTensorDiagonal(double E, int id = 0, double B = 0) const; 

    /**
        add a break to the diffusion tensor
        @param index    power law slope after the break
        @param E        Energy of the break
    */
    void addBreak(double index, double E);


};

class DiffusionCRINGE : public DiffusionTensor {
private: 
    const double norm = 5.18e24; // norm of the diffusion tensor in m2/s
    const double slopes[5] = { 0.0116, 0.566, 0.159, 0.453, 1.050 }; // spectral indicies of the energy dependence
    const double breaks[4] = { pow(10,0.711), pow(10, 2.571), pow(10, 4.23), pow(10, 5.89)}; // rigidity [GV], position of the breaks
    const double softness[4] = { 0.0630, 0.75, 0.167, 0.022 }; // softness of the spectral breaks

public: 
    Vector3d getTensorDiagonal(double Energy, int id = 0, double B = 0) const;

    std::string getDescription() {
        return "Diffusion Tensor with the best fit of the CRINGE model.\n";
    }

};

class DiffusionBPLsoft : public DiffusionTensor {
  private: 
    double norm = 7.03e25; // norm of the diffusion tensor in m2/s
    double breakEnergy = 63 * GeV; // energy of the break
    double index1 = 0.334; // spectral index below the break
    double index2 = -0.324; // spectral index above the break
    double epsilon = 0.1; // anisotropy of the tensor

  public:

    DiffusionBPLsoft(double norm = 7.03e25, double breakEnergy = 63 * GeV, double index1 = 0.334, double index2 = -0.324, double epsilon = 0.1);

    Vector3d getTensorDiagonal(double Energy, int id = 0, double B = 0) const;

    std::string getDescription() {
        return "Diffusion Tensor with a broken power law with a soft break.\n";
    }

    void setEpsilon(double eps);
    double getEpsilon() const;

    void setIndex1(double a);
    double getIndex1() const;

    void setIndex2(double a);
    double getIndex2() const;

    void setBreakEnergy(double E);
    double getBreakEnergy() const;

    void setNorm(double n);
    double getNorm() const;
};

} // namespace crpropa

#endif // CRPROPA_DIFFUSIONTENSOR_H