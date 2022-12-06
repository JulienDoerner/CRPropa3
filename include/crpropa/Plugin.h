#ifndef CRPROPA_PLUGIN_H
#define CRPROPA_PLUGIN_H

#include <crpropa/Source.h>
#include <crpropa/Module.h>
#include <crpropa/module/Observer.h>
#include <crpropa/advectionField/AdvectionField.h>
#include <crpropa/magneticField/turbulentField/ModulatedTurbulentField.h>
#include <crpropa/Random.h>

using namespace crpropa;

/**
 * @class SourceRadiusFromSFR
 * @brief Sample the radial position according to the measured SFR. 
 * 
 * The angle is choosen randomly uniform within 0 and 2 pi. The normalisation of the SFR is arbitrary.
 */
class SourceRadiusFromSFR: public crpropa::SourceFeature{
private:
	std::vector<double> radius;
	std::vector<double> sfr;
	double sfrMax, rMax;  // maximal values for eficient sampling
public:
	SourceRadiusFromSFR(std::string filename);
	void prepareParticle(ParticleState& state) const;
	void loadData(std::string filename);
	double getSFRMax();
	double getRMax();
};
/**
 * @class SourceZpositionGauss
 * @brief sampling z-Position of the source from a gausian-distribution 
 * 
 * Input is a (radial constant) scaleheight and a maximal z-position. Recomanded to be 5*z.
 */
class SourceZpositionGauss: public crpropa::SourceFeature {
private:
	double zMax;	// maximal range to sample
	double h;		// scale-height
public:
	SourceZpositionGauss(double zMax, double h);
	void prepareParticle(ParticleState& state) const;

	double getZMax() const;
	double getH() const;
	void setZMax(double zMax);
	void setH(double h);
};

/**
 * @class SynchrotronSelfCompton
 * @brief Energy loss due to synchrotron radiation and IC scattering
 * 
 * This module calculates the energy loss of electrons due to synchrotron radiation and inverse compton scattering 
 * as a continues energy loss process. All secondaries are neglected. This approach can be found e.g. Mulcahy et al. DOI: 10.1051/0004-6361/201628446
 */
class SynchrotronSelfCompton: public crpropa::Module {
private:
	double uRad;
	ref_ptr<crpropa::ModulatedTurbulentField> field;
public:
	SynchrotronSelfCompton(ref_ptr<crpropa::ModulatedTurbulentField> field, double Urad);

	void process(Candidate *cand) const;

	void setURad(double uRad);
	void setMagneticField(ref_ptr<crpropa::ModulatedTurbulentField> field);

	double getURad() const;

	double energyLoss(Vector3d pos, double E) const;
	std::string getDescription() const;
};

/**
 * @class ObserverScaleHeight
 * @brief Variable scale height for M51
 * 
 * Data for the radial dependence of the scaleheight are read from a file and interpolated in between. 
 * 
 */
class ObserverScaleHeight: public crpropa::ObserverFeature {
private:
	std::vector<double> rBin, zBin;
	double Rmax;
public:
	ObserverScaleHeight(double Rmax, std::string path);
	void loadData(std::string path);
	DetectionState checkDetection(Candidate *candidate) const;
};

/**
 * @class AdvectionFieldFromList
 * @brief Advection field in z-direction (away from source) with radial dependence read from a table
 */

// class AdvectionFieldFromList: public AdvectionField {
// private:
//     std::vector<double> radius, velocity;
// public:
//     AdvectionFieldFromList(std::string path);
//     Vector3d getField(Vector3d &pos) const;
//     void loadData(std::string path);    
// };



/*class SourceSpiralArm: public SourceFeature {
private:
	double r3; // parameter that allow a cusp in the inner Galaxy
	double sigmaRinner, sigmaRouter, sigmaPhi, sigmaZ; // parameter for scalelength of the spiral arm
	std::vector<double> aI, alphaI; // parameter for each spiral arm;

public:
	SourceSpiralArm();
	void prepareParticle(ParticleState &particle) const;
	double phiArm(double r, int i);
	double sourceDensity(double r, double phi, double z);
};*/

// -----------------------------------------------------------------------------------------------------------------------

class Array4D: public Referenced {
public:
	Array4D(int nX, int nY, int nZ, int nE);
	int nX, nY, nZ, nE;
	std::vector< std::vector< std::vector< std::vector<double>>>> data;
	void fillGrid();
	void fillGrid(int nX, int nY, int nZ, int nE);
	void setValue(double value, int iX, int iY, int iZ, int iE);
	std::vector<double> xGrid, yGrid, zGrid, eGrid;
	void printHistSize();
};

class ObserverHistogram4D: public Module {
public:
	
	std::vector<double> xGrid, yGrid, zGrid, eGrid;
	ref_ptr<Array4D> hist;
	double weightSum = 0;
	void getIndexFromValue(double value, std::vector<double> grid, int &i) const;

	ObserverHistogram4D();
	void saveHistogram(std::string path);
	void fillCandidateInHistogram(const double x, const double y, const double z, const double energy, const double w, Array4D &data) const;
	void process(Candidate *cand) const;
	
	void setXGrid(double min, double max, int nX);
	void setYGrid(double min, double max, int nY);
	void setZGrid(double min, double max, int nZ);
	void setEGrid(double min, double max, int nE, bool scaleLog = true);
};


// class ObserverHistogram: public Module {
// private:
// 	int nSide, npface, ncap, npix;
// 	std::string folder;
// 	std::string saveName;
// 	bool onEarth = false;
// 	std::vector<std::ostream* > file_outs;

// public:
// 	ObserverHistogram(int nSide, std::string folder, std::string saveName);
	
// 	void process(Candidate* cand) const;

// 	void initFiles();
// 	void setNSide(int n);

// 	int getNSide() const;
// 	int getBinNumber(Vector3d &position) const;

// 	void closeAll();
// };


#endif // CRPROPA_PLUGIN_H