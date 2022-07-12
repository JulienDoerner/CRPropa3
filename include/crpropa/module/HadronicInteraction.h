#ifndef CRPROPA_HADRONICINTERACTION_H
#define CRPROPA_HADRONICINTERACTION_H

#include "crpropa/Module.h"
#include "crpropa/Vector3.h"
#include "crpropa/Random.h"
#include "crpropa/massDistribution/Density.h"
#include "crpropa/ParticleMass.h"
namespace crpropa {
/**
 * \addtogroup EnergyLosses
 * @{
 */

/**
 @class HadronicInteraction
 @brief interactions of nuclei with background protons. Work based on Kelner et al. 2006
 */
class HadronicInteraction: public Module {
protected:
	ref_ptr<Density> massDensity;
	double limit;
	double const L_MAX[7] = {0.96, 0.96, 0.94, 0.98, 0.98, 0.94, 0.98};
	double const W_NDL[7] = {15.0, 20.0, 15.0, 15.0, 15.0, 20.0, 15.0};
	double const W_NDH[7] = {44.0, 45.0, 47.0, 42.0, 40.0, 45.0, 40.0};
	double const mProton = nuclearMass(1,1);
	double const mPion = 0.135 * GeV/c_squared;
  	double const c1[9] = { 0.152322227732,0.807220022742,
      	2.005135155619,3.783473973331,6.204956777877,
      	9.372985251688,13.466236911092,18.833597788992,
      	26.374071890927};
	double const d1[9] = {0.336126421798,0.411213980424,
		0.199287525371,0.474605627657e-1,0.5599626610793e-2,
		0.305249767093e-3,0.659212302608e-5,0.411076933035e-7,
		0.329087403035e-10};

	int flag_Function=0; // 0: Kellner, 1: Kamae, 2: Dermer, 3: Kamae+Kellner, 4: Dermer+Kellner, 5: Dermer+Kamae+Kellner
	double EnergySplit = 100*GeV;
	bool havePhotons = true;
	bool havePions = false;
public:
	HadronicInteraction(
		ref_ptr<Density> density,
		double limit = 0.1, int flag = 0);
	void setMassDensity(ref_ptr<Density> density);
	void setLimit(double limit);
	void process(Candidate *candidate) const;
	void performInteraction(Candidate *candidate) const;

	void setFlagFunction(int flag){
		flag_Function = flag;
	}
	void setEnergySplit(double Energy){
		EnergySplit = Energy;
	}
	void setHavePhotons(bool photon) {
		havePhotons = photon;
	}
	void setHavePions(bool pion) {
		havePions = pion;
	}

	// general functions
	double xSection(double ePrimary) const;
	double spectrumPion(double x, double ePrimary) const;
	double spectrumPhoton(double x, double ePrimary) const;

	// model Kellner 06
	double KellnerXSection(double ePrimary) const;
	double KellnerSpectrumPion(double x, double ePrimary) const;
	double KellnerSpectrumPhoton(double x, double ePrimary) const;

	// model Kamae 06
	double KamaeXSection(double ePrimary) const;
	//double KamaeSpectrumPion(double x, double ePrimary) const;
	double KamaeSpectrumPhoton(double Esec, double Tp) const;
	double KamaeGammaND(double Esec, double Tp) const;
	double KamaeGammaDiff(double Esec, double Tp) const;
	double KamaeGamma1232(double Esec, double Tp) const;
	double KamaeGamma1600(double Esec, double Tp) const;

	// model Dermer 86
	double DermerXSection(double ePrimary) const;
	double DermerSpectrumPhoton(double Esec, double Tp) const;
	double pionb_pi0(double proton_gamma, double total_pion_energy, double cos_theta_max) const;
	double bw(double md) const;




	template<typename SpectrumFunction>
	double sampleParticleEnergy(SpectrumFunction&& spectrumFunction, double ePrimary) const {
		const double xMin = 1. / 1000.;
		const double xMax = 1.;
		const double stepSize = 1. / 1000.;

		double Fmax = 0.;
		double x = xMin;
		while (x < xMax) {
			const double F = spectrumFunction(x);
			if (F > Fmax)
				Fmax = F;
			x += stepSize;
		}

		Random &random = Random::instance();
		double F = 0.;
		do {
			x = std::pow(10, -3 * random.rand());
			F = spectrumFunction(x);
		} while (F < random.rand() * Fmax);
		return x * ePrimary;
	}
};

} // namespace crpropa

#endif // CRPROPA_HADRONICINTERACTION_H
