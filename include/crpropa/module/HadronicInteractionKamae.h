#ifndef CRPROPA_HADRONIC_KAMAE_H
#define CRPROPA_HADRONIC_KAMAE_H

#include <crpropa/Module.h>
#include <crpropa/Units.h>
#include <crpropa/massDistribution/Massdistribution.h>

using namespace crpropa;

class HadronicInteractionKamae: public Module {
private:
    ref_ptr<Density> massDensity;
	std::vector<double> a, b, c, d, e, f, g; // parameter for crossection
	double eq6(double x, std::vector<double> a) const;
	double eq9(double x, std::vector<double> b) const;
	double eq12(double x, std::vector<double> c) const;
    double limit;
	
public:
    HadronicInteractionKamae(ref_ptr<Density> density, double limit = 0.01);
	void process(Candidate* cand) const;
    void performInteraction(Candidate* cand) const;

    double crossection(double Tp) const;
    
    // table 3
    double SpectrumPhoton(double Esec, double Tp) const;
    double PhotonND(double Esec, double Tp) const;
    double PhotonDiff(double Esec, double Tp) const;
    double Photon1232(double Esec, double Tp) const;
    double Photon1600(double Esec, double Tp) const;

    // table 4
    double SpectrumElectron(double Esec, double Tp) const;
    double ElectronND(double Esec, double Tp) const;
    double ElectronDiff(double Esec, double Tp) const;
    double Electron1232(double Esec, double Tp) const;
    double Electron1600(double Esec, double Tp) const;

    // table 5
    double SpectrumPositron(double Esec, double Tp) const;
    double PositronND(double Esec, double Tp) const;
    double PositronDiff(double Esec, double Tp) const;
    double Positron1232(double Esec, double Tp) const;
    double Positron1600(double Esec, double Tp) const;

    // table 6
    double SpectrumElectronNeutrino(double Esec, double Tp) const;
    double ElectronNeutrinoND(double Esec, double Tp) const;
    double ElectronNeutrinoDiff(double Esec, double Tp) const;
    double ElectronNeutrino1232(double Esec, double Tp) const;
    double ElectronNeutrino1600(double Esec, double Tp) const;

    // // table 7
    double SpectrumAntiElectronNeutrino(double Esec, double Tp) const;
    double AntiElectronNeutrinoND(double Esec, double Tp) const;
    double AntiElectronNeutrinoDiff(double Esec, double Tp) const;
    double AntiElectronNeutrino1232(double Esec, double Tp) const;
    double AntiElectronNeutrino1600(double Esec, double Tp) const;

    // // table 8
    double SpectrumMuonNeutrino(double Esec, double Tp) const;
    double MuonNeutrinoND(double Esec, double Tp) const;
    double MuonNeutrinoDiff(double Esec, double Tp) const;
    double MuonNeutrino1232(double Esec, double Tp) const;
    double MuonNeutrino1600(double Esec, double Tp) const;
    
    // // table 9
    double SpectrumAntiMuonNeutrino(double Esec, double Tp) const;
    double AntiMuonNeutrinoND(double Esec, double Tp) const;
    double AntiMuonNeutrinoDiff(double Esec, double Tp) const;
    double AntiMuonNeutrino1232(double Esec, double Tp) const;
    double AntiMuonNeutrino1600(double Esec, double Tp) const;
};

#endif // CRPROPA_HADRONIC_KAMAE_H
