#ifndef CRPROPA_KISSMANN_H
#define CRPROPA_KISSMANN_H

#include <crpropa/Module.h>
#include <crpropa/Source.h>
#include <crpropa/Random.h>
#include <crpropa/Units.h>
#include <crpropa/Vector3.h>

#include <cmath>

using namespace crpropa;

class SourceSNRKissmann: public SourceFeature {
private:
	double R_earth = 8.5 * kpc;
	double R_max = 20 * kpc;
	double alpha = 0.475063;
    double beta = 2.16570;
	double frMax = pow(alpha/beta, alpha)*exp(beta - alpha);
public:
	SourceSNRKissmann();
	double f_r(double r) const;
	void prepareParticle(ParticleState& particle) const;
};

class SourceSpiralArm: public SourceFeature {
private:
    double alpha[4] = {0.242, 0.279, 0.249, 0.240};
    double a[4] = {0.246 * kpc, 0.608 * kpc, 0.449 * kpc, 0.378 * kpc};
    double r3 = 2.9 * kpc;
    double sRouter = 3.1 * kpc;
    double sRinner = 0.7 * kpc;
    double sPhi = 15 * deg;
    double sZ = 70 * pc;

    double rMax = 20 * kpc;
    double zMax = 2 * kpc;
public:
    SourceSpiralArm();
    double sourceDensity(Vector3d pos) const;
    void prepareParticle(ParticleState &particle) const;

    void setRMax(double r);
    void setZMax(double z);

    double getRMax();
    double getZMax();
};

class CheckNaNModule: public Module {
    public:
    void process(Candidate *cand) const;
};

class SourceEnergyNaN: public SourceFeature{
    public:
    void prepareParticle(ParticleState &particle) const {
        particle.setEnergy(NAN);
    }
};

#endif // CRPROPA_KISSMANN_H