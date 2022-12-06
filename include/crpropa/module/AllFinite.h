#ifndef CRPROPA_ALLFINITE_H
#define CRPROPA_ALLFINITE_H

#include <crpropa/Module.h>
#include <crpropa/Vector3.h>
#include <crpropa/ParticleState.h>
#include <crpropa/Source.h>
#include <cmath>

using namespace crpropa;

class CheckNaNModule: public Module {
    public:
    void process(Candidate *cand) const;
};

class SourceEnergyNaN: public SourceFeature {
    public:
    void prepareParticle(ParticleState &particle) const {
        particle.setEnergy(NAN);
    }
};

#endif // CRPROPA_ALLFINITE_H