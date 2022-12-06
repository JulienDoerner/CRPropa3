#include <crpropa/Plugin_Kissmann.h>
#include <crpropa/Common.h>

#include <kiss/logger.h>
#include <cmath>
 
using namespace crpropa;
// -----------------------------------------------------------------------------------------------------------------------

SourceSNRKissmann::SourceSNRKissmann() : SourceFeature(){
}

double SourceSNRKissmann::f_r(double r) const{
	if(r>15*kpc)
		return 0;
	if (r > 10 * kpc)
		r = 10 * kpc;
		
	return pow(r / R_earth, alpha) * std::exp(- beta * (r - R_earth) / R_earth); 
}

void SourceSNRKissmann::prepareParticle(ParticleState& particle) const {
  	Random &random = Random::instance();
	double RPos;
	while (true){
		RPos = random.rand()*R_max;
		double fTest = random.rand()*frMax;
		double fR=f_r(RPos);
		if (fTest<=fR) {
			break;
		}
	}
	double phi = random.rand()*2*M_PI;
	Vector3d pos(cos(phi)*RPos, sin(phi)*RPos, 0.);
	particle.setPosition(pos);
}

// -----------------------------------------------------------------------------------------------------------------------


SourceSpiralArm::SourceSpiralArm() : SourceFeature() {
}

double SourceSpiralArm::sourceDensity(Vector3d pos) const {
    double r = std::sqrt(pos.x * pos.x + pos.y * pos.y);
    double phi = std::atan2(pos.y, pos.x); // - 0.5 * M_PI; // convert coordinate system to Steiman model

    if(!std::isfinite(phi)) {
        KISS_LOG_WARNING << "detect phi as not finite. Set to zero \n";
        phi = 0;
    }

    double nSource = 0;
    for(size_t i = 0; i < 4; i++) { // loop over all arms
        double arm = 0;
        if (r < r3) {
            arm = std::exp(- (r3 - r) / sRinner);
        }
        else {
            arm = std::exp( - (r - r3) / sRouter);
        }

        double phiArm = std::log(r / a[i]) / alpha[i];
        double dPhi = fabs(phi - phiArm);
        dPhi = fmod(dPhi, 2 * M_PI);
        
        // check for finit value of dPhi 
        if(!std::isfinite(dPhi)) {
            KISS_LOG_WARNING << "dPhi is not finite, set to zero!\n";
            dPhi = 0;
        }

        while (dPhi > M_PI) {
            dPhi -= 2 * M_PI;
        }

        arm *= std::exp(- pow_integer<2>(dPhi / sPhi));
        arm *= std::exp(- 0.5 * pow_integer<2>(pos.z / sZ));
        nSource += arm;
    }
    return nSource;
}


void SourceSpiralArm::prepareParticle(ParticleState &particle) const {
    Random &random = Random::instance();
    Vector3d pos;
    double counter = 0;
    while (true) {
        double r = random.rand(rMax);
        double phi = random.rand() * 2 * M_PI;
        if(! std::isfinite(phi)) {
            KISS_LOG_WARNING << "detect phi not finite. Set to zero \n";
            phi = 0;
        }
        pos.x = r * std::cos(phi);
        pos.y = r * std::sin(phi);
        pos.z = 2 * (0.5 - random.rand()) * zMax;

        double test = random.rand();
        double crit = sourceDensity(pos);
        if(test <= crit) {
            break;
        }
        counter += 1;
        if (counter > 1e8)
            KISS_LOG_WARNING << "not possible to find a position within 1e10 tries\n";
            counter = 0;
    } 
    // test nan in position
    if(std::isnan(pos.getR())){
        KISS_LOG_WARNING << "detected NaN in a Candidate. set position to Vector3d(0., 0., 0.)\n";
        pos = Vector3d(0.);
    }
    particle.setPosition(pos);
}

void SourceSpiralArm::setRMax(double r) {
    rMax = r;
}

void SourceSpiralArm::setZMax(double z) {
    zMax = z;
}

double SourceSpiralArm::getRMax() {
    return rMax;
}

double SourceSpiralArm::getZMax() {
    return zMax;
}


// -----------------------------------------------------------------------------------------------------------------------

