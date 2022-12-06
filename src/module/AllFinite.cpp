#include <crpropa/module/AllFinite.h>
#include <kiss/logger.h>

using namespace crpropa;

void CheckNaNModule::process(Candidate *cand) const {
    bool isFinit = true;
    
    // check position information
    isFinit = isFinit & std::isfinite(cand -> current.getPosition().getR2());

    // check direction information
    isFinit = isFinit & std::isfinite(cand -> current.getDirection().getR2());

    // check energy
    isFinit = isFinit & std::isfinite(cand -> current.getEnergy());

    // check Trajectory Lenght
    isFinit = isFinit & std::isfinite(cand -> getTrajectoryLength());

    // // check Time
    // isFinit = isFinit & std::isfinite(cand -> getTime());
    
    if(!isFinit) {
        KISS_LOG_WARNING << "deactivated candidate not finite vaule. \n " << cand -> getDescription() << "\n";
        cand -> setActive(false);
    }

    if(cand -> current.getEnergy() <= 0) {
        KISS_LOG_WARNING << "candidate with Energy 0 or lower is found. Candidate will be deactivated. \n" << cand -> getDescription() << "\n";
        cand -> setActive(false);
    }
}