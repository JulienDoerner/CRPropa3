#include <crpropa/massDistribution/Density.h>
#include <crpropa/Module.h>

using namespace crpropa;

class PiZeroLoss : public Module {
  private: 
	ref_ptr<Density> density;
  public:
	PiZeroLoss(ref_ptr<Density> dens);
	void process(Candidate *candidate) const;
};