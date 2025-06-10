#ifndef CRPROPA_COLUMNDENSITY_H
#define CRPROPA_COLUMNDENSITY_H

#include <crpropa/Module.h>
#include <crpropa/massDistribution/Density.h>

using namespace crpropa;

/*
    @class ColumnDensity
    @brief Module to calculate the column density of a particle. 

    The calculated column density is stored as a property of the Candidate.
*/
class ColumnDensity : public Module {
  private: 
    std::string key = "CD"; //< key for the property
    ref_ptr<Density> density; //< reference to the density object

  public: 
    ColumnDensity(ref_ptr<Density> density);

    void process(Candidate *candidate) const;

    void setKey(std::string key);
    std::string getKey() const;
};

#endif // CRPROPA_COLUMNDENSITY_H