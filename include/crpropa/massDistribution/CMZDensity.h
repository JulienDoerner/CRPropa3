#ifndef CRPROPA_CMZDENSITY_H
#define CRPROPA_CMZDENSITY_H

#include "crpropa/massDistribution/Density.h"
#include "crpropa/Vector3.h"
#include "crpropa/Units.h"

#include <math.h>

namespace crpropa {
    /*
     ..
     */
    
    class CMZDensity: public Density {
    private: 
        bool useClouds, useDiffuse;

        // parameters for the molecular clouds
        std::vector<double> MC_density = {1.8e10, 1e10, 1e10, 1.3e10, 7e9, 9e9, 8e9, 2.8e10, 3.2e10, 3e9, 6e9 };
        std::vector<double> MC_radius = {1.7, 9.4, 4.5, 2.4, 1.9, 1.9, 3.3, 3.5, 2.4, 10, 1.8}; // in pc
        std::vector<Vector3d> MC_orig;

    public:

        double cloudDensity(const Vector3d &position, const Vector3d &center, double radius, double density) const;

        CMZDensity(bool useClouds = true, bool useDiffuse = true);
        double getDensity(const Vector3d &position) const;
        double getHIDensity(const Vector3d &position) const;
        double getH2Density(const Vector3d &position) const;
        double getNucleonDensity(const Vector3d &position) const;
        
        void printPositions() {
            for(int i = 0; i < MC_orig.size(); i++) {
                std::cout << "number i = "<< i << " at position " << MC_orig[i] << "\n";
            }
        }
    };
    
} // namespace

#endif //CRPROPA_CMZDENSITY_H

