#include "crpropa/massDistribution/CMZDensity.h"
#include "crpropa/Common.h"
namespace crpropa {
    
CMZDensity::CMZDensity(bool clouds, bool diffuse) {
    useClouds = clouds; 
    useDiffuse = diffuse;
    
    // add positions for molecular clouds 
    MC_orig.push_back(Vector3d(0, 81.6, -13.32) * pc);     // SgrC
    MC_orig.push_back(Vector3d(0, 19.29, -11.87) * pc);    // 20 km/s
    MC_orig.push_back(Vector3d(0, 2.97, -10.38) * pc);     // 50 km/s
    MC_orig.push_back(Vector3d(0, -37.53, 2.37) * pc);     // Ridge A
    MC_orig.push_back(Vector3d(0, -50.44, 8.16) * pc);     // Ridge B
    MC_orig.push_back(Vector3d(0, -56.37, 7.42) * pc);     // Ridge C
    MC_orig.push_back(Vector3d(0, -61.12, 7.71) * pc);     // Ridge D
    MC_orig.push_back(Vector3d(0, -70.91, 0.74) * pc);     // Ridge E
    MC_orig.push_back(Vector3d(0, -73.58, 2.97) * pc);     // Ridge F
    MC_orig.push_back(Vector3d(0, -97.92, 5.93) * pc);     // SgrB2
    MC_orig.push_back(Vector3d(0, -167.66, -16.32) * pc);  // SgrD
}

double CMZDensity::cloudDensity(const Vector3d &pos, const Vector3d &orig, double R, double n) const {
    // check if position is inside of the cloud 
    bool inside = (pos - orig).getR() < R;
    return inside ? n : 0; // return n inside and 0 outside
}

double CMZDensity::getDensity(const Vector3d &position) const {
    return getH2Density(position) + getHIDensity(position);
}

double CMZDensity::getNucleonDensity(const Vector3d &position) const {
    return getHIDensity(position) +  2 * getH2Density(position);
}

double CMZDensity::getH2Density(const Vector3d &position) const {
        //VALUES FOR n PER M^3
        double n=0;

        if (useClouds) {
            for(int i = 0; i < MC_density.size(); i++)
                n += cloudDensity(position, MC_orig[i], MC_radius[i] * pc, MC_density[i]);
        }
        
        if (useDiffuse) {
            // intercloud medium in CMZ from Ferriere 2007
            double x = - position.x;
            double y = - position.y;
            double z = position.z;
            double xC = - 50*pc;        //offset
            double yC = 50*pc;
            double ThettaC = 70./180. * M_PI; //radiant
            double xX = (x - xC)*cos(ThettaC) + (y -yC)*sin(ThettaC);
            double yY = -(x -xC)*sin(ThettaC) + (y - yC)*cos(ThettaC);
            
            double X_C = 125 * pc;
            double Lc = 137 * pc;
            double Hc = 30 * pc;
            
            double arg1 = pow_integer<4>( (std::sqrt(xX * xX + (2.5*yY * 2.5 * yY))  - X_C ) / Lc);
            double arg2 = pow_integer<2>( z / Hc);

            n += 150 / ccm * std::exp(- arg1 - arg2);  
        }
        
        return n;
}


double CMZDensity::getHIDensity(const Vector3d &position) const {
        //VALUES FOR n PER M^3
        double n=0;
        double pi=3.1415926535;
        double x = - position.x;    // - because CRPropas framework is toward the Galactic Center
        double y = - position.y;
        double z = position.z;
        double cmm=1e6; 
        
        // intercloud medium in CMZ from Ferriere 2007
        double xC = -50*pc;        //offset
        double yC = 50*pc;
        
        double ThettaC = 70./180*M_PI; //radiant

        double xX = (x - xC)*cos(ThettaC) + (y -yC)*sin(ThettaC);
        double yY = -(x -xC)*sin(ThettaC) + (y - yC)*cos(ThettaC);
        double zZ=z;

        double X_C = 125 * pc;
        double Lc = 137 * pc;
        double Hc = 54 * pc;
        double arg1 = pow_integer<4>( (std::sqrt(xX * xX + 2.5*yY * 2.5 * yY)  - X_C ) / Lc);
        double arg2 = pow_integer<2>( z / Hc);
        return 8.8 * cmm * std::exp(-arg1) * std::exp(-arg2);

        n+=8.8 * std::exp(-pow((std::sqrt(xX*xX+(2.5*yY)*(2.5*yY))-125.*pc)*1./(137.*pc), 4.)) * std::exp(-pow((zZ/(54.*pc)),2.)) * cmm ;
   
        return n;
}

    

    
} //END NAMESPACE CRPROPA

