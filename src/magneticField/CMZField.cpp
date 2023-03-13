#include "crpropa/magneticField/CMZField.h"
#include "crpropa/Units.h"

namespace crpropa {

CMZField::CMZField() {
    useMCField = false;
    useICField = true;
    useNTFField = false;
    useRadioArc = false;
}

double CMZField::getA(double a1) const {
    return 4*log(2)/a1/a1;  
}

double CMZField::getL(double a2) const {
    return a2/2*log(2);
}

Vector3d CMZField::BPol(const Vector3d& position,const Vector3d& mid, double B1, double a, double L) const{
    // cylindircal coordinates
    Vector3d pos = position - mid;
    double r = sqrt(pos.x*pos.x + pos.y*pos.y);
    double phi = std::atan2(pos.y, pos.x);

    double r1 = 1/(1+a*pos.z*pos.z);
    double Bs = B1*exp(-r1*r/L);
    double Br = 2*a*pow_integer<3>(r1)*r*pos.z*Bs;
    
    Vector3d b = Vector3d(0.);
    b.z = r1*r1*Bs;
    // transform to cartesian coordinates
    b.x = Br*cos(phi);
    b.y = Br*sin(phi);
    return b;
}

Vector3d CMZField::BAz(const Vector3d& position, const Vector3d& mid, double B1, double eta, double R) const {
    // cylindrical coordinates
    Vector3d pos = position - mid;
    double r = sqrt(pos.x*pos.x + pos.y*pos.y);
    double phi = std::atan2(pos.y,pos.x);

    Vector3d bVec(0.);
    double Hc = R/sqrt(log(2));
    double b = 1.;
    double m = 1;
    double r1 = R/10;
    double v = m/eta*log((r+b)/(R+b));
    double cosV = cos(v + m*phi);

    double Br=0;
    double Bphi=0;
    
    if(r>r1){
        double Pre = B1*cosV*exp(-pos.z*pos.z/Hc/Hc);
        Br = Pre*R/r;
        Bphi=-Pre/eta*R/(r+b);
    }
    else{
        double Pre = B1*exp(-pos.z*pos.z/Hc/Hc)*R/r1*(3*r/r1 - 2*r*r/r1/r1)*cosV;
        Br = Pre;
        Bphi = 1 + 6*(r-r1)/(2*r-3*r1)*(sin(v+m*phi)-sin(v))/cosV;
        Bphi *= -Pre*r/eta/(r+b);
    }

    bVec.x = Br*cos(phi) - Bphi*sin(phi);
    bVec.y = Br*sin(phi) + Bphi*cos(phi);

    return bVec;
}

bool CMZField::getUseMCField() const {
    return useMCField;
}
bool CMZField::getUseICField() const {
    return useICField;
}
bool CMZField::getUseNTFField() const {
    return useNTFField;
}
bool CMZField::getUseRadioArc() const {
    return useRadioArc;
}

void CMZField::setUseMCField(bool use) {
    useMCField = use;
}
void CMZField::setUseICField(bool use) {
    useICField = use;
}
void CMZField::setUseNTFField(bool use) {
    useNTFField = use;
}
void CMZField::setUseRadioArc(bool use) {
    useRadioArc = use;
}

Vector3d CMZField::getMCField(const Vector3d& pos) const {//Field in molecular clouds
    Vector3d b(0.); // calculation in milli Gauss
    double eta=0.01;
    double N=59; // normalization factor, depends on eta

	// azimuthal component in dense clouds
    //A=SgrC 
    Vector3d mid(0, 81.60 * pc, -16.32 * pc);
    double R = 1.7*pc; 
    double B1 = 0.6;  
    b += BAz(pos, mid, B1 / N, eta, R);

    // A=G0.253+0.016 Dust Ridge A
    mid = Vector3d(0, - 37.53 * pc, 2.37 * pc);
    R = 2.4*pc; 
    B1 = 0.5;
    b += BAz(pos, mid, B1 / N, eta, R);

    //A=Dust Ridge B
    mid = Vector3d(0, -50.44, 8.16) * pc;
    R = 1.9 * pc; 
    B1 = 1.3;
    b += BAz(pos, mid, B1 / N, eta, R);
    
    //A=Dust Ridge C
    mid = Vector3d(0, -56.37, 7.42) * pc;
    R = 1.9 * pc; 
    B1 = 0.4;
    b += BAz(pos, mid, B1 / N, eta, R);

    //A=Dust Ridge D
    mid = Vector3d(0, -61.12, 7.71) * pc;
    R = 3.3 * pc; 
    B1 = 0.4;
    b += BAz(pos, mid, B1 / N, eta, R);
    
    //A=Dust Ridge E
    mid = Vector3d(0, -70.91, 0.74) * pc;
    R = 3.5 * pc; 
    B1 = 0.7;
    b += BAz(pos, mid, B1 / N, eta, R);
        
    //A=Dust Ridge F
    mid = Vector3d(0, -73.58, 2.97) * pc;
    R = 2.4 * pc; 
    B1 = 0.9;
    b += BAz(pos, mid, B1 / N, eta, R);

    //Sgr D
    mid = Vector3d(0, -167.66, -16.32) * pc;
    R = 1.8 * pc;
    B1= 0.2;
    b += BAz(pos, mid, B1 / N, eta, R);

    //Sgr B2
    mid = Vector3d(0, -97.92, -5.93) * pc;
    R = 10 * pc; 
    B1 = 1.3;
    b += BAz(pos, mid, B1 / N, eta, R);
        
    //A=Inner R=5pc
    mid = Vector3d(0, 8.31, -6.82) * pc;
    R = 5 * pc; 
    B1 = 8.4;
    b += BAz(pos, mid, B1 / N, eta, R);

    //20 km s^-1
    mid = Vector3d(0, 19.29, -11.87) * pc;
    R = 9.4 * pc; 
    B1 = 0.9;
    b += BAz(pos, mid, B1 / N, eta, R);    
    
    //50 km s^-1
    mid = Vector3d(0, 2.97, -10.38) * pc;
    R = 4.5 * pc;
    B1= 1.3;
    b += BAz(pos, mid, B1 / N, eta, R);    
    
    //SgrA* is different orrientated! 
    //only phi component
    double x = pos.x;
    double y = pos.y - 8.31 * pc;
    double z = pos.z + 6.82 * pc;
    R=1.2e12 * meter; 
    B1=65./3.07;
    double Hc = R/sqrt(log(2));
    double r = sqrt(x*x + y*y);
    double r1 = R/10;
    double phi= std::atan2(y,x);
    double Bphi;

    if(r>r1){
        Bphi = - B1*exp(-z*z/Hc/Hc)*R/r;
    }
    else{
        - B1*exp(-z*z/Hc/Hc)*R/r1*(3*r/r1- 2*r*r/r1*r1);
    }

    b.x -= Bphi*sin(phi);
    b.y += Bphi*cos(phi);

    return b * milli * gauss;
} 

Vector3d CMZField::getICField(const Vector3d& pos) const {//Field in intercloud medium--> poloidal field
    Vector3d mid(0., 8.31 * pc, -6.82 * pc); // centered in SgrA*

    double eta = 0.85;
    double B1 = 1e-5*gauss;
    double B2 = B1/eta;
    double a = 4*log(2)/pow(70*pc, 2); 
    double L = 158*pc/log(2);

    return BPol(pos, mid, B2, a, L);
}
  
Vector3d CMZField::getNTFField(const Vector3d& pos) const {//Field in the non-thermal filaments--> predominantly poloidal field (except "pelical"-> azimuthal)
    Vector3d b(0.); // units in muG
    Vector3d mid(0.);

    //A=SgrC
    mid = Vector3d(0., 81.60, -1.48) * pc;
    double a1=27.4 * pc;
    double a2=1.7 * pc;
    double eta=0.48;
    double B1=100;
    b += BPol(pos, mid, B1/eta, getA(a1), getL(a2));

    //A=G359.15-0.2 The Snake
    mid = Vector3d(0, 126.11, -25.22) * pc;
    a1=12.9 * pc;
    a2=1.5 * pc;
    B1=88;
    b += BPol(pos, mid, B1/eta, getA(a1), getL(a2));

    //A=G359.54+0.18 Nonthermal Filament (Ripple)
    mid = Vector3d(0, 66.76, 25.22) * pc;
    a1=15.0 * pc;
    a2=2.7 * pc;
    B1=1000;
    b += BPol(pos, mid, B1/eta, getA(a1), getL(a2));

    //A=G359.79 +17 Nonthermal Filament
    mid = Vector3d(0, 29.67, 23.74) * pc;
    a1=16.2 * pc;
    a2=3.5 * pc;
    B1=1000;
    b += BPol(pos, mid, B1/eta, getA(a1), getL(a2));

    //A=G359.96 +0.09  Nonthermal Filament Southern Thread
    mid = Vector3d(0, 5.93, 16.32) * pc;
    a1=28.7 * pc;
    a2=1.7 * pc;
    B1=100;
    b += BPol(pos, mid, B1/eta, getA(a1), getL(a2));

    //A=G0.09 +0.17  Nonthermal Filament Northern thread
    mid = Vector3d(0, -13.35, 25.22) * pc;
    a1=29.4 * pc;
    a2=2.0 * pc;
    B1=140;
    b += BPol(pos, mid, B1/eta, getA(a1), getL(a2));

    //A=G359.85+0.47  Nonthermal Filament The Pelican  is not poloidal but azimuthal
    mid = Vector3d(0, 22.25, 69.73) * pc;
    a1=13.4 * pc;
    a2=2.2 * pc;
    B1=70; 
    // by and bz switched because pelican is differently oriented
    Vector3d bPelican = BPol(pos, mid, B1/eta, getA(a1), getL(a2));
    b.x += bPelican.x;
    b.y += bPelican.z;
    b.z += bPelican.y;
    
	return b*muG;
}
  
Vector3d CMZField::getRadioArcField(const Vector3d& pos) const {//Field in the non-thermal filaments--> predominantly poloidal field
    //poloidal field in the non-thermal filament region A=RadioArc
    double eta=0.48;
    Vector3d mid(0, -26.70 * pc, -10.38 * pc);
    double a1=70.5 * pc;
    double a2=9.9 * pc;
    double B1=1000;
    return BPol(pos, mid, B1/eta, getA(a1), getL(a2)) * muG;
}

Vector3d CMZField::getField(const Vector3d& pos) const{
    Vector3d b(0.);

    if(useMCField){
        b += getMCField(pos);
    }
    if(useICField){
        b += getICField(pos);
    }
    if(useNTFField){
        b += getNTFField(pos);
    }
    if(useRadioArc){
        b += getRadioArcField(pos);
    }

    return b;
}

} // namespace crpropa
