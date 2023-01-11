#include <limits>
#include <cmath>
#include "kiss/logger.h"
#include "crpropa/Geometry.h"

#include <iostream>

namespace crpropa {

// Plane ------------------------------------------------------------------
Plane::Plane(const Vector3d& _x0, const Vector3d& _n) : x0(_x0), n(_n) {
};

Plane::Plane(const Vector3d& _x0, const Vector3d& v1,const Vector3d& v2) : x0(_x0), n(0,0,0) {
	n = v1.cross(v2);
	n /= n.getR();
};

double Plane::distance(const Vector3d &x) const {
	Vector3d dX = x - x0;
	return n.dot(dX);
};

std::string Plane::getDescription() const {
	std::stringstream ss;
	ss << "Plane: " << std::endl
	   << "   x0: " << x0 << std::endl
	   << "    n: " << n << std::endl;
	return ss.str();
};

Vector3d Plane::normal(const Vector3d& point) const {
	return n;
}


// Sphere ------------------------------------------------------------------
Sphere::Sphere(const Vector3d& _center, double _radius) : center(_center), radius(_radius) {
};

double Sphere::distance(const Vector3d &point) const {
	Vector3d dR = point - center;
	return dR.getR() - radius;
}

Vector3d Sphere::normal(const Vector3d& point) const {
	Vector3d d = point - center;
	return d.getUnitVector();
}

std::string Sphere::getDescription() const {
	std::stringstream ss;
	ss << "Sphere: " << std::endl
			<< "   Center: " << center << std::endl
			<< "   Radius: " << radius << std::endl;
	return ss.str();
};


// ParaxialBox -------------------------------------------------------------
ParaxialBox::ParaxialBox(const Vector3d& _corner, const Vector3d& _size) : corner(_corner), size(_size) {
};

double ParaxialBox::distance(const Vector3d &point) const {
	Vector3d X = point - corner - size/2.;

	// inside the cube
	if ((fabs(X.x) <= size.x/2.) and (fabs(X.y) <= size.y/2.) and (fabs(X.z) <= size.z/2.)) { 
		Vector3d Xp = size/2. - X.abs();
		double d = std::min(Xp.x, std::min(Xp.y, Xp.z));

		return -1. * d;
	}

	double a = std::max(0., fabs(X.x) - size.x/2.);
	double b = std::max(0., fabs(X.y) - size.y/2.);
	double c = std::max(0., fabs(X.z) - size.z/2.);

	return sqrt(a*a + b*b +c*c);
}

Vector3d ParaxialBox::normal(const Vector3d& point) const {
	Vector3d d = (point - corner).abs();
	Vector3d d2 = d + size;
	Vector3d n;
	double dmin = std::numeric_limits<double>::infinity();
	if (d.x < dmin) {
		dmin = d.x;
		n = Vector3d(-1, 0, 0);
	}
	if (d.y < dmin) {
		dmin = d.y;
		n = Vector3d(0, -1, 0);
	}
	if (d.z < dmin) {
		dmin = d.z;
		n = Vector3d(0, 0, -1);
	}
	if (d2.x < dmin) {
		dmin = d2.x;
		n = Vector3d(1, 0, 0);
	}
	if (d2.y < dmin) {
		dmin = d2.y;
		n = Vector3d(0, 1, 0);
	}
	if (d2.z < dmin) {
		// dmin = d2.z;
		n = Vector3d(0, 0, 1);
	}

	return n;
}

std::string ParaxialBox::getDescription() const {
	std::stringstream ss;
	ss << "ParaxialBox: " << std::endl
		 << "   corner: " << corner << std::endl
		 << "     size: " << size << std::endl;
	return ss.str();
};


// Cylinder -------------------------------------------------------------------
Cylinder::Cylinder(const Vector3d& _mid, const double _R, const double _h): mid(_mid), R(_R), h(_h) {};

double Cylinder::distance(const Vector3d &point) const {
  Vector3d pos = point-mid;
  double r = sqrt(pos.x*pos.x + pos.y*pos.y);
  
  double dR = r-R;
  double dZ = fabs(pos.z)-h/2;

  if((dR>0)&&(dZ>0)){
    return sqrt(dR*dR + dZ*dZ);
  }
  else{
    return std::max(dR, dZ);
  }
}

Vector3d Cylinder::normal(const Vector3d &point) const {
  Vector3d pos = point - mid;
  double r = sqrt(pos.x*pos.x + pos.y*pos.y);

  double dR = r-R;
  double dZ = fabs(pos.z) - h/2;

  if(fabs(dR)<fabs(dZ)){
    pos.z = 0;
    return pos.getUnitVector()*dR/fabs(dR);
  }
  if(pos.z>h/2){
    return Vector3d(0,0,1);
  }
  if(pos.z < -h/2){
    return Vector3d(0,0,-1);
  }
  if(pos.z>0){
    return Vector3d(0,0,-1);
  }
  if(pos.z<0){
    return Vector3d(0,0,1);
  }
}
std::string Cylinder::getDescription() const
{
	std::stringstream ss;
	ss << "Cylinder: " << std::endl
		 << "   mid: " << mid << std::endl
		 << "   radius: " << R << std::endl
     << "   height: " << h << std::endl;
	return ss.str();
};

// HollowCylinder -------------------------------------------------------------
HollowCylinder::HollowCylinder(const Vector3d& _mid, const double Rmin, const double Rmax, const double h): mid(_mid), Rmin(Rmin), Rmax(Rmax), h(h) {};

double HollowCylinder::distance(const Vector3d& point) const {
  Vector3d pos = point - mid;
  double r = std::sqrt(pos.x*pos.x+pos.y*pos.y);
  double dR1 = Rmin -r;
  double dR2 = r-Rmax;
  double dZ = fabs(pos.z)- h/2;

  // outside Rmax
  if(dR2>0){
    if(dZ<0)
      return dR2;
    else
      return sqrt(dR2*dR2+dZ*dZ);
  }
  
  // inside Rmin
  if(dR1>0) {
    if(dZ<0)
      return dR1;
    else
      return sqrt(dR1*dR1 + dZ*dZ);
  }

  return std::max({dR1, dR2, dZ}); 
}

Vector3d HollowCylinder::normal(const Vector3d& point) const {
  Vector3d pos = point - mid;
  double r = std::sqrt(pos.x*pos.x+pos.y*pos.y);
  double dR1 = Rmin -r;
  double dR2 = r-Rmax;
  double dZ = fabs(pos.z)- h/2;

  if(dZ > std::min(fabs(dR1), fabs(dR2))){
    Vector3d eR = Vector3d(pos.x, pos.y,0).getUnitVector();
    double sign = dR1*dR2 / fabs(dR1*dR2);
    return sign*eR; 
  }
  else
  {
    if(pos.z>h/2)
      return Vector3d(0,0,1);
    if((pos.z<-h/2) or (pos.z>0))
      return Vector3d(0,0,-1);
    else
      return Vector3d(0,0,1);
  } 
}
std::string HollowCylinder::getDescription() const
{
	std::stringstream ss;
	ss << "HollowCylinder: " << std::endl
		 << "   mid: " << mid << std::endl
		 << "   inner radius: " << Rmin << std::endl
     << "   outer radius: " << Rmax << std::endl
     << "   height: " << h << std::endl;
	return ss.str();
};

} // namespace
