#include "crpropa/PhotonBackground.h"
#include "crpropa/Units.h"
#include "crpropa/Random.h"

#include <fstream>
#include <stdexcept>
#include <limits>
#include <cmath>

namespace crpropa {

// SpaceDependendPhotonField::SpaceDependendPhotonField(std::string fieldName) : PhotonField() {
// 	this->fieldName = fieldName;
// 	this->isRedshiftDependent = false;
// 	// this->isSpaceDependend = true;
// 	// this -> grid = grid;

// 	readPhotonEnergy(getDataPath("") + "Scaling/" + this->fieldName + "_photonEnergy.txt");
// 	readPhotonDensity(getDataPath("") + "Scaling/" + this->fieldName + "_photonDensity.txt");
// }

// double SpaceDependendPhotonField::getPhotonDensity(double eps) const {
// 	return interpolate(eps, this->photonEnergies, this->photonDensity);
// }


// double SpaceDependendPhotonField::getMinimumPhotonEnergy(double z) const{
// 	return photonEnergies[0];
// }

// double SpaceDependendPhotonField::getMaximumPhotonEnergy(double z) const{
// 	return photonEnergies[photonEnergies.size() -1];
// }

// void SpaceDependendPhotonField::readPhotonEnergy(std::string filePath) {
// 	std::ifstream infile(filePath.c_str());
// 	if (!infile.good())
// 		throw std::runtime_error("SpaceDependendPhotonField::readPhotonEnergy: could not open " + filePath);

// 	std::string line;
// 	while (std::getline(infile, line)) {
// 		if ((line.size() > 0) & (line[0] != '#') )
// 			this->photonEnergies.push_back(std::stod(line));
// 	}
// 	infile.close();
// }

// void SpaceDependendPhotonField::readPhotonDensity(std::string filePath) {
// 	std::ifstream infile(filePath.c_str());
// 	if (!infile.good())
// 		throw std::runtime_error("SpaceDependendPhotonField::readPhotonDensity: could not open " + filePath);

// 	std::string line;
// 	while (std::getline(infile, line)) {
// 		if ((line.size() > 0) & (line[0] != '#') )
// 			this->photonDensity.push_back(std::stod(line));
// 	}
// 	infile.close();
// }

// double SpaceDependendPhotonField::getSpaceScale(Vector3d &pos) const {
// 	return grid -> interpolate(pos);
// }

// void SpaceDependendPhotonField::setSpaceGrid(ref_ptr<Grid1f> grid) {
// 	this->grid = grid;
// }

// --------------------------------------------------------------------------------------------------------

TabularPhotonField::TabularPhotonField(std::string fieldName, bool isRedshiftDependent) : PhotonField() {
	this->fieldName = fieldName;
	this->isRedshiftDependent = isRedshiftDependent;

	initAll();
}

void TabularPhotonField::initAll() {
	readPhotonEnergy(getDataPath("") + "Scaling/" + this->fieldName + "_photonEnergy.txt");
	readPhotonDensity(getDataPath("") + "Scaling/" + this->fieldName + "_photonDensity.txt");
	if (this->isRedshiftDependent)
		readRedshift(getDataPath("") + "Scaling/" + this->fieldName + "_redshift.txt");

	checkInputData();

	if (this->isRedshiftDependent)
		initRedshiftScaling();	
}


double TabularPhotonField::getPhotonDensity(double Ephoton, double z) const {
	if (this->isRedshiftDependent) {
		return interpolate2d(Ephoton, z, this->photonEnergies, this->redshifts, this->photonDensity);
	} else {
		return interpolate(Ephoton, this->photonEnergies, this->photonDensity);
	}
}


double TabularPhotonField::getRedshiftScaling(double z) const {
	if (!this->isRedshiftDependent)
		return 1.;
 
	if (z < this->redshifts.front())
		return 1.;
 
	if (z > this->redshifts.back())
		return 0.;
 
	return interpolate(z, this->redshifts, this->redshiftScalings);
}

double TabularPhotonField::getMinimumPhotonEnergy(double z) const{
	return photonEnergies[0];
}

double TabularPhotonField::getMaximumPhotonEnergy(double z) const{
	return photonEnergies[photonEnergies.size() -1];
}

void TabularPhotonField::readPhotonEnergy(std::string filePath) {
	std::ifstream infile(filePath.c_str());
	if (!infile.good())
		throw std::runtime_error("TabularPhotonField::readPhotonEnergy: could not open " + filePath);

	std::string line;
	while (std::getline(infile, line)) {
		if ((line.size() > 0) & (line[0] != '#') )
			this->photonEnergies.push_back(std::stod(line));
	}
	infile.close();
}

void TabularPhotonField::readPhotonDensity(std::string filePath) {
	std::ifstream infile(filePath.c_str());
	if (!infile.good())
		throw std::runtime_error("TabularPhotonField::readPhotonDensity: could not open " + filePath);

	std::string line;
	while (std::getline(infile, line)) {
		if ((line.size() > 0) & (line[0] != '#') )
			this->photonDensity.push_back(std::stod(line));
	}
	infile.close();
}

void TabularPhotonField::readRedshift(std::string filePath) {
	std::ifstream infile(filePath.c_str());
	if (!infile.good())
		throw std::runtime_error("TabularPhotonField::initRedshift: could not open " + filePath);

	std::string line;
	while (std::getline(infile, line)) {
		if ((line.size() > 0) & (line[0] != '#') )
			this->redshifts.push_back(std::stod(line));
	}
	infile.close();
}

void TabularPhotonField::initRedshiftScaling() {
	double n0 = 0.;
	for (int i = 0; i < this->redshifts.size(); ++i) {
		double z = this->redshifts[i];
		double n = 0.;
		for (int j = 0; j < this->photonEnergies.size()-1; ++j) {
			double e_j = this->photonEnergies[j];
			double e_j1 = this->photonEnergies[j+1];
			double deltaLogE = std::log10(e_j1) - std::log10(e_j);
			if (z == 0.)
				n0 += (getPhotonDensity(e_j, 0) + getPhotonDensity(e_j1, 0)) / 2. * deltaLogE;
			n += (getPhotonDensity(e_j, z) + getPhotonDensity(e_j1, z)) / 2. * deltaLogE;
		}
		this->redshiftScalings.push_back(n / n0);
	}
}

void TabularPhotonField::checkInputData() const {
	if (this->isRedshiftDependent) {
		if (this->photonDensity.size() != this->photonEnergies.size() * this-> redshifts.size())
			throw std::runtime_error("TabularPhotonField::checkInputData: length of photon density input is unequal to length of photon energy input times length of redshift input");
	} else {
		if (this->photonEnergies.size() != this->photonDensity.size())
			throw std::runtime_error("TabularPhotonField::checkInputData: length of photon energy input is unequal to length of photon density input");
	}

	for (int i = 0; i < this->photonEnergies.size(); ++i) {
		double ePrevious = 0.;
		double e = this->photonEnergies[i];
		if (e <= 0.)
			throw std::runtime_error("TabularPhotonField::checkInputData: a value in the photon energy input is not positive");
		if (e <= ePrevious)
			throw std::runtime_error("TabularPhotonField::checkInputData: photon energy values are not strictly increasing");
		ePrevious = e;
	}

	for (int i = 0; i < this->photonDensity.size(); ++i) {
		if (this->photonDensity[i] < 0.)
			throw std::runtime_error("TabularPhotonField::checkInputData: a value in the photon density input is negative");
	}

	if (this->isRedshiftDependent) {
		if (this->redshifts[0] != 0.)
			throw std::runtime_error("TabularPhotonField::checkInputData: redshift input must start with zero");

		for (int i = 0; i < this->redshifts.size(); ++i) {
			double zPrevious = -1.;
			double z = this->redshifts[i];
			if (z < 0.)
				throw std::runtime_error("TabularPhotonField::checkInputData: a value in the redshift input is negative");
			if (z <= zPrevious)
				throw std::runtime_error("TabularPhotonField::checkInputData: redshift values are not strictly increasing");
			zPrevious = z;
		}

		for (int i = 0; i < this->redshiftScalings.size(); ++i) {
			double scalingFactor = this->redshiftScalings[i];
			if (scalingFactor <= 0.)
				throw std::runtime_error("TabularPhotonField::checkInputData: initRedshiftScaling has created a non-positive scaling factor");
		}
	}
}

// --------------------------------------------------------


SpaceDependendPhotonField::SpaceDependendPhotonField(std::string fieldName, bool isRedshiftDependent) : PhotonField() {
	this->fieldName = fieldName;
	this->isRedshiftDependent = isRedshiftDependent;

	initAll();
}

void SpaceDependendPhotonField::initAll() {
	readPhotonEnergy(getDataPath("") + "Scaling/" + this->fieldName + "_photonEnergy.txt");
	readPhotonDensity(getDataPath("") + "Scaling/" + this->fieldName + "_photonDensity.txt");
	if (this->isRedshiftDependent)
		readRedshift(getDataPath("") + "Scaling/" + this->fieldName + "_redshift.txt");

	checkInputData();

	if (this->isRedshiftDependent)
		initRedshiftScaling();	
}


double SpaceDependendPhotonField::getPhotonDensity(double Ephoton, double z) const {
	if((Ephoton < getMinimumPhotonEnergy(z))||(Ephoton > getMaximumPhotonEnergy(z)))
		return 0;

	if (this->isRedshiftDependent) {
		return interpolate2d(Ephoton, z, this->photonEnergies, this->redshifts, this->photonDensity);
	} else {
		return interpolate(Ephoton, this->photonEnergies, this->photonDensity);
	}
}


double SpaceDependendPhotonField::getRedshiftScaling(double z) const {
	if (!this->isRedshiftDependent)
		return 1.;
 
	if (z < this->redshifts.front())
		return 1.;
 
	if (z > this->redshifts.back())
		return 0.;
 
	return interpolate(z, this->redshifts, this->redshiftScalings);
}

double SpaceDependendPhotonField::getMinimumPhotonEnergy(double z) const{
	return photonEnergies[0];
}

double SpaceDependendPhotonField::getMaximumPhotonEnergy(double z) const{
	return photonEnergies[photonEnergies.size() -1];
}

void SpaceDependendPhotonField::readPhotonEnergy(std::string filePath) {
	std::ifstream infile(filePath.c_str());
	if (!infile.good())
		throw std::runtime_error("SpaceDependendPhotonField::readPhotonEnergy: could not open " + filePath);

	std::string line;
	while (std::getline(infile, line)) {
		if ((line.size() > 0) & (line[0] != '#') )
			this->photonEnergies.push_back(std::stod(line));
	}
	infile.close();
}

void SpaceDependendPhotonField::readPhotonDensity(std::string filePath) {
	std::ifstream infile(filePath.c_str());
	if (!infile.good())
		throw std::runtime_error("SpaceDependendPhotonField::readPhotonDensity: could not open " + filePath);

	std::string line;
	while (std::getline(infile, line)) {
		if ((line.size() > 0) & (line[0] != '#') )
			this->photonDensity.push_back(std::stod(line));
	}
	infile.close();
}

void SpaceDependendPhotonField::readRedshift(std::string filePath) {
	std::ifstream infile(filePath.c_str());
	if (!infile.good())
		throw std::runtime_error("SpaceDependendPhotonField::initRedshift: could not open " + filePath);

	std::string line;
	while (std::getline(infile, line)) {
		if ((line.size() > 0) & (line[0] != '#') )
			this->redshifts.push_back(std::stod(line));
	}
	infile.close();
}

void SpaceDependendPhotonField::initRedshiftScaling() {
	double n0 = 0.;
	for (int i = 0; i < this->redshifts.size(); ++i) {
		double z = this->redshifts[i];
		double n = 0.;
		for (int j = 0; j < this->photonEnergies.size()-1; ++j) {
			double e_j = this->photonEnergies[j];
			double e_j1 = this->photonEnergies[j+1];
			double deltaLogE = std::log10(e_j1) - std::log10(e_j);
			if (z == 0.)
				n0 += (getPhotonDensity(e_j, 0) + getPhotonDensity(e_j1, 0)) / 2. * deltaLogE;
			n += (getPhotonDensity(e_j, z) + getPhotonDensity(e_j1, z)) / 2. * deltaLogE;
		}
		this->redshiftScalings.push_back(n / n0);
	}
}

void SpaceDependendPhotonField::checkInputData() const {
	if (this->isRedshiftDependent) {
		if (this->photonDensity.size() != this->photonEnergies.size() * this-> redshifts.size())
			throw std::runtime_error("SpaceDependendPhotonField::checkInputData: length of photon density input is unequal to length of photon energy input times length of redshift input");
	} else {
		if (this->photonEnergies.size() != this->photonDensity.size())
			throw std::runtime_error("SpaceDependendPhotonField::checkInputData: length of photon energy input is unequal to length of photon density input");
	}

	for (int i = 0; i < this->photonEnergies.size(); ++i) {
		double ePrevious = 0.;
		double e = this->photonEnergies[i];
		if (e <= 0.)
			throw std::runtime_error("SpaceDependendPhotonField::checkInputData: a value in the photon energy input is not positive");
		if (e <= ePrevious)
			throw std::runtime_error("SpaceDependendPhotonField::checkInputData: photon energy values are not strictly increasing");
		ePrevious = e;
	}

	for (int i = 0; i < this->photonDensity.size(); ++i) {
		if (this->photonDensity[i] < 0.)
			throw std::runtime_error("SpaceDependendPhotonField::checkInputData: a value in the photon density input is negative");
	}

	if (this->isRedshiftDependent) {
		if (this->redshifts[0] != 0.)
			throw std::runtime_error("SpaceDependendPhotonField::checkInputData: redshift input must start with zero");

		for (int i = 0; i < this->redshifts.size(); ++i) {
			double zPrevious = -1.;
			double z = this->redshifts[i];
			if (z < 0.)
				throw std::runtime_error("SpaceDependendPhotonField::checkInputData: a value in the redshift input is negative");
			if (z <= zPrevious)
				throw std::runtime_error("SpaceDependendPhotonField::checkInputData: redshift values are not strictly increasing");
			zPrevious = z;
		}

		for (int i = 0; i < this->redshiftScalings.size(); ++i) {
			double scalingFactor = this->redshiftScalings[i];
			if (scalingFactor <= 0.)
				throw std::runtime_error("SpaceDependendPhotonField::checkInputData: initRedshiftScaling has created a non-positive scaling factor");
		}
	}
}

double SpaceDependendPhotonField::getSpaceScale(Vector3d &pos) {
	if(isSpaceDependend) 
		return grid.interpolate(pos);
	return 1.;
}

void SpaceDependendPhotonField::readGrid(std::string filePath, Vector3d &orig, size_t N, double spacing) {
	GridProperties prop(orig, N, spacing);
	this -> grid = Grid1f(prop);
	this -> isSpaceDependend = true;

	std::ifstream infile(filePath.c_str());
	if (!infile.good())
		throw std::runtime_error("SpaceDependendPhotonField::readGrid: could not open " + filePath);

	std::string line;
	int iX = 0, iY = 0, iZ = 0;
	while (std::getline(infile, line)) {
		if ((line.size() > 0) & (line[0] != '#') ){
			if(iZ == N) {
				iY += 1;
				iZ = 0;
			}
			if(iY == N) {
				iX += 1;
				iY = 0;
			}
			if(iX == N) {
				throw std::runtime_error("SpaceDependendPhotonField::readGrid: to many lines in file for grid size.");
			}
			iZ ++;
			grid.setValue(iX, iY, iZ, std::stof(line));
		}
			
	}
	infile.close();
}

void SpaceDependendPhotonField::setGrid(Grid1f grid) {
	this->grid = grid;
}

BlackbodyPhotonField::BlackbodyPhotonField(std::string fieldName, double blackbodyTemperature) : PhotonField() {
	this->fieldName = fieldName;
	this->blackbodyTemperature = blackbodyTemperature;
	this->quantile = 0.0001; // tested to be sufficient, only used for extreme values of primary energy or temperature
	this->isSpaceDependend = false;
}

double BlackbodyPhotonField::getPhotonDensity(double Ephoton, double z) const {
	return 8 * M_PI * pow_integer<3>(Ephoton / (h_planck * c_light)) / std::expm1(Ephoton / (k_boltzmann * this->blackbodyTemperature));
}

double BlackbodyPhotonField::getMinimumPhotonEnergy(double z) const {
	double A;
	int quantile_int = 10000 * quantile;
	switch (quantile_int)
	{
	case 1:	// 0.01 % percentil
		A = 1.093586e-5 * eV / kelvin;
		break;
	case 10:		// 0.1 % percentil
		A = 2.402189e-5 * eV / kelvin;
		break;
	case 100:		// 1 % percentil
		A = 5.417942e-5 * eV / kelvin;
		break;
	default:
		throw std::runtime_error("Quantile not understood. Please use 0.01 (1%), 0.001 (0.1%) or 0.0001 (0.01%) \n");
		break;
	}
	return A * this -> blackbodyTemperature;
}

double BlackbodyPhotonField::getMaximumPhotonEnergy(double z) const {
	double factor = std::max(1., blackbodyTemperature / 2.73);
	return 0.1 * factor * eV; // T dependent scaling, starting at 0.1 eV as suitable for CMB
}

void BlackbodyPhotonField::setQuantile(double q) {
	if(not ((q == 0.0001) or (q == 0.001) or (q == 0.01)))
		throw std::runtime_error("Quantile not understood. Please use 0.01 (1%), 0.001 (0.1%) or 0.0001 (0.01%) \n");
	this -> quantile = q;
}

} // namespace crpropa
