#include <crpropa/Plugin.h>
#include <crpropa/Vector3.h>
#include <crpropa/Units.h>
#include <crpropa/Common.h>
#include <crpropa/Version.h>

#include <string>
#include <fstream>
#include <sstream>
#include <iostream>

using namespace crpropa;

// ----------------------------------------------------------------------------

SourceRadiusFromSFR::SourceRadiusFromSFR(std::string filename) : SourceFeature() {
	loadData(filename);
}

void SourceRadiusFromSFR::loadData(std::string filename) {
	std::ifstream infile(filename.c_str());

	if(!infile.good())
		throw std::runtime_error("SourceRadiusFromSFR could not read file: " + filename);
	
	std::istream *in;
	std::string line;
	in = &infile;
	double r, rMax, s, sMax, _;

	while (std::getline(*in, line))
	{
		std::stringstream stream(line);
		if (stream.peek() == '#')
			continue;
		
		stream >> r >> s >> _;
		if (s > sMax) 
			sMax = s;
		
		if ((r>rMax) && (s>0))
			rMax = r;

		radius.push_back(r * kpc);
		sfr.push_back(s);
	}
	this -> sfrMax = sMax;
	this -> rMax = rMax * kpc;
	infile.close();
}

void SourceRadiusFromSFR::prepareParticle(ParticleState& state) const {
	Random &random = Random::instance();
	Vector3d pos = state.getPosition();
	double Rpos;
	int counter = 0;
	while (true)
	{	
		Rpos = random.rand() * rMax;
		double sPos = interpolate(Rpos, radius, sfr);
		double sTest = random.rand() * sfrMax;
		if(sTest <= sPos / 1.6){
			break;
		}
		counter++;
	}
	double phi = random.rand() * 2 * M_PI;
	pos.x = cos(phi) * Rpos;
	pos.y = sin(phi) * Rpos;
	state.setPosition(pos);
}

double SourceRadiusFromSFR::getSFRMax() {
	return sfrMax;
}

double SourceRadiusFromSFR::getRMax() {
	return rMax;
}

// ------------------------------------------------------------------------------

SourceZpositionGauss::SourceZpositionGauss(double zMax, double h) : SourceFeature(), zMax(zMax), h(h) {}

void SourceZpositionGauss::prepareParticle(ParticleState& state) const {
	Vector3d pos = state.getPosition();
	double zPos;
	Random &random = Random::instance();
	double fPos, fTest;
	while (true)
	{
		zPos = (random.rand() - 0.5) * 2 * zMax;
		fPos = exp(- zPos * zPos / h / h);
		fTest = random.rand();
		if (fTest <= fPos) {
			break;
		}
	}
	pos.z = zPos;
	state.setPosition(pos);
}

double SourceZpositionGauss::getZMax() const {
	return zMax;
}

double SourceZpositionGauss::getH() const {
	return h;
}

void SourceZpositionGauss::setZMax(double z) {
	this->zMax = z;
}

void SourceZpositionGauss::setH(double h) {
	this->h = h;
}

// -----------------------------------------------------------------------------------------------------------------------

SynchrotronSelfCompton::SynchrotronSelfCompton(ref_ptr<ModulatedTurbulentField> field, double uRad)  : Module(){
	setMagneticField(field);
	setURad(uRad);
}

void SynchrotronSelfCompton::process(Candidate *c) const {
	int id = fabs(c -> current.getId());
	
	if (id != 11)
		return; // only for electrons
	
	double E = c -> current.getEnergy();
	Vector3d pos = c -> current.getPosition();
	double dT = c -> getCurrentTimeStep(); 

	double dEdT = energyLoss(pos, E);
	double dE = dEdT * dT;
	c -> current.setEnergy(E - dE);
}

double SynchrotronSelfCompton::energyLoss(Vector3d pos, double E) const {
	double B = field -> getBrmsAtPosition(pos) / gauss;
	double beta = 8e-17 * (uRad + 6e11 * pow_integer<2>(B) / 8 / M_PI) / GeV * second;
	return beta * E * E;
}

void SynchrotronSelfCompton::setMagneticField(ref_ptr<ModulatedTurbulentField> field) {
	this -> field = field;
}

void SynchrotronSelfCompton::setURad(double u) {
	uRad = u;
}

double SynchrotronSelfCompton::getURad() const {
	return uRad;
}

std::string SynchrotronSelfCompton::getDescription() const {
	std::stringstream ss;
	ss << "continues energy loss due to Synchrotron and Inverse Compton \n";
	ss << "using a photon field with energy density uRad = " << uRad / eV * ccm << "eV / ccm \n"; 
	return ss.str();
}

// -----------------------------------------------------------------------------------------------------------------------

ObserverScaleHeight::ObserverScaleHeight(double Rmax, std::string path) : ObserverFeature(){
	this -> Rmax = Rmax;
	loadData(path);
}

void ObserverScaleHeight::loadData(std::string path) {
	std::ifstream infile(path.c_str());
	if(!infile.good())
		throw std::runtime_error("could not open file: " + path);

	double r, z;
	std::istream *in;
	std::string line;
	in = &infile;
	while(std::getline(*in, line)){
		std::stringstream stream(line);
		if(stream.peek() == '#')
			continue; // coments start with "#"

		stream >> r >> z;
		rBin.push_back(r);
		zBin.push_back(z);
	}
	infile.close();
}

DetectionState ObserverScaleHeight::checkDetection(Candidate *candidate) const {
	Vector3d pos = candidate -> current.getPosition();

	double R = std::sqrt(pos.x * pos.x + pos.y * pos.y);
	if (R > Rmax)
		return DETECTED;

	double zMax = interpolate(R, rBin, zBin);

	if (fabs(pos.z) > zMax)
		return DETECTED;

	return NOTHING;
}

// -----------------------------------------------------------------------------------------------------------------------
/*
AdvectionFieldFromList::AdvectionFieldFromList(std::string path) {
    loadData(path);
}

void AdvectionFieldFromList::loadData(std::string path) {
    std::ifstream infile(path.c_str());

    if(!infile.good())
        throw std::runtime_error("AdvectionFieldFromList could not open " + path);

    std::istream *in;
    std::string line;
    in = &infile;

    double r, v;
    while (std::getline(*in, line))
    {
        std::stringstream stream(line);
        if (stream.peek() == '#')
            continue;
        
        stream >> r >> v;
        radius.push_back(r * kpc);
        velocity.push_back(v);
    }
    infile.close();
}

Vector3d AdvectionFieldFromList::getField(Vector3d &pos) const {
    double r = std::sqrt(pos.x * pos.x + pos.y * pos.y);
    double sgnZ = pos.z / fabs(pos.z);  

    double v = interpolate(r, radius, velocity) * sgnZ;
    return Vector3d(0., 0., v);
}*/



// -----------------------------------------------------------------------------------------------------------------------
/*
SourceSpiralArm::SourceSpiralArm() : SourceFeature() {
	alphaI = {0.242, 0.279, 0.249, 0.240};
	aI = {0.246 * kpc, 0.608 * kpc, 0.449 * kpc, 0.378 * kpc};
	sigmaRinner = 0.7 * kpc;
	sigmaRouter = 3.1 * kpc;
	sigmaPhi = 15 * deg;
	sigmaZ = 70 * pc;
	r3 = 2.9 * kpc;
}


double SourceSpiralArm::sourceDensity(double r, double phi, double z) {
	double result = 0;
	for(size_t i = 0; i < 4; i++){
		double phiI = phiArm(r, i);
		double arm = 1;
		if(r < r3){
			arm *= std::exp(- (r3 - r)/sigmaRinner);
		}
		else{
			arm *= std::exp(- (r - r3)/sigmaRouter);
		}
		arm *= std::exp(- pow_integer<2>((phi - phiI)/sigmaPhi));
		arm *= std::exp(- pow_integer<2>(z / sigmaZ) / 2);
		result += arm;
	}
	return result;
}

double SourceSpiralArm::phiArm(double r, int i) {
	return std::log(r / aI[i]) / alphaI[i];
}
*/

// -----------------------------------------------------------------------------------------------------------------------

Array4D::Array4D(int nX, int nY, int nZ, int nE){
	this -> nX = nX;
	this -> nY = nY;
	this -> nZ = nZ;
	this -> nE = nE;
	fillGrid(); 
	// std::cout << "created a new Array4D \n";
}

void Array4D::printHistSize() {
	std::cout << "level 0:\t size: nX = " << nX << " vector size: " << data.size() << "\n";
	std::cout << "level 1:\t size: nY = " << nY << " vector size: " << data[0].size() << "\n";
	std::cout << "level 2:\t size: nZ = " << nZ << " vector size: " << data[0][0].size() << "\n";
	std::cout << "level 3:\t size: nE = " << nE << " vector size: " << data[0][0][0].size() << "\n";
} 

void Array4D::fillGrid(){
	data.clear();
	std::vector<double> eVec;
	// loop over inverse order
	// innermost loop over energy
	for(int iE = 0; iE < nE; iE ++){
		eVec.push_back(0.);
	}

	// second inner loop for z-axis
	std::vector< std::vector<double>> zVec;
	for(int iZ = 0; iZ < nZ; iZ++) {
		zVec.push_back(eVec);
	}

	// third inner loop for y-axis
	std::vector< std::vector< std::vector<double>>> yVec;
	for(int iY = 0; iY < nY; iY ++){
		yVec.push_back(zVec);
	}

	// outer loop for x-axis
	std::vector< std::vector< std::vector< std::vector<double>>>> xVec;
	for(int iX = 0; iX < nX; iX++) {
		xVec.push_back(yVec);
	}
	this-> data = xVec;
}

void Array4D::fillGrid(int nX, int nY, int nZ, int nE){
	this -> nX = nX;
	this -> nY = nY;
	this -> nZ = nZ;
	this -> nE = nE;
	fillGrid();
}

void Array4D::setValue(double v, int iX, int iY, int iZ, int iE){
	data[iX][iY][iZ][iE] = v;
}

ObserverHistogram4D::ObserverHistogram4D() {
	hist = new Array4D(1, 1, 1, 1);
	// std::cout << "observerHistogram4D started \n";
}

void ObserverHistogram4D::getIndexFromValue(double value, std::vector<double> grid, int &I) const {
	// std::cout << "start finding position in Grid:\n";
	for(int i = 0; i < grid.size() - 1; i++) {
		if ((value <= grid[i+1]) & (value > grid[i])) {
			I = i;
			// std::cout << "\t" << "go to return statement at i = " << i << " for value: " << value 
			// 	<< " and grid value: " << grid[i] << "," << grid[i+1] << "\n";
			return;
		}
	}
	std::cout << "\t" << "go to return after for loop and use: i = " << grid.size() - 1 <<"\n";
	I =  grid.size() - 1;
}

void ObserverHistogram4D::saveHistogram(std::string path) {
	std::ofstream outputFile;
	outputFile.open(path);

	// write header information
	outputFile << "# Output of a 4D array. Axes contains X, Y, Z, E column\n";
	outputFile << "# The fastes changing index is for the energy column (last column).";
	outputFile << "# The slowest changing column has a row for each bin \n";
	outputFile << "# column length are: \n";
	outputFile << "# nX: " << hist -> nX << "\t xmin: " << xGrid[0] << "\t xmax: " << xGrid[xGrid.size() - 1] << "\n";
	outputFile << "# nY: " << hist -> nY << "\t ymin: " << yGrid[0] << "\t ymax: " << yGrid[yGrid.size() - 1] << "\n";
	outputFile << "# nZ: " << hist -> nZ << "\t zmin: " << zGrid[0] << "\t zmax: " << zGrid[zGrid.size() - 1] << "\n";
	outputFile << "# nE: " << hist -> nE << "\n";
	outputFile << "# Egrid:";
	for (int i = 0; i < eGrid.size(); i++) {
		outputFile << "\t" << eGrid[i];
	}
	outputFile << "\n";

	// write all data
	for (int iX = 0; iX <( hist -> nX); iX++) {
		for (int iY = 0; iY < (hist -> nY); iY++){
			for (int iZ = 0; iZ < (hist -> nZ); iZ++){
				for (int iE = 0; iE < (hist -> nE); iE++){
					outputFile << hist -> data[iX][iY][iZ][iE] << " ";
				}
			}
		}
		outputFile << "\n";
	}

	outputFile.close();
}


void ObserverHistogram4D::fillCandidateInHistogram(const double x, const double y, const double z, const double e, const double w, Array4D &data) const {
	int iX, iY, iZ, iE;
	getIndexFromValue(x, xGrid, iX);
	getIndexFromValue(y, yGrid, iY);
	getIndexFromValue(z, zGrid, iZ);
	getIndexFromValue(e, eGrid, iE);

	hist->data[iX][iY][iZ][iE] += w;
}

void ObserverHistogram4D::setXGrid(double min, double max, int n){
	hist -> fillGrid(n, hist -> nY, hist -> nZ, hist -> nE);
	xGrid.clear();
	for(int i = 0; i <= n; i++){ 
		xGrid.push_back(min + (max - min) / n * i);
	}
}

void ObserverHistogram4D::setYGrid(double min, double max, int n){
	hist -> fillGrid(hist -> nX, n, hist -> nZ, hist -> nE);
	yGrid.clear();
	for(int i = 0; i <= n; i++){
		yGrid.push_back(min + (max - min) / n * i);
	}
}
void ObserverHistogram4D::setZGrid(double min, double max, int n){
	hist -> fillGrid(hist -> nX, hist -> nY, n, hist -> nE);
	zGrid.clear();
	for(int i = 0; i <= n; i++){
		zGrid.push_back(min + (max - min) / n * i);
	}
}
void ObserverHistogram4D::setEGrid(double min, double max, int n, bool scaleLog){
	// std::cout << "create eGrid: \n";
	hist -> fillGrid(hist -> nX, hist -> nY, hist -> nZ, n);
	eGrid.clear();
	if (scaleLog) { // log scaling
		double lgMin = log10(min);
		double lgMax = log10(max);
		// std::cout << "log scaling in range: " << lgMin << " " << lgMax << "\n";
		for(int i = 0; i <= n; i++) {
			double lg = lgMin + (lgMax - lgMin) / n * i;
			// std::cout << lg << "\t";
			eGrid.push_back(pow(10, lg));
		}
	}
	else { // linear scaling
		for(int i = 0; i <= n; i++){
			xGrid.push_back(min + (max - min) / n * i);
		}
	}
}

void ObserverHistogram4D::process(Candidate *cand) const {
	// std::cout << "start process in hist module \n";
	const double E = cand -> current.getEnergy();
	Vector3d pos = cand -> current.getPosition();
	// check for min and max of grid: 
	if( (E < eGrid[0]) || (E > eGrid[eGrid.size() -1])) {return;}
	if( (pos.x < xGrid[0]) || (pos.x > xGrid[xGrid.size() -1])) {return;}
	if( (pos.y < yGrid[0]) || (pos.y > yGrid[yGrid.size() -1])) {return;}
	if( (pos.z < zGrid[0]) || (pos.z > xGrid[zGrid.size() -1])) {return;}
	


	const double w = cand -> getWeight();
	// std::cout << "get properties from candidate: -- done \n";
	int iX, iY, iZ, iE;
	getIndexFromValue(pos.x, xGrid, iX);
	getIndexFromValue(pos.y, yGrid, iY);
	getIndexFromValue(pos.z, zGrid, iZ);
	getIndexFromValue(E, eGrid, iE);

	// std::cout << "get index done: (" << iX << "," << iY << "," << iZ << "," << iE << ")\n";
	// hist -> printHistSize();
	
	// if((iX > (hist->nX)) || (iY > (hist -> nY)) || (iZ > (hist -> nZ)) || (iE >(hist -> nE))) {
	// 	std::stringstream ss;
	// 	ss << "index not okay. Index is (" << iX << "," << iY << "," << iZ << "," << iE<< ")\n";
	// 	ss << "maximal array size is (" << hist-> nX << "," << hist -> nY << "," << (hist -> nZ) << "," << (hist-> nE) << ")\n";
	// 	ss << "candidate: " << cand -> getDescription()<<"\n";
	// 	ss << "\t zGrid: ";
	// 	for(int i = 0; i < zGrid.size(); i++){
	// 		ss << "  " << zGrid[i];
	// 	}
	// 	ss << "\n";
	// 	throw std::runtime_error(ss.str());
	// }

	double w0 = hist -> data[iX][iY][iZ][iE];
	hist -> setValue(w0 + w, iX, iY, iZ, iE);
}

// ---------------------------------------------------------------------------------------------

// inline double fmodulo (double v1, double v2) {
// 	// modulo operation from planck math utils
// 	using namespace std;
// 	if (v1>=0)
//     return (v1<v2) ? v1 : fmod(v1,v2);
// 	double tmp=fmod(v1,v2)+v2;
// 	return (tmp==v2) ? 0. : tmp;
// }

// ObserverHistogram::ObserverHistogram(int nSide, std::string folder, std::string saveName): folder(folder), saveName(saveName) {
// 	setNSide(nSide);
// }

// void ObserverHistogram::process(Candidate *cand) const {
// 	Vector3d pos = cand -> current.getPosition();
// 	if (onEarth) {
// 		pos += Vector3d(8.5 * kpc, 0., 0.);
// 	}
// 	int bin = getBinNumber(pos);
// 	int SN0 = cand -> getSourceSerialNumber();

// #pragma omp critical
// 	{	
// 		std::cout << "start writing \n";
// 		(*file_outs[bin]) << SN0 << "\n";
// 	}

// 	// deactivate after detection
// 	cand -> setActive(false);
// }

// void ObserverHistogram::initFiles() {
// 	for (int i = 0; i < npix; i++){
// 		std::stringstream file;
// 		file << folder << "/" << i << "_" << saveName;

// 		// std::shared_ptr<std::ofstream> outfile(new std::ofstream);
// 		// outfile -> open(file.str());
// 		file_outs.push_back(new std::ostream(file.str(), std::ios::binary));
// 		std::stringsteam out;
// 		out << "# ObserverHistogram Bin: " << i << "\n";
// 		out << "# List of serial number at source \n";
// 		out << "# CRPropa version: " << g_GIT_DESC << "\n#\n";
// 		file_outs[i] -> write(out.str());
// 	}
// }

// void ObserverHistogram::setNSide(int n) {
// 	nSide = n;
// 	npface = nSide * nSide;
// 	ncap = (npface - nSide)<<1;
// 	npix = 12 * npface;
// 	initFiles();
// }

// int ObserverHistogram::getBinNumber(Vector3d &position) const {
// 	/**
// 	 * @brief bin position for a position in the healpix ordering
// 	 * 
// 	 * Healpix bin position for a given position in centered coordinates. 
// 	 * Only Healpix Ordering RING is supported. 
// 	 * This subroutine is extracted from healpix_base.cc (Version 3.82 of healpix)
// 	 */
// 	double theta = position.getTheta();
// 	double z = std::cos(theta);
// 	double za = abs(z);
// 	double phi = position.getPhi();

// 	double tt = fmodulo(phi / M_PI_2, 4.0); // in [0,4)
	
// 	if(za < 2./3.) { // Equatorial region
// 		int nl4 = 4 * nSide;
// 		double temp1 = nSide * (0.5 + tt);
// 		double temp2 = nSide * z * 0.75;
// 		int jp = int(temp1 - temp2); // index of   ascending edge line
// 		int jm = int(temp1 + temp2); // index of descendinge edge line

// 		// ring number counted from z = 2/3
// 		int ir = nSide + 1 + jp - jm; // in {1, 2n+1}
// 		int kshift = 1 - (ir&1); // kshift = 1 if ir even, 0 otherwise

// 		int t1 = jp + jm - nSide + kshift + 1 + nl4 * 2;
// 		int ip = (t1>>1) & (nl4 - 1);	// in {ß0, 4n -1}, case order_>0 from healpix_base.cc
// 		return ncap + (ir - 1) * nl4 + ip;
// 	}
// 	else { // North & South polar caps
// 		bool have_sth = false;
// 		double sth = 0.;
// 		if((theta < 0.01) or (theta > M_PI - 0.01)) {
// 			have_sth = true;
// 			sth = std::sin(theta);
// 		}
// 		double tp = tt - int(tt);
// 		double tmp = ((za < 0.99)||(!have_sth)) ? (nSide * std::sqrt(3 * (1 - za))) : nSide * sth / std::sqrt((1 + za) / 3.);
// 		int jp = int(tp * tmp); // increasing edge line index
// 		int jm = int((1. - tp) * tmp); // decreasing edge line index

// 		int ir = jp + jm + 1; // ring number counted from closest pole
// 		int ip = int(tt * ir); // in {0, 4 * ir - 1}

// 		return (z>0) ? 2 * ir * (ir - 1) + ip : npix - 2 * ir * (ir - 1) + ip;
// 	}
// }

// int ObserverHistogram::getNSide() const {
// 	return nSide;
// }

// void ObserverHistogram::closeAll() {
// 	for (int i = 0; i < npix; i++)
// 	{
// 		file_outs[i] -> close();
// 	}
	
// }
