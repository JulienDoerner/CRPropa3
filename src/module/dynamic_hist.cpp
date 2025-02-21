#include "crpropa/module/dynamic_hist.h"
#include "crpropa/Vector3.h"
#include "crpropa/GridTools.h"

#include <string>
#include <fstream>
#include <thread> 


using namespace crpropa;

const std::string myPropertyName = "counter";

DynamicHistogram::DynamicHistogram() : Module() { }

DynamicHistogram::DynamicHistogram(double Emin, double Emax, int nBins, GridProperties &gridProperties, std::string filename) : Module(), Emin(Emin), Emax(Emax), nBins(nBins), gridProperties(&gridProperties) {
	// init the spatial grids
	initGrid();
	nSpatialBins = gridProperties.Nx * gridProperties.Ny * gridProperties.Nz;
}

void DynamicHistogram::initGrid() {
	// init bins and grid
	energyBins.push_back(Emin);
	for (int i = 0; i < nBins; i++) {
		// center of the energy bin
		double E = Emin * pow(Emax / Emin, (i + 0.5)/ nBins);
		energyCenter.push_back(E);

		// upper edge of the energy bin
		energyBins.push_back(Emin * pow(Emax / Emin, (i + 1.) / nBins));

		// create histogram
		histograms[energyCenter[i]] = new Grid1f(*gridProperties);
	}
}

void DynamicHistogram::process(Candidate *cand) const {

	// get index of the particle, add to position 
	int index = indexFromPosition(cand->current.getPosition());
	int bin = getEnergyBin(cand->current.getEnergy());
	Vector3d pos_index = indiciesFromPosition(cand->current.getPosition());
	int iX = pos_index.x;
	int iY = pos_index.y;
	int iZ = pos_index.z;
	double w0 = candidate -> getWeight();

	// check bounds of indicies
	if (pos_index.x < 0 || pos_index.x >= gridProperties->Nx ||
		pos_index.y < 0 || pos_index.y >= gridProperties->Ny ||
		pos_index.z < 0 || pos_index.z >= gridProperties->Nz ||
		bin < 0 || bin >= energyCenter.size() ||
		index < 0 || index >= nSpatialBins) {

		std::cerr << "Error: particle out of grid" << std::endl;
		std::cerr << "Particle : " << cand -> getDescription() << std::endl;
		std::cerr << "Energy bin: " << bin << " Energy: " << cand -> current.getEnergy() << std::endl;
		std::cerr << "Position: " << cand -> current.getPosition() << std::endl;
		std::cerr << "Index: " << pos_index << std::endl;
		return;
	}

#pragma omp critical (updateHistogram)
{
	double w0 = histograms[energyCenter[bin]] -> get(iX, iY, iZ) += w0;
}

	// add the current bin to the candidate property 
	if (!)
}

int DynamicHistogram::getEnergyBin(double E) const {
	
	// out of range
	if (E >= energyBins.back())
		return energyBins.size() - 1;
	if (E <= energyBins[0])
		return 0;

	// find bin
	for (int i = 0; i < energyBins.size(); i++) {
		if (E < energyBins[i])
			return i - 1;
	}
}

Vector3d DynamicHistogram::indiciesFromPosition(Vector3d pos) const {
	// assumes the position is within the grid 
	Vector3d r = (pos - gridProperties->origin) / gridProperties->spacing;
	return r.floor();
}

int DynamicHistogram::indexFromPosition(Vector3d pos) const {
	Vector3<int> ind = indiciesFromPosition(pos);
	int index =  ind.x * gridProperties -> Ny * gridProperties -> Nz + ind.y * gridProperties -> Nz + ind.z;
	// boundaries: clip between 0 and nx * ny * nz - 1
	index = std::max(0, index);
	index = std::min(index, int(gridProperties -> Nx * gridProperties -> Ny * gridProperties -> Nz - 1));
	return index;
}

void DynamicHistogram::storeGrid(std::string directory, int i) const {
	std::string filename = directory + "/grid_" + std::to_string(i) + ".dat";
	std::cout << "save Grid for energy " << energyCenter[i] << " to " << filename << std::endl;

	// dumping without any conversion but including the grid properties
	dumpGridToTxt(histograms[energyCenter[i]], filename, 1, true);
}

void DynamicHistogram::storeAllGrids(std::string directory) const {
	for (int i = 0; i < energyCenter.size(); i++) {
		storeGrid(directory, i);
	}
}

hid_t DynamicHistogram::getDset(int index) const {
	if (SN_dset.find(index) == SN_dset.end()) {
		std::cerr << "Error: dataset not found" << std::endl;
	}
	return SN_dset[index];
}

void DynamicHistogram::flushBuffer() const {
	const_cast<DynamicHistogram*>(this) -> nParticles = 0; // reset the number of particles
	std::cout << "flush buffer "<< std::endl;

	if(buffer.size() == 0) {
		std::cerr << "Error: no particles in the buffer" << std::endl;
		return;
	}

	// collect all particles for each spatial bin 
	std::vector<int> keys;
	for (auto it = buffer.begin(); it != buffer.end(); ++it) {
		keys.push_back(it -> first);
	}

#pragma omp parallel for
	for (int iK = 0; iK < keys.size(); iK++) {
		int i = keys[iK];
		std::cout << "flush spatial bin " << i << " from thread PID = " << std::this_thread::get_id() << std::endl;
		if (buffer.find(i) == buffer.end())
			continue; // index not in buffer

		SNList &SNs = buffer[i];
		hid_t dset = getDset(i);

		int n = SNs.size();
		
		if (n == 0) {
			std::cerr << "Error: no particles in the buffer" << std::endl;
			continue; // no particles to flush from this buffer 
		}

		hid_t file_space = H5Dget_space(dset);
		hsize_t count = H5Sget_simple_extent_npoints(file_space);

		// resize the dataset
		hsize_t new_count[1] = {count + n};
		H5Dset_extent(dset, new_count);
		
		// get updated space
		H5Sclose(file_space);
		file_space = H5Dget_space(dset);

		hsize_t offset[1] = {count};
		hsize_t cnt[1] = {n};

		// write the data
		H5Sselect_hyperslab(file_space, H5S_SELECT_SET, offset, NULL, cnt, NULL);
		hid_t mspace_id = H5Screate_simple(1, cnt, NULL);
		
		H5Dwrite(dset, H5T_NATIVE_INT, mspace_id, file_space, H5P_DEFAULT, SNs.data());

		H5Sclose(mspace_id);
		H5Sclose(file_space);

		H5Fflush(file, H5F_SCOPE_GLOBAL);
	}
	
	// reset the buffer
	buffer.clear();
}

void DynamicHistogram::close() {
	// close all datasets
	for (int i = 0; i < SN_dset.size(); i++) {
		H5Dclose(SN_dset[i]);
	}
	H5Fclose(file);
}


// -----	DEBUGING functions

void DynamicHistogram::printEnergyBins() const {
	std::cout << "Energy bins: (total  " << energyBins.size() << " bin edges)" << std::endl;
	for (int i = 0; i < energyBins.size(); i++) {
		std::cout << energyBins[i] << " ";
	}
	std::cout << std::endl;
}

void DynamicHistogram::printEnergyCenter() const {
	std::cout << "Energy center: (total  " << energyCenter.size() << " positions)" << std::endl;
	for (int i = 0; i < energyCenter.size(); i++) {
		std::cout << energyCenter[i] << " ";
	}
	std::cout << std::endl;
}

void DynamicHistogram::printGridProperties() const {
	std::cout << gridProperties-> getDescription();
}


void DynamicHistogram::printGrid(int i) const {
	if (i == -1) {
		std::cout << "print all histograms" << std::endl;
		for(int j = 0; j < energyCenter.size(); j++) {
			printGrid(j);
		}
		return;
	}

	std::cout << "print the grid j = " << i << " for energy " << energyCenter[i] << std::endl;
	
	for (int iX = 0; iX < gridProperties->Nx; iX++) {
		for(int iY = 0; iY < gridProperties->Ny; iY++) {
			for(int iZ = 0; iZ < gridProperties->Nz; iZ++) {
				std::cout << "iX/iY/iZ " << iX << "/" << iY <<"/" << iZ << " : " << histograms[energyCenter[i]] -> get(iX, iY, iZ) << std::endl;
			}
		}
	}
}