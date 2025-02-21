#ifndef CRPROPA_DYNAMIC_HIST_H
#define CRPROPA_DYNAMIC_HIST_H

#include <crpropa/Module.h>
#include <crpropa/module/Output.h>
#include <crpropa/module/HDF5Output.h>
#include <crpropa/Units.h>
#include <crpropa/Grid.h>

#include <vector>
#include <list>
#include <map>
#include <omp.h>


using namespace crpropa;

/**
	@class DynamicHistogram
	@brief Module to create a dynamic histogram in energy and space

	This module creates a weighted histogram of the Candidate Population based on the current position and energy of the particle at the detection / process. 
	The final histogram is stored in a HDF5 file. 

	Addionally the current bin in space and energy is stored in the candidate as a property. The bin indices are stored as a vector with different keys. 
 */
class DynamicHistogram : public Module {
  private: 

	// energy histogram 
	double Emin, Emax; //< minimum and maximum energy of the histogram
	int nBins; //< number of bins in energy 
	std::vector<double> energyCenter; //< mid point of the energy bins for the histograms
	std::vector<double> energyBins; //< energy bins for the histograms
	
	// spatial grid for weighted sum 
	mutable std::map<double, ref_ptr<Grid1f> > histograms; //< maps energy and histogram
	GridProperties* gridProperties; //< properties of the spatial grid
	int nSpatialBins; //< number of spatial bins
	
	// enabel debug mode 
	bool debug = false; 

	// storing the current bins in the particle 
	std::string myPropertyName = "DH"; //< (pre)name of the property 

  public:

	/** Default Constructor */
	DynamicHistogram();

	/**  Constructor
	* @param Emin minimum energy of the histogram
	* @param Emax maximum energy of the histogram
	* @param nBins number of bins in energy
	* @param gridProperties properties of the spaitial grid	
	* @param filename name of the file to store the histograms
	*/
	DynamicHistogram(double Emin, double Emax, int nBins, GridProperties &gridProperties, std::string filename);

	/* initilize the grid from the parameters*/
	void initGrid(); 

	/* main process function */
	void process(Candidate *candidate) const; 

	// energy bin for a given energy
	int getEnergyBin(double E) const;

	// spatial bin for a given position
	int indexFromPosition(Vector3d pos) const;
	Vector3d indiciesFromPosition(Vector3d pos) const;

	// store the grid to a file 
	void storeGridsTXT(std::string filename) const;
	void storeGridsHDF(std::string filename) const;

	// debuging / printing function 
	void printEnergyBins() const;
	void printEnergyCenter() const;
	void printGridProperties() const; 
	void printGrid(int i = -1) const;
	void enableDebug() {debug = true;}
};

#endif // CRPROPA_DYNAMIC_HIST_H