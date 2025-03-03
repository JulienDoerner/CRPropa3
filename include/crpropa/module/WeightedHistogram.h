#ifndef CRPROPA_WEIGHEDHISTOGRAM_H
#define CRPROPA_WEIGHEDHISTOGRAM_H

#include "crpropa/module/Output.h" 
#include "crpropa/Vector3.h"
#include "crpropa/Grid.h"

#include <vector>
#include <string>
#include <fstream>

using namespace crpropa;

class WeightedHistogram : public Output {
  private: 
    // spatial and energy grid
    GridProperties gridProperties;           //< information about the spatial grid
    std::vector<double> energyBins;         //< energy bins (center) 
    bool logspace = true;                   //< logspace or linspace energy bins
    double dE;                              //< energy bin width (or log width)
    std::vector<ref_ptr<Grid1f> > gridList; //< list of grids for each energy bin

    // SN List
    typedef std::vector<int> tSNList; 
    mutable std::vector<std::vector<std::vector<tSNList> > > SNList; //< list of SN for each energy bin and each cell
    ref_ptr<Output> fullOutput; //< full output of all candidates to map the SN 
    mutable tSNList allSN;     //< list of all SN for the full output

    // buffer for candidates
    int bufferSize = 1e8;   //< number of candidates before writing. Default less than 1 GB in memory
    mutable int bufferCounter = 0;          //< counter for candidates

    // output of the SN List
    std::string directory;                  //< directory for output
    std::string filename = "HistOutput";    //< filename (Prefix) for output. Output file will be <Prefix>_<outputCounter>.txt
    mutable int outputCounter = 0;                  //< number of files written

  public: 
    WeightedHistogram(GridProperties gp, double Emin, double Emax, int nEbins, std::string dir, ref_ptr<Output> fullOutput = NULL);
    
    // destructor
    ~WeightedHistogram() {
        flushBuffer();
        saveGrid();    
    }

    // set the energy bins in a logspace or linsspace grid 
    void initEnergyRange(double Emin, double Emax, int nEbins, bool logspace = true);
    void initBuffer() const;
    void setBufferSize(int n);

    // process the candidate
    void process(Candidate *candidate) const;

    // get spatial index in grid
    Vector3d getGridIndex(const Vector3d &pos) const; 

    // get energy index in grid (gives -1 for out of range)
    int getEnergyIndex(double E) const; 

    // saving results
    void flushBuffer() const;
    void saveGrid() const;

    // debuging information 
    void printEnergyGrid(double scale = 1) const; 
    void printGrid(int i, double scale = 1) const;
    void printAllSN() const;
    void printBuffer() const; 
    void printDirectory() const;

};


#endif // CRPROPA_WEIGHEDHISTOGRAM_H