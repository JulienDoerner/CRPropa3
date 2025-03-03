#include "crpropa/module/WeightedHistogram.h"
#include "crpropa/module/TextOutput.h"
#include "crpropa/GridTools.h"

#include <sstream>
#include <iomanip>

using namespace crpropa;

WeightedHistogram::WeightedHistogram(
    GridProperties gp, double Emin, double Emax, int nEbins, std::string dir, ref_ptr<Output> fullOutput) :
    gridProperties(gp), directory(dir), Output() {
    
    // set full output
    if(fullOutput.valid())
        this->fullOutput = fullOutput;
    else
        this->fullOutput = new TextOutput(dir + "/FullOutput.txt");
        
    initEnergyRange(Emin, Emax, nEbins, true);
    initBuffer();
};

void WeightedHistogram::initEnergyRange(double Emin, double Emax, int nBins, bool logspace) {
    this -> logspace = logspace;
    // set energy bins 
    if (logspace) {
        double logEmin = log10(Emin);
        double logEmax = log10(Emax);
        dE = (logEmax - logEmin) / (nBins - 1);
        for (int i = 0; i < nBins; i++) {
            energyBins.push_back(logEmin + i * dE);
        }
    } else {
        dE = (Emax - Emin) / (nBins  - 1);
        for (int i = 0; i < nBins; i++) {
            energyBins.push_back(Emin + i  * dE);
        }
    }

    // set the gridList 
    for (int i = 0; i < nBins; i++) {
        gridList.push_back(new Grid1f(gridProperties));
    }
}

void WeightedHistogram::initBuffer() const {
    SNList.clear();
    // init the spatial (buffer) grid
    for(int iX = 0; iX < gridProperties.Nx; iX++) {
        std::vector<std::vector<tSNList> > yzList;
        for(int iY = 0; iY < gridProperties.Ny; iY++) {
            std::vector<tSNList> zList;
            for(int iZ = 0; iZ < gridProperties.Nz; iZ++) {
                zList.push_back(tSNList());
            }
            yzList.push_back(zList);
        }
        SNList.push_back(yzList);
    }
}

void WeightedHistogram::setBufferSize(int n) {
    bufferSize = n;
}

void WeightedHistogram::process(Candidate *cand) const { 
    // check if candidate is in allSN
    int SN = cand->getSerialNumber();
    if (std::find(allSN.begin(), allSN.end(), SN) == allSN.end()) {
#pragma omp critical (WH_fullOutput)
        allSN.push_back(SN);
        fullOutput->process(cand);
    }

    // get energy index 
    int iE = getEnergyIndex(cand->current.getEnergy());
    if (iE == -1) 
        return; // energy out of range

    // get grid index
    Vector3d index = getGridIndex(cand->current.getPosition());
    int iX = index.x;
    int iY = index.y;
    int iZ = index.z;

    // add weight to the grid 
    double w = cand->getWeight();
#pragma omp atomic
    gridList[iE]->get(iX, iY, iZ) += w;

    // add SN to the buffer 
#pragma omp critical(WH_writeBuffer)
    SNList[iX][iY][iZ].push_back(SN);

#pragma omp atomic 
    bufferCounter++;

    // flush buffer if necessary
    if (bufferCounter >= bufferSize) {
        flushBuffer();
        bufferCounter = 0;
    }

}

Vector3d WeightedHistogram::getGridIndex(const Vector3d &pos) const {
    return (pos - gridProperties.origin) / gridProperties.spacing; 
}

int WeightedHistogram::getEnergyIndex(double Energy) const {
    double E = Energy;
    if (logspace) 
        E = std::log10(Energy);

    // lower bondary 
    if (E < energyBins[0] - dE / 2.)
        return -1;

    // in grid    
    for (int i = 0; i < energyBins.size(); i++) {
        if (E < energyBins[i] + dE / 2.)
            return i;
    }

    // upper bondary
    return -1;
}

void WeightedHistogram::flushBuffer() const {
    // filename and open file
    std::stringstream ss; 
    ss << directory << "/" << filename << "_";
    ss << std::setw(4) << std::setfill('0') << outputCounter << ".txt";
    std::cout << ss.str() << std::endl;

    std::ofstream file(ss.str().c_str());
    if (!file.is_open()) {
        std::cerr << "Error: Could not open file " << ss.str() << std::endl;
        return;
    }

    // write header
    file << "# Serial number list for all particles in the grid" << std::endl;
    file << "# The grid uses the following properties: " << std::endl;
    file << "# " << gridProperties.getDescription();
    file << "# Each row reads as (iX, iY, iZ): SN,SN,SN,... " << std::endl;

#pragma omp critical(WH_writeBuffer) 
{
    // write buffer
    for(int iX = 0; iX < gridProperties.Nx; iX++) {
        for(int iY = 0; iY < gridProperties.Ny; iY++) {
            for(int iZ = 0; iZ < gridProperties.Nz; iZ++) {
                tSNList sn = SNList[iX][iY][iZ];
                file << "("<< iX << "," << iY << "," << iZ << "): ";
                if (sn.size() > 0) {
                    for (int i = 0; i < sn.size(); i++) {
                        file << sn[i];
                        if (i < sn.size() - 1)
                            file << ",";
                    }
                }
                file << std::endl;
            }
        }
    }

    // reset buffer
    initBuffer();
    outputCounter++;
    bufferCounter = 0;
}

    // finally close file 
    file.close();
}

void WeightedHistogram::saveGrid() const {
    // store the grid as text files 
    for(int iE = 0; iE < energyBins.size(); iE++) {
        double E = energyBins[iE];
        if (logspace)
            E = pow(10, E);
        std::stringstream ss;
        ss << directory << "/Grid_" << iE << "_" << E / GeV <<"_GeV.txt";
        std::string file = ss.str();

        // save grid with properties
        dumpGridToTxt(gridList[iE], file, 1., true); 
    }
}

// ------------------- Debugging Information -------------------

void WeightedHistogram::printEnergyGrid(double scale) const {
    std::cout << "Energy Grid: ";
    for (int i = 0; i < energyBins.size(); i++) {
        if (logspace)
            std::cout << pow(10, energyBins[i]) / scale << " ";
        else
            std::cout << energyBins[i] / scale << " ";
    }
    std::cout << std::endl;
}

void WeightedHistogram::printGrid(int i, double scale) const {
    std::cout << "Grid " << i << ": ";
    std::vector<float> grid = gridList[i]->getGrid();
    for (int j = 0; j < grid.size(); j++) {
        std::cout << grid[j] / scale << " ";
    }
    std::cout << std::endl;
}

void WeightedHistogram::printAllSN() const {
    std::cout << "All SN: ";
    for (int i = 0; i < allSN.size(); i++) {
        std::cout << allSN[i] << " ";
    }
    std::cout << std::endl;
}

void WeightedHistogram::printBuffer() const {
    std::cout << "Buffer: " << bufferCounter << std::endl;

    for(int iX = 0; iX < SNList.size(); iX++) {
        for(int iY = 0; iY < SNList[iX].size(); iY++) {
            for(int iZ = 0; iZ < SNList[iX][iY].size(); iZ++){
                tSNList sn = SNList[iX][iY][iZ];
                if (sn.size() > 0) {
                    std::cout << "( " << iX << ", " << iY << ", " << iZ << "): ";
                    for (int i = 0; i < sn.size(); i++) {
                        std::cout << sn[i] << " ";
                    }
                    std::cout << std::endl;
                }
            }   
        }
    }
}

void WeightedHistogram::printDirectory() const {
    std::cout << "Directory: " << directory << std::endl;
    std::cout << "Filename: " << filename << std::endl;
}