#include <crpropa/module/HadronicInteractionKamae.h>
#include <crpropa/Common.h>
#include <crpropa/ParticleID.h>
#include <crpropa/Random.h>

using namespace crpropa;

HadronicInteractionKamae::HadronicInteractionKamae(ref_ptr<Density> dens, double limit) : massDensity(dens), limit(limit) { 
    writeOut = new std::ofstream("Kamae_Module_writeout.txt");
    *writeOut << "process" << "\t" << "Ep [TeV]" << "\t" << "Esec [TeV]" << "\t" << "ND" << "\t" <<  "Diff" << "\t" << "1232" << "\t" << "1600" << "\n";

    write2 = new std::ofstream("Kamae_eq9.txt");
    *write2 << "x" << "\t" << "arg1" << "\t" << "arg2" << "\t"<<"x_br" << "\n";

}   


// ---------------------- Fit functions --------------------------------------------------------------

double HadronicInteractionKamae::eq6(double x, std::vector<double> a) const {
    double f1 = a[0] * std::exp(- a[1] * pow_integer<2>(x - a[3] + a[2] * pow_integer<2>(x - a[3])));
    double pow2 = - a[5] * pow_integer<2>(x - a[8] + a[6] * pow_integer<2>(x - a[8]) + a[7] * pow_integer<3>(x - a[8]));
    double f2 = a[4] * std::exp(pow2);
    // writeOut << f1 << " \t" << f2 << "(" << pow2 << ")\n";
    return std::max(f1 + f2, 0.);
}

double HadronicInteractionKamae::eq9(double x, std::vector<double> b) const {
    double arg1 = - b[1] * pow_integer<2>((x - b[2]) / (1 + b[3] * (x - b[2])));
    double arg2 = - b[5] * pow_integer<2>((x - b[6]) / (1 + b[7] * (x - b[6])));
    *write2<< x <<"\t" << b[0] * std::exp(arg1) <<"\t" << b[4] * std::exp(arg2) << "\t" << b[6] - 1. / b[7] << "\n";
    return std::max(b[0] * std::exp(arg1) + b[4] * std::exp(arg2), 0.);
}

double HadronicInteractionKamae::eq12(double x, std::vector<double> c) const {
    double inner = (x - c[2]) / (1 + c[3] * (x - c[2]) + c[4] * pow_integer<2>(x - c[2]));
    return std::max(c[0] * std::exp( -c[1] * pow_integer<2>(inner)), 0.);
}

// ---------------------- General --------------------------------------------------------------

double HadronicInteractionKamae::crossection(double Tp) const {
    // see erratum: ApJ 662:779(2007)

    double m0c2 = mass_proton * c_squared;
    // double p = std::sqrt(pow_integer<2>(Tp + m0c2) - m0c2 * m0c2) / c_light;
    // p = p / GeV * c_light;
    double Ep = (Tp + m0c2) / GeV; // total energy 
    double p = std::sqrt(pow_integer<2>(Ep * GeV) - m0c2 * m0c2) / c_light; 
    p = p / GeV * c_light; // unit conversion GeV/c
    double x = std::log10(p);
    double x2 = x * x;
    double x3 = x2 * x;
    double sigma = 0.;
    
    // parameter from table 1
    double a[8] = {0.1176, 0.3829, 23.10, 6.454, -5.764, -23.63, 94.75, 0.02667};
    double b[2] = {11.34, 23.72};
    double c[3] = {28.5, -6.133, 1.464};
    double d[7] = {0.3522, 0.1530, 1.498, 2.0, 30., 3.155, 1.042};
    double e[2] = {5.922, 1.632};
    double f[5] = {0.0834, 9.5, -5.5, 1.68, 3134};
    double g[5] = {0.0004257, 4.5, -7.0, 2.1, 503.5};

    // non diffractive
    double sigma_nd = 0.;
    if(p < 1)
        sigma_nd = 0.;
    else if(p < 1.3){
        double pow1 = - a[6] * pow_integer<2>(x + a[7]);
        double s1 = a[2];
        double s2 = a[3] * x2;
        double s3 = a[4] * x3;
        double s4 = a[5] * std::exp(pow1);
        double s5 = 0.57 * pow(x / a[0], 1.2); 
        sigma_nd = s5 * (s1 + s2 + s3 + s4);
    }
    else if(p < 2.4) 
        sigma_nd = (b[0] * fabs(a[1] - x) + b[1] * fabs(a[0] - x)) / (a[1] - a[0]);
    else if(p < 10)
        sigma_nd = a[2] + a[3] * x2 + a[4] * x3 + a[5] * std::exp(-a[6] * pow_integer<2>(x + a[7]));
    else 
        sigma_nd = c[0] + c[1] * x + c[2] * x2;
    
    // diffractive
    double sigma_diff = 0.;
    if(p < 2.25)
        sigma_diff = 0;
    else if(p < 3.2)
        sigma_diff = std::sqrt((x - d[0])/d[1]) * (d[2] + d[3] * std::log10(d[4] * (x - 0.25)) + d[5] * x2 - d[6] * x3);
    else if(p < 100) 
        sigma_diff = d[2] + d[3] * std::log10(d[4] * (x - 0.25)) + d[5] * x2 - d[6] * x3;
    else
        sigma_diff = e[0] + e[1] * x;
    
    // res(1232)
    double sigma_1232 = 0.;
    if(Ep < 1.4) 
        sigma_1232 = 0;
    else if(Ep < 1.6) 
        sigma_1232 = f[0] * pow_integer<10>(Ep);
    else if(Ep < 1.8)
        sigma_1232 = f[1] * std::exp(f[2] * pow_integer<2>(Ep - f[3]));
    else if(Ep < 10)
        sigma_1232 = f[4] / pow_integer<10>(Ep);
    else
        sigma_1232 = 0;
    
    // res(1600)
    double sigma_1600 = 0.;
    if(Ep < 1.6) 
        sigma_1600 = 0.;
    else if(Ep < 1.9)
        sigma_1600 = g[0] * pow_integer<14>(Ep);
    else if(Ep < 2.3)
        sigma_1600 = g[1] * std::exp( g[2] * pow_integer<2>(Ep - g[3]));
    else if(Ep < 20)
        sigma_1600 = g[4] / pow_integer<6>(Ep);
    else
        sigma_1600 = 0;
    
    sigma = sigma_nd + sigma_diff + sigma_1232 + sigma_1600;
    return sigma;
}

void HadronicInteractionKamae::process(Candidate* cand) const {
    // extract candidate properties
    int id = cand -> current.getId();
    if(!isNucleus(id))
        return; // only for nucleons
    
    double E = cand -> current.getEnergy();
    double restEnergy = cand -> current.getMass() * c_squared;
    double Tp = E - restEnergy;
    int N = massNumber(id); // scaling crossection for individual nucleous
    double step = cand -> getNextStep();

    // decide if interaction happens
    double cross = this->crossection(Tp);
    Random &random = Random::instance();
    double nH = massDensity -> getNucleonDensity(cand-> current.getPosition());
    double meanFreePath = cross * nH / N;
    
    if(random.rand() > step / meanFreePath){
        cand -> limitNextStep(limit * meanFreePath);
        return;
    }

    cand -> setActive(false);
    performInteraction(cand);
    return;
}

void HadronicInteractionKamae::performInteraction(Candidate* cand) const {
    // return;

    double E = cand -> current.getEnergy();
    double Tp = E - cand -> current.getMass() * c_squared;
    Random &rand = Random::instance();

    double epsMin = 1e-6 * Tp;
    double epsMax = Tp;
    double nPhoton = gaussInt([this, Tp](double x) {return this->SpectrumPhoton(x, Tp); }, epsMin, epsMax);
    double nEminus = gaussInt([this, Tp](double x) {return this -> SpectrumElectron(x, Tp);}, epsMin, epsMax);
    double nEplus = gaussInt([this, Tp](double x) {return this -> SpectrumPositron(x, Tp);}, epsMin, epsMax);
    double nNuE = gaussInt([this, Tp](double x) {return this -> SpectrumElectronNeutrino(x, Tp);}, epsMin, epsMax);
    double nAntiNuE = gaussInt([this, Tp](double x) {return this -> SpectrumAntiElectronNeutrino(x, Tp);}, epsMin, epsMax);
    double nNuMu = gaussInt([this, Tp](double x) {return this -> SpectrumMuonNeutrino(x, Tp);}, epsMin, epsMax);
    double nAntiNuMu = gaussInt([this, Tp](double x) {return this -> SpectrumAntiMuonNeutrino(x, Tp);}, epsMin, epsMax);

    double nTotal = nPhoton + nEminus + nEplus + nNuE + nAntiNuE + nNuMu + nAntiNuMu;

    // sample secondaries while energy left
    double Eleft = Tp;

    // find maximum propability
    int nStep = 1000;
    double pMaxPhoton, pMaxEminus, pMaxEplus, pMaxNuE, pMaxAnitNuE, pMaxNuMu, pMaxAntiNuMu, x;
    for(int i = 0; i < nStep; i++){
        x = epsMin + i * (epsMax - epsMin) / nStep;
        double pPhoton = SpectrumPhoton(x, Tp);
        if(pPhoton > pMaxPhoton)
            pMaxPhoton = pPhoton;
        
        double pEminus = SpectrumElectron(x, Tp);
        if(pEminus > pMaxEminus)
            pMaxEminus = pEminus;
        
        double pEplus = SpectrumPositron(x, Tp);
        if(pEplus > pMaxEminus)
            pMaxEplus = pEplus;
        
        double pNuE = SpectrumElectronNeutrino(x, Tp);
        if(pNuE > pMaxNuE) 
            pMaxNuE = pNuE;
        
        double pAntiNuE = SpectrumAntiElectronNeutrino(x, Tp);
        if(pAntiNuE > pMaxAnitNuE)
            pMaxAnitNuE = pAntiNuE;
        
        double pNuMu = SpectrumMuonNeutrino(x, Tp);
        if(pNuMu > pMaxNuMu)
            pMaxNuMu = pNuMu;
        
        double pAntiNuMu = SpectrumAntiMuonNeutrino(x, Tp);
        if(pAntiNuMu > pMaxAntiNuMu)
            pMaxAntiNuMu = pAntiNuMu;
    }

    while(true) {
        // choose random which particle to sample
        double nPart = rand.rand() * nTotal;
        double Esample = 0.;
        int Id = 0;
        double pMax = 0.;
        std::function<const double(double)> sampleFrom;
        if(nPart < nPhoton){
            // sample Photon   
            sampleFrom = [this, Tp](double x) {return this -> SpectrumPhoton(x, Tp);};
            Id = 22;
            pMax = pMaxPhoton;
        }
        nPart -= nPhoton;
        if(nPart < nEminus) {
            // sample electron
            sampleFrom = [this, Tp](double x) {return this -> SpectrumElectron(x, Tp);};
            Id = 11;
            pMax = pMaxEminus;
        }
        nPart -= nEminus;
        if(nPart < nEplus) {
            // sample positron;
            sampleFrom = [this, Tp](double x) {return this -> SpectrumPositron(x, Tp);};
            Id = -11;
            pMax = pMaxEplus;
        }
        nPart -= nEplus;
        if(nPart < nNuE) {
            // sample electron neutrino
            sampleFrom = [this, Tp](double x) {return this -> SpectrumElectronNeutrino(x, Tp);};
            Id = 12;
            pMax = pMaxNuE;
        }
        nPart -= nNuE;
        if(nPart < nAntiNuE) {
            // sample anti electron neutrino
            sampleFrom = [this, Tp](double x) {return this -> SpectrumAntiElectronNeutrino(x, Tp);};
            Id = -12;
            pMax = pMaxAnitNuE;
        }
        nPart -= nAntiNuE;
        if(nPart < nNuMu) {
            // sample muon neutrino
            sampleFrom = [this, Tp](double x) {return this -> SpectrumMuonNeutrino(x, Tp);};
            Id = 14;
            pMax = pMaxNuMu;
        }
        nPart -= nNuMu;
        if(nPart < nAntiNuMu) {
            // sample anti muon neutrino
            sampleFrom = [this, Tp](double x) {return this -> SpectrumAntiMuonNeutrino(x, Tp);};
            Id = -14;
            pMax = pMaxAntiNuMu;
        }

        // sampling 
        double epsProbe, pFunction, pTest;
        int counter = 0;
        while (true) {
            epsProbe = rand.rand() * epsMax;
            pFunction = sampleFrom(epsProbe);
            pTest = rand.rand() * pMax;
            counter++;
            if(pTest < pFunction)
                break;
            if((counter / 10 )== 0) {
                // *writeOut << counter;
            }
        }
        
        // check if remaining energy is positiv and create secondary
        Eleft -= epsProbe;
        if(Eleft > 0) {
            cand -> addSecondary(Id, epsProbe);
        }
        else {
            // decide randomly to keep the secondary and break sampling
            pTest = rand.rand();
            if(pTest > abs(Eleft) / epsProbe) {
                cand -> addSecondary(Id, epsProbe);
            }
            return;
        }
    }
}

// ---------------------- PHOTON --------------------------------------------------------------

double HadronicInteractionKamae::SpectrumPhoton(double Esec, double Tp) const {
    if (Esec >= Tp) 
        return 0.; // nothing left
    
    double s1 = PhotonND(Esec, Tp);
    double s2 = PhotonDiff(Esec, Tp);
    double s3 = Photon1232(Esec, Tp);
    double s4 = Photon1600(Esec, Tp);

    *writeOut << "photon" << "\t" << Tp / TeV << "\t" << Esec / GeV << "\t" << s1 << "\t" << s2 << "\t" << s3 << "\t" << s4 << "\n";
    return s1 + s2 + s3 + s4;
}

double HadronicInteractionKamae::PhotonND(double Esec, double Tp) const {
    if (Esec >= Tp) 
        return 0.; // nothing left
    
    double x = std::log10(Esec / GeV); 
    double y = std::log10(Tp / TeV);

    std::vector<double> a;
    for(int i = 0; i < 9; i++) {a.push_back(0.);}

    double z = y + 3.25;
    double y2 = y * y;
    double y3 = y2 * y;
    double y4 = y3 * y;
    
    double rescale = 1.01;
    if(Tp > 1.95 * GeV) 
        rescale = 3.05 * exp(-107 * pow_integer<2>(z / (1. + 8.08 * z)));

    // tabel 3
    z = y + 3.3;
    a[0] = -0.51187 * z + 7.6179 * pow_integer<2>(z) - 2.1332 * pow_integer<3>(z) + 0.22184 * pow_integer<4>(z);
    a[1] = -1.2592e-5 + 1.4439e-5 * std::exp(- 0.29360 * (y + 3.4)); 
    a[1] += 5.9363e-5 / (y + 4.1485) + 2.264e-6 * y - 3.3723e-7 * y2;
    a[2] = - 174.83 + 152.78 * std::log10(1.5682 * (y + 3.4)) - 808.74 / (y + 4.6157);
    a[3] = 0.81177 + 0.56385 * y + 0.0040031 * y2 - 0.0057658 * y3 + 0.00012057 * y4;
    a[4] = 0.68631 * (y + 3.32) + 10.145 * pow_integer<2>(y + 3.32) - 4.6176 * pow_integer<3>(y + 3.32) + 0.86824 * pow_integer<4>(y + 3.32) - 0.053741 * pow_integer<5>(y + 3.32);
    a[5] = 9.0466e-7 +1.4539e-6 * std::log10(0.0152074 * (y + 3.4)); 
    a[5] += 1.3253e-4 / pow_integer<2>(y + 4.7171) - 4.1228e-7 * y + 2.2036e-7 * y2;
    a[6] = -339.45 + 618.73 * std::log10(0.31595 * (y + 3.9)) + 250.20 / pow_integer<2>(y + 4.4395);
    a[7] = -35.105 + 36.167 * y - 9.3575 * y2 + 0.33717 * y3;
    a[8] = 0.17554 + 0.37300 * y - 0.014938 * y2 + 0.0032314 * y3 + 0.0025579 * y4;

    double sigma = eq6(x, a);

    // kinematic limits from table 2
    double W_lo = 20;
    double W_hi = 45;
    double Lmin = -2.6;
    double Lmax = 0.96 * (y + 3);
    double fND_kl = 1 / (exp(W_lo * (Lmin - x)) + 1) / (exp(W_hi * (x - Lmax)) + 1);

    return sigma * rescale * fND_kl; 
}

double HadronicInteractionKamae::PhotonDiff(double Esec, double Tp) const {
    if (Esec >= Tp) 
        return 0.; // nothing left
    
    double x = std::log10(Esec / GeV); 
    double y = std::log10(Tp / TeV);

    std::vector<double> b; 
    for(int i = 0; i < 8; i++) {b.push_back(0.);}

    double y2 = y * y;
    double y3 = y2 * y;
    double y4 = y3 * y;
    
    // table 3
    // diffraction (eq. 9)
    if(Tp > 5.52 * GeV) {
        b[0] = 60.142 * std::tanh(-0.37555 * (y + 2.2)) - 5.9564 * pow_integer<2>(y + 0.59913) + 6.0162e-3 * pow_integer<4>(y + 9.4773);
        b[1] = 35.322 + 3.8026 * std::tanh(-2.5979 * (y + 1.9)) - 2.1870e-4 * pow_integer<2>(y + 369.13);
        b[2] = -15.732 -0.082064 * std::tanh(-1.9621 * (y + 2.1)) + 2.3355e-4 * pow_integer<2>(y + 252.43);
        b[3] = -0.086827 + 0.37646 * std::exp(- 0.53053 * pow_integer<2>((y + 1.0444) / (1.0 + 0.27437 * (y + 1.0444))));
    }
    b[4] = 2.5982 + 0.39131 * pow_integer<2>(y + 2.95) - 0.0049693 * pow_integer<4>(y + 2.95); 
    b[4] += 0.94131 * std::exp(-24.347 * pow_integer<2>(y + 2.45 - 0.19717 * pow_integer<2>(y + 2.45)));
    b[5] = 0.11198 - 0.64582 * y + 0.16114 * y2 + 2.2853 * std::exp(-0.0032432* pow_integer<2>((y - 0.83562)/(1 + 0.33933 * (y - 0.83562))));
    b[6] = 1.7843 + 0.91914 * y + 0.050118 * y2 + 0.038096 * y3 - 0.027334 * y4 - 0.0035556 * y4 * y + 0.0025742 * y4 * y2;
    b[7] = -0.19870 - 0.071003 * y + 0.019328 * y2 - 0.28321 * std::exp(- 6.0516 * pow_integer<2>(y + 1.8441));

    double f_kl = 1 / (std::exp(75 * (x - (y + 3))) + 1); 
    double f_diff = eq9(x, b);
    return f_diff * f_kl;
}

double HadronicInteractionKamae::Photon1232(double Esec, double Tp) const {
    if (Esec >= Tp) 
        return 0.; // nothing left
    
    double x = std::log10(Esec / GeV); 
    double y = std::log10(Tp / TeV);
    double y2 = y * y;

    std::vector<double> c;
    for(int i = 0; i < 5; i++) {c.push_back(0.);}

    // table 3
    // res(1232), eq. 12
    c[0] = 2.4316 * std::exp( - 69.484 * pow_integer<2>((y + 3.1301) / (1. + 1.24921 * (y + 3.1301)))) - 6.3003 - 9.5349 / y + 0.38121 * y2;
    c[1] = 56.872 + 40.627 * y + 7.7528 * y2;
    c[2] = -5.4918 - 6.7872 * std::tanh(4.7128 * (y + 2.1)) + 0.68048 * y;
    c[3] = - 0.36414 + 0.039777 * y;
    c[4] = -0.72807 - 0.48828 * y - 0.092876 * y2;

    double f_res = eq12(x, c);
    double f_kl = 1 / (std::exp(75 * (x - (y + 3))) + 1); 
    
    return f_res * f_kl;
}

double HadronicInteractionKamae::Photon1600(double Esec, double Tp) const {
    if (Esec >= Tp) 
        return 0.; // nothing left
    
    double x = std::log10(Esec / GeV); 
    double y = std::log10(Tp / TeV);
    double y2 = y * y;
    double y3 = y2 * y;

    std::vector<double> d;
    for(int i = 0; i < 5; i++) {d.push_back(0.);}

    // table 3
    // res(1600), eq. 12
    d[0] = 3.2433 * std::exp(-57.133 * pow_integer<2>((y + 2.9507) / (1. + 1.2912 * (y + 2.9507)))) - 1.0640 - 0.43925 * y;
    d[1] = 16.901 + 5.9539 * y - 2.1257 * y2 - 0.92057 * y3;
    d[2] = -6.6638 - 7.5010 * std::tanh(30.322 * (y + 2.1)) + 0.54662 * y;
    d[3] = -1.50648 - 0.87211 * y - 0.17097 * y2;
    d[4] = 0.42795 + 0.55136 * y + 0.20707 * y2 + 0.027552 * y3;

    double f_res = eq12(x, d);
    double f_kl = 1 / (std::exp(75 * (x - (y + 3))) + 1); 
    return f_res * f_kl;
}


// ---------------------- Electron --------------------------------------------------------------
double HadronicInteractionKamae::SpectrumElectron(double Esec, double Tp) const {
    if (Esec >= Tp) 
        return 0;

    double s1 = ElectronND(Esec, Tp);
    double s2 = ElectronDiff(Esec, Tp);
    double s3 = Electron1232(Esec, Tp);
    double s4 = Electron1600(Esec, Tp);

    *writeOut << "e-" << "\t" << Tp / TeV << "\t" << Esec / GeV << "\t" << s1 << "\t" << s2 << "\t" << s3 << "\t" << s4 << "\n";
    return s1 + s2 + s3 + s4;
}

double HadronicInteractionKamae::ElectronND(double Esec, double Tp) const {
    std::vector<double> a;
    for(short i = 0; i < 9; i++)
        a.push_back(0.);
    
    double y = std::log10(Tp / TeV);
    double y2 = y * y;
    double y3 = y2 * y;
    double y4 = y3 * y;
    double z; // helper variable
    double x = std::log10(Esec / GeV);

    double rescale = 1.01; // for Tp > 1.95 GeV
    if(Tp < 1.95 * GeV) {
        z = y + 3.25;
        double pow = - 107 * pow_integer<2>(z / (1 + 8.8 * z));
        rescale = 3.05 * std::exp(pow);
    }

    // table 4
    z = y + 3.3;
    a[0] = -0.018639 * z + 2.4315 * z * z - 0.57719 * pow_integer<3>(z) + 0.063435 * pow_integer<4>(z);
    a[1] = 7.1827e-6 - 3.5067e-6 * y + 1.3264e-6 * y2 - 3.3481e-7 * y3 + 2.3551e-8 * y4 + 3.4297e-8 * y4 *y;
    z = y + 3.4;
    a[2] = 563.91 - 362.18 * log10(2.7187 * z) - 2.8924e4 / pow_integer<2>(y + 7.9031);
    a[3] = 0.52684 + 0.57717 * y + 0.0045336 * y2 - 0.0089066 * y3;
    z = y + 3.32;
    a[4] = 0.36108 * z + 1.6963 * z * z - 0.074456 * pow_integer<3>(z) - 0.071455 * pow_integer<4>(z) + 0.010473 * pow_integer<5>(z);
    a[5] = 9.7387e-5 + 7.8573e-5 * log10(0.0036055 * (y + 4.3)) + 0.0002466 / (y + 4.9390) - 3.8097e-7 * y2;
    a[6] = -273. - 106.22 * log10(0.341 * (y + 3.4)) + 89.037 * y - 12.546 * y2;
    a[7] = 432.53 - 883.99 * log10(0.19737 * (y + 3.9)) - 4.1938e4 / pow_integer<2>(y + 8.5518);
    a[8] = -0.12756 + 0.43478 * y - 0.0027797 * y2 - 0.0083074 * y3;
    
    // writeOut << x << "\t";
    // for(short i = 0; i < 9; i++) 
        // writeOut << "a" << i<< " "<< a[i] << "\t";
    

    double sigma = eq6(x, a);

    // kinematic limits
    double W_lo = 20;
    double W_hi = 45;
    double Lmin = -2.6;
    double Lmax = 0.96 * (y + 3);
    double fND_kl = 1 / (exp(W_lo * (Lmin - x)) + 1) / (exp(W_hi * (x - Lmax)) + 1);

    return sigma * rescale * fND_kl;  
}

double HadronicInteractionKamae::ElectronDiff(double Esec, double Tp) const {
    std::vector<double> b;
    for(short i = 0; i < 8; i++)
        b.push_back(0.);
    
    double y = std::log10(Tp / TeV);
    double y2 = y * y;
    double y3 = y2 * y;
    double y4 = y3 * y;
    double z; // helper variable
    double x = std::log10(Esec / GeV);

    // table 4
    if(Tp > 5.52 * GeV) {
        b[0] = 0.20463 * std::tanh(-6.237 * (y + 2.2)) - 0.16362 * pow_integer<2>(y + 1.6878) + 3.5183e-4 * pow_integer<4>(y + 9.64);
        z = y + 2.0154;
        b[1] = 1.6537 + 3.853 * std::exp(-3.2027 * pow_integer<2>(z / (1 + 0.62779 * z)));
        b[2] = -10.722 - 0.082672 * std::tanh(-1.8879 * (y + 2.1)) + 1.4895e-4 * pow_integer<2>(y + 256.63);
        z = y + 1.9877;
        b[3] = -0.023752 - 0.51734 * std::exp(-3.3087 * pow_integer<2>(z / ( 1 + 0.403 * z)));
    }
    z = y + 2.9;
    b[4] = 0.94921 + 0.1228 * z * z - 7.1585e-4 * pow_integer<4>(z) + 0.5213 * std::log10(z);
    b[5] = -4.2295 - 1.0025 * std::tanh(9.0733 * (y + 1.9)) - 0.11452 * (y - 62.382);
    b[6] = 1.4862 + 0.99544 * y - 0.042763 * y2 - 0.0040065 * y3 + 0.0057987 * y4;
    b[7] = 6.2629 + 6.9517 * std::tanh(-0.3648 * (y + 2.1)) - 0.026033 * (y - 2.8542);

    // kinematic limits
    double Wdiff = 75.;
    double Lmax = y + 3;
    double fKL = 1 / (exp(Wdiff / (x - Lmax)) + 1);

    return eq9(x, b) * fKL;
}

double HadronicInteractionKamae::Electron1232(double Esec, double Tp) const {
    // see table 4
    return 0.;
}

double HadronicInteractionKamae::Electron1600(double Esec, double Tp) const {
    std::vector<double> d;
    for(short i = 0; i < 5; i++)
        d.push_back(0.);
    
    double y = std::log10(Tp / TeV);
    double y2 = y * y;
    double y3 = y2 * y;
    double z = y + 2.9537; // helper variable
    double x = std::log10(Esec / GeV);

    // table 4
    d[0] = 0.37790 * exp(-56.826 * pow_integer<2>(z / (1 + 1.5221 * z))) - 0.059458 + 0.0096583 * y2;
    d[1] = -5.5135 - 3.3988 * y;
    d[2] = -7.1209 - 7.185 * tanh(30.801 * (y + 2.1)) + 0.35108 * y;
    d[3] = -6.7841 - 4.8385 * y - 0.91523 * y2;
    d[4] = -134.03 - 139.63 * y - 48.316 *  y2 - 5.5526 * y3;

    // kinematic limits
    double Wdiff = 75.;
    double Lmax = y + 3;
    double fKL = 1 / (exp(Wdiff / (x - Lmax)) + 1);

    double sigma = eq12(x, d);
    // writeOut << x << " " << sigma << " " << fKL << "\n";
    return sigma * fKL;
}

// ---------------------- Positron --------------------------------------------------------------

double HadronicInteractionKamae::SpectrumPositron(double Esec, double Tp) const {
    if (Esec >= Tp) 
        return 0;

    double s1, s2, s3, s4;
    s1 = PositronND(Esec, Tp);
    s2 = PositronDiff(Esec, Tp);
    s3 = Positron1232(Esec, Tp);
    s4 = Positron1600(Esec, Tp);


    *writeOut << "e+" << "\t" << Tp / TeV << "\t" << Esec / GeV << "\t" << s1 << "\t" << s2 << "\t" << s3 << "\t" << s4 << "\n";
    return s1 + s2 + s3 + s4;
}

double HadronicInteractionKamae::PositronND(double Esec, double Tp) const {
    std::vector<double> a;
    for(short i = 0; i < 9; i++)
        a.push_back(0.);
    
    double y = std::log10(Tp / TeV);
    double y2 = y * y;
    double y3 = y2 * y;
    double y4 = y3 * y;
    double z = y + 3.25; // helper variable
    double x = std::log10(Esec / GeV);

    double rescale = 1.;
    if (Tp < 5.52 * GeV) {
        rescale = 2.2 * exp( - 98.9 * pow_integer<2>(z / (1 + 10.4 * z)));
    }

    // table 5
    z = y + 3.3;
    a[0] = z * (-0.79606 + z * (7.7496 + z * (-3.9326 + z * (0.80202 - 0.054994 * z))));
    a[1] = 6.7943e-6 - 3.5345e-6 * y + 6.0927e-7 * y2 + 2.0219e-7 * y3 + 5.1005e-8 * y4 - 4.2622e-8 * y4 * y;
    a[2] = 44.827 - 81.378 * log10(0.027733 * (y + 3.5)) - 1.3886e4 / (y + 8.4417);
    a[3] = 0.52010 + 0.59336 * y + 0.012032 * y2 - 0.0064242 * y3;
    z = y + 3.32;
    a[4] = z *( 2.1361 + z * (1.8514 + z *(-0.47872 + z * (0.0032043 + z * 0.0082955))));
    a[5] = 1.0845e-6 + 1.4336e-6 * log10(0.0077255 * (y + 4.3)) + 1.3018e-4 / pow_integer<2>(y + 4.8188) + 9.3601e-8 * y;
    a[6] = -267.74 + 14.175 * log10(0.3591 * (y + 3.4)) + 64.669 / pow_integer<2>(y - 7.7036); 
    a[7] = 138.26 - 539.84 * log10(0.12467 * (y + 3.9)) - 1.9869e4 / pow_integer<2>(y + 7.6884) + 1.0675 * y2;
    a[8] = -0.14707 + 0.40135 * y + 0.0039899 * y2 - 0.0016602 * y3;

    double sigma = eq6(x, a);

    // kinematic limits
    double W_lo = 15;
    double W_hi = 47;
    double Lmin = -2.6;
    double Lmax = 0.94 * (y + 3);
    double fND_kl = 1 / (exp(W_lo * (Lmin - x)) + 1) / (exp(W_hi * (x - Lmax)) + 1);

    return rescale * sigma * fND_kl;
}

double HadronicInteractionKamae::PositronDiff(double Esec, double Tp) const {
    std::vector<double> b;
    for(short i = 0; i < 9; i++)
        b.push_back(0.);
    
    double y = std::log10(Tp / TeV);
    double y2 = y * y;
    double y3 = y2 * y;
    double y4 = y3 * y;
    double z = y + 1.8781; // helper variable
    double x = std::log10(Esec / GeV);    

    // table 5
    if(Tp > 11.05 * GeV) {
        b[0] = 29.192 * tanh(-0.37879 * (y + 2.2)) - 3.2196 * pow_integer<2>(y + 0.675) + 0.0036687 * pow_integer<4>(y + 9.0824);
        b[1] = -142.97 + 147.86 * exp(-0.37194 * pow_integer<2>(z / (1. + 3.8389 * z)));
        b[2] = -14.487 - 4.2223 * tanh(-13.546 * (y + 2.2)) + 1.6988e-4 * pow_integer<2>(y + 234.65);
        z = y + 1.8194;
        b[3] = -0.0036974 - 0.41976 * exp(-6.1527 * pow_integer<2>(z / (1. + 0.99946 * z)));
    }
    z = y + 2.9;
    b[4] = 1.8108 + z * z * (0.18545  - 0.0020049 * z * z) + 0.85084 * exp(- 14.987 * pow_integer<2>(y + 2.29 - 0.18967 * pow_integer<2>(y + 2.29)));
    b[5] = 2.0404 - 0.51548 * tanh(2.2758 * (y + 1.9)) - 0.035009 / (y - 6.6555);
    b[6] = 1.5258 + 1.0132 * y - 0.064388 * y2 - 0.0040209 * y3 - 0.0082772 * y4;
    b[7] = 3.0551 + 3.5240 * tanh(-0.36739 * (y + 2.1)) - 0.13382 * (y - 2.7718);

    double sigma = eq9(x, b);

    // kinematic limits (for all species, diff and res)
    double Wdiff = 75.;
    double Lmax = y + 3;
    double fKL = 1 / (exp(Wdiff / (x - Lmax)) + 1);    

    return sigma * fKL;
}

double HadronicInteractionKamae::Positron1232(double Esec, double Tp) const {
    std::vector<double> c;
    for(short i = 0; i < 5; i++)
        c.push_back(0.);
    
    double y = std::log10(Tp / TeV);
    double y2 = y * y;
    double z = y + 3.1272; // helper variable
    double x = std::log10(Esec / GeV);    
    
    // table 5
    c[0] = 2.9841 * exp(-67.857 * pow_integer<2>(z / (1. + 0.22831 * z))) - 6.5855 - 9.6984/y + 0.41256 * y2;
    c[1] = 6.8276 + 5.2236 * y + 1.4630 * y2;
    c[2] = -6.0291 - 6.4581 * tanh(5.083 * (y + 2.1)) + 0.46352 * y;
    c[3] = 0.593 + 0.36093 * y;
    c[4] = 0.77368 + 0.44776 * y + 0.056409 * y2;

    double sigma = eq12(x, c);

    // kinematic limits (for all species, diff and res)
    double Wdiff = 75.;
    double Lmax = y + 3;
    double fKL = 1 / (exp(Wdiff / (x - Lmax)) + 1);       

    return sigma * fKL;
}

double HadronicInteractionKamae::Positron1600(double Esec, double Tp) const {
    std::vector<double> d;
    for(short i = 0; i < 5; i++)
        d.push_back(0.);
    
    double y = std::log10(Tp / TeV);
    double y2 = y * y;
    double y3 = y2 * y;
    double z = y + 2.9485; // helper variable
    double x = std::log10(Esec / GeV);    
    
    // table 5
    d[0] = 1.9186 * exp(- 56.544 * pow_integer<2>(z / (1. + 1.2892 * z))) - 0.2372 + 0.041315 * y2;
    d[1] = -4.9866 - 3.1435 * y;
    d[2] = -7.055 - 7.2165 * tanh(31.033 * (y + 2.1)) + 0.38541 * y;
    d[3] = -2.8915 - 2.1495 * y - 0.45006 * y2;
    d[4] = -1.297 - 0.13947 * y - 0.41197 * y2 - 0.10641 * y3;
    
    double sigma = eq12(x, d);

    // kinematic limits (for all species, diff and res)
    double Wdiff = 75.;
    double Lmax = y + 3;
    double fKL = 1 / (exp(Wdiff / (x - Lmax)) + 1);       

    return sigma * fKL;
}

// ---------------------- Electron Neutrino --------------------------------------------------------------

double HadronicInteractionKamae::SpectrumElectronNeutrino(double Esec, double Tp) const{
    if(Esec >= Tp) 
        return 0.;
    
    double s1 = ElectronNeutrinoND(Esec, Tp);
    double s2 = ElectronNeutrinoDiff(Esec, Tp);
    double s3 = ElectronNeutrino1232(Esec, Tp);
    double s4 = ElectronNeutrino1600(Esec, Tp);

    *writeOut << "nue" << "\t" << Tp / TeV << "\t" << Esec / GeV << "\t" << s1 << "\t" << s2 << "\t" << s3 << "\t" << s4 << "\n";
    return s1 + s2 + s3 + s4;
}

double HadronicInteractionKamae::ElectronNeutrinoND(double Esec, double Tp) const {
    std::vector<double> a;
    for(short i = 0; i < 9; i++)
        a.push_back(0.);
    
    double y = std::log10(Tp / TeV);
    double y2 = y * y;
    double y3 = y2 * y;
    double y4 = y3 * y;
    double z = y + 3.26; // helper variable
    double x = std::log10(Esec / GeV);

    double rescale = 1.;
    if(Tp < 7.81 * GeV) 
        rescale = 0.329 * exp(-247 * pow_integer<2>(z / (1. + 6.56 * z))) - 0.957 * y - 0.229 * y2;
    
    // table 6
    z = y + 3.31;
    a[0] = 0.0074087 + z * (2.9161 + z *(0.99061 + z * (-0.28694 + z * 0.038799)));
    a[1] = -3.248e-5 + 7.1944e-5 * exp(-0.21814 * (y + 3.4)) + 2.0467e-5 / (y + 4.1640) + 5.6954e-6 * y - 3.4105e-7 * y2;
    a[2] = -230.50 + 58.802 * y -  9.9393 * y2 + 1.2473 * y3 - 0.26322 * y4;
    a[3] = 0.45064 + 0.56930 * y + 0.012428 * y2 - 0.0070889 * y3;
    z = y + 3.32;
    a[4] = -0.011883 + z * (1.7992 + z * (3.5264 + z * (-1.7478 + z * (0.32077 - z * 0.017667))));
    z = y + 3.4;
    a[5] = -1.6238e-7 + 1.8116e-6 * exp(-0.30111 * z) + 9.6112e-5 / pow_integer<2>(y + 4.8229);
    a[6] = -261.30 - 43.351 * log10(0.35298 * z) + 70.925 / pow_integer<2>(y - 8.7147);
    a[7] = 184.45 - 1473.6 / (y + 6.8788) - 4.0536 * y2;
    a[8] = -0.24019 + 0.38504 * y + 0.0096869 * y2 - 0.0015046 * y3;

    double sigma = eq6(x, a);

    // kinematic limits (non diffractive - depends on particle type)
    double W_lo = 15;
    double W_hi = 42;
    double Lmin = -2.6;
    double Lmax = 0.98 * (y + 3);
    double fND_kl = 1 / (exp(W_lo * (Lmin - x)) + 1) / (exp(W_hi * (x - Lmax)) + 1);

    return rescale * sigma * fND_kl;
}

double HadronicInteractionKamae::ElectronNeutrinoDiff(double Esec, double Tp) const{
    std::vector<double> b; 
    for(short i = 0; i < 8; i++)
        b.push_back(0.);
    
    double y = std::log10(Tp / TeV);
    double y2 = y * y;
    double y3 = y2 * y;
    double y4 = y3 * y;
    double z = y + 3.26; // helper variable
    double x = std::log10(Esec / GeV);

    // table 6
    if(Tp > 11.05 * GeV) {
        b[0] = 53.809 * tanh(-0.41421 * (y + 2.2)) - 6.7538 * pow_integer<2>(y + 0.7601) + 0.0088080 * pow_integer<4>(y + 8.5075);
        z = y + 1.8901;
        b[1] = -50.211 + 55.131 * exp(1.3651 * pow_integer<2>(z / (1. + 4.444 * z)));
        b[2] = -17.231 + 0.041100 * tanh(7.9638 * (y + 1.9)) - 0.055449 * y + 2.5866e-4 * pow_integer<2>(y + 250.68);
        z = y + 1.8998;
        b[3] = 12.335 - 12.893 * exp(-1.4412 * pow_integer<2>(z / (1. + 5.5969 * z)));
    }
    z = y + 2.2;
    b[4] = 1.3558 + 0.46601 * (y + 2.95) + 0.052978 * z * z + 0.79575 * exp(-5.4007 * pow_integer<2>(z + 4.6121 * z * z));
    b[5] = 1.8756 - 0.42169 * tanh(1.6100 * (y + 1.9)) - 0.051026 * (y - 3.9573);
    b[6] = 1.5016 + 1.0118 * y - 0.072787 * y2 - 0.0038858 * y3 + 0.009365 * y4;
    b[7] = 4.9735 + 5.5674 * tanh(-0.36249 * (y + 2.1)) -0.20660 * (y - 2.8604);

    double sigma = eq9(x, b);

    // kinematic limits (for all species, diff and res)
    double Wdiff = 75.;
    double Lmax = y + 3;
    double fKL = 1 / (exp(Wdiff / (x - Lmax)) + 1);    

    return sigma * fKL;
}

double HadronicInteractionKamae::ElectronNeutrino1232(double Esec, double Tp) const {
    std::vector<double> c;
    for(short i = 0; i < 5; i++) 
        c.push_back(0.);
    
    double y = std::log10(Tp / TeV);
    double y2 = y * y;
    double z = y + 3.1282; // helper variable
    double x = std::log10(Esec / GeV);

    // table 6
    c[0] = 2.8290 * exp(-71.339 * pow_integer<2>(z / (1. + 0.48420 * z))) - 9.6339 - 15.733 / y + 0.52413 * y2;
    c[1] = -24.571 -15.831 * y - 2.12 * y2;
    c[2] = -5.9593 - 6.4695 * tanh(4.7225 * (y + 2.1)) + 0.50003 * y;
    c[3] = 0.26022 + 0.24545 * y;
    c[4] = 0.076498 + 0.061678 * y + 0.0040028 * y2;

    double sigma = eq12(x, c);

    // kinematic limits (for all species, diff and res)
    double Wdiff = 75.;
    double Lmax = y + 3;
    double fKL = 1 / (exp(Wdiff / (x - Lmax)) + 1);    

    return sigma * fKL;
}    

double HadronicInteractionKamae::ElectronNeutrino1600(double Esec, double Tp) const {
    std::vector<double> d;
    for(short i = 0; i < 5; i++) 
        d.push_back(0.);
    
    double y = std::log10(Tp / TeV);
    double y2 = y * y;
    double y3 = y2 * y;
    double z = y + 2.9509; // helper variable
    double x = std::log10(Esec / GeV);

    // table 6
    d[0] = 1.7951 * exp(-57.26 * pow_integer<2>(z / (1. + 1.4101 * z))) - 0.58604 - 0.23868 * y;
    d[1] = -2.6395 -1.5105 * y + 0.22174 * y2;
    d[2] = -7.0512 - 7.1970 * tanh(31.074 * (y + 2.1)) + 0.39007 * y;
    d[3] = -1.4271 - 1.0399 * y - 0.24179 * y2;
    d[4] = 0.74875 + 0.63616 * y + 0.17396 * y2 + 0.017636 * y3;

    double sigma = eq12(x, d);

    // kinematic limits (for all species, diff and res)
    double Wdiff = 75.;
    double Lmax = y + 3;
    double fKL = 1 / (exp(Wdiff / (x - Lmax)) + 1);    

    return sigma * fKL;
}        

// ---------------------- Anti Electron Neutrino --------------------------------------------------------------

double HadronicInteractionKamae::SpectrumAntiElectronNeutrino(double Esec, double Tp) const {
    
    // the formulars for the anti electron neutrino allow a infinit expresseion for the diffractive process. This happens for example
    // when Tp = 1 TeV and Esec = 9.65 GeV is used. Therefore we assume that the spectra of the electron neutrino and the anti electron 
    // neutrino are the same. 
    *writeOut << "Bar";
    // return SpectrumElectronNeutrino(Esec, Tp); 

    if (Esec >= Tp) 
        return 0.;
    
    double s1 = AntiElectronNeutrinoND(Esec, Tp);
    double s2 = AntiElectronNeutrinoDiff(Esec, Tp);
    double s3 = AntiElectronNeutrino1232(Esec, Tp);
    double s4 = AntiElectronNeutrino1600(Esec, Tp);


    *writeOut << "nueBar" << "\t" << Tp / TeV << "\t" << Esec / GeV << "\t" << s1 << "\t" << s2 << "\t" << s3 << "\t" << s4 << "\n";    
    return s1 + s2 + s3 + s4;
}

double HadronicInteractionKamae::AntiElectronNeutrinoND(double Esec, double Tp) const {
    std::vector<double> a;
    for(short i = 0; i < 9; i++)
        a.push_back(0.);
    
    double y = std::log10(Tp / TeV);
    double y2 = y * y;
    double y3 = y2 * y;
    double y4 = y3 * y;
    double z = y + 3.27; // helper variable
    double x = std::log10(Esec / GeV);

    double rescale = 1.;
    if(Tp > 15.6 * GeV) 
        rescale = 2.67 * exp(-45.7 * pow_integer<2>(z / (1. + 6.59 * z))) - 0.301 * y - 0.208 * y2;
    
    z = y + 3.31;
    a[0] = 0.0013113 + z * (0.36538 + z * (1.5178 + z *(-0.20668 + z * 0.024255)));
    a[1] = -4.7833e-6 + 4.5837e-5 * exp( -0.42980 * (y + 3.4)) + 6.1559e-6 / (y + 4.1731) + 1.1928e-6 * y;
    a[2] = -245.22 + 73.223 * y - 19.652 * y2 + 0.083138 * y3 + 0.71561 * y4;
    a[3] = 0.45232 + 0.52934 * y + 0.010078 * y2 - 0.0017092 * y3;
    z = y + 3.32;
    a[4] = -0.0025734 + z * (0.38424 + z * (1.5517 + z * (0.17336 + z * (-0.17160 + z * 0.021059))));
    a[5] = 4.7673e-5 + 5.4936e-5 * log10(0.0067905 * (y + 4.3)) + 0.00020740 / (y + 4.9772);
    a[6] = -270.30 - 114.47 * log10(0.34352 * (y + 3.4)) + 80.085 * y - 7.924 * y2;
    a[7] = 3272.9 - 2.9161e5 / (y + 87.847) - 6.2330 * y2;
    a[8] = -0.17787 + 0.36771 * y - 0.025397 * y2 + 0.0019238 * y3 + 0.0032725 * y4;

    double sigma = eq6(x, a);

    // kinematic limits (non diffractive - depends on particle type)
    double Lmax = 0.98 * (y + 3);
    double W_lo = 15;
    double W_hi = 40;
    double Lmin = -2.6;
    double fND_kl = 1 / (exp(W_lo * (Lmin - x)) + 1) / (exp(W_hi * (x - Lmax)) + 1);

    return rescale * sigma * fND_kl;
}

double HadronicInteractionKamae::AntiElectronNeutrinoDiff(double Esec, double Tp) const {
    std::vector<double> b; 
    for(short i = 0; i < 8; i++)
        b.push_back(0.);
    
    double y = std::log10(Tp / TeV);
    double y2 = y * y;
    double y3 = y2 * y;
    double y4 = y3 * y;
    double z = y + 3.26; // helper variable
    double x = std::log10(Esec / GeV);

    // table 7
    if(Tp > 11.05 * GeV) {
        b[0] = 41.307 * tanh(-0.37411 * (y + 2.2)) - 4.1223 * pow_integer<2>(y + 0.55505) + 0.0042652 * pow_integer<4>(y + 9.2685);
        z = y + 1.9196;
        b[1] = -132.50 + 142.12 * exp(-8.0289 * pow_integer<2>(z / (1. + 11.530 * z))); 
        b[2] = -17.223 + 0.011285 * tanh(-69.746 * (y + 1.9)) - 0.048233 * y + 2.5881e-4 * pow_integer<2>(y + 250.77);
        z = y + 1.9292;
        b[3] = 8.1991 - 9.6437 * exp(-45.261 * pow_integer<2>(z / (1. + 16.682 * z))); 
    }
    z = pow_integer<2>(y + 2.95);
    b[4] = 0.55919 + z * (0.36647 + z * 0.056194); 
    z = y + 2.2;
    b[4] += 0.49957 * exp(-5.5317 * pow_integer<2>(z + 0.43867 * z * z));
    b[5] = 1.2544 - 0.52362 * tanh(2.7638 * (y + 1.9)) + 0.055837 * (y - 17.638);
    b[6] = 1.4788 + 1.0278 * y - 0.092852 * y2 - 0.0062734 * y3 + 0.01192 * y4;
    b[7] = 5.1651 + 5.7398 * tanh(-0.37356 * (y + 2.1)) - 0.22234 * (y - 2.7889);

    // std::cout << b[6] - 1. / b[7] <<"\t" << b[6] << " , " << b[7] << " , " << y << " , " << y2 << " , " << y3 << " , " << y4 << "\n";
    double sigma = eq9(x, b);

    // kinematic limits (for all species, diff and res)
    double Wdiff = 75.;
    double Lmax = y + 3;
    double fKL = 1 / (exp(Wdiff / (x - Lmax)) + 1);    

    return sigma * fKL;
}

double HadronicInteractionKamae::AntiElectronNeutrino1232(double Esec, double Tp) const{
    return 0.;
}

double HadronicInteractionKamae::AntiElectronNeutrino1600(double Esec, double Tp) const {
    std::vector<double> d;
    for(short i = 0; i < 5; i++) 
        d.push_back(0.);
    
    double y = std::log10(Tp / TeV);
    double y2 = y * y;
    double y3 = y2 * y;
    double z = y + 2.9537; // helper variable
    double x = std::log10(Esec / GeV);

    // table 7
    d[0] = 0.36459 * exp(-58.210 * pow_integer<2>(z / (1. + 1.4320 * z))) - 0.11283 - 0.046244 * y;
    d[1] = -9.5066 - 5.4655 * y - 0.31769 * y2;
    d[2] = -7.1831 - 7.1551 * tanh(30.354 * (y + 2.1)) + 0.33757 * y;
    d[3] = 2.7938 + 16992 * y + 0.20161 * y2;
    d[4] = 0.61878 + 0.62371 * y + 0.18913 * y2 + 0.019118 * y3;

    double sigma = eq12(x, d);

    // kinematic limits (for all species, diff and res)
    double Wdiff = 75.;
    double Lmax = y + 3;
    double fKL = 1 / (exp(Wdiff / (x - Lmax)) + 1);    

    return sigma * fKL;
}

// ---------------------- Muon Neutrino --------------------------------------------------------------

double HadronicInteractionKamae::SpectrumMuonNeutrino(double Esec, double Tp) const {
    if (Esec >= Tp)
        return 0.;
    
    double s1 = MuonNeutrinoND(Esec, Tp);
    double s2 = MuonNeutrinoDiff(Esec, Tp);
    double s3 = MuonNeutrino1232(Esec, Tp);
    double s4 = MuonNeutrino1600(Esec, Tp);

    *writeOut << "numu" << "\t" << Tp / TeV << "\t" << Esec / GeV << "\t" << s1 << "\t" << s2 << "\t" << s3 << "\t" << s4 << "\n";
    return s1 + s2 + s3 + s4;
}

double HadronicInteractionKamae::MuonNeutrinoND(double Esec, double Tp) const {
    std::vector<double> a;
    for(short i = 0; i < 9; i++)
        a.push_back(0.);
    
    double y = std::log10(Tp / TeV);
    double y2 = y * y;
    double y3 = y2 * y;
    double y4 = y3 * y;
    double z = y + 3.25; // helper variable
    double x = std::log10(Esec / GeV);

    double rescale = 1.;
    if(Tp > 15.6 * GeV) 
        rescale = 2.23 * exp(-93.4 * pow_integer<2>(z / (1. + 8.38 * z))) - 0.376 * y - 0.121 * y2;
    
    // table 8 
    z = y + 3.3;
    a[0] = z * (-0.63611  + z * (9.9015 + z * (-4.5897 + z * (0.91778 + z * - 0.060724))));
    a[1] = 6.87e-6 - 2.8245e-6 * y + 7.6032e-7 * y2 - 3.2953e-7 * y3 + 7.4292e-8 * y4;
    a[2] = -240.46 + 58.405 * y - 9.8556 * y2 + 3.1401 * y3 - 0.88932 * y4;
    a[3] = 0.49935 + 0.60919 * y + 0.0024963 * y2 - 0.0099910 * y3;
    z = y + 3.32;
    a[4] = z * (2.5094 + z * (4.1350 + z * (-0.89534 + z * (-0.0027577 + z * 0.014511))));
    a[5] = 8.2046e-7 + 1.4085e-6 * log10(0.016793 * (y + 4.3)) + 0.0001334 / pow_integer<2>(y + 4.7136);
    a[6] = -267.55 - 0.21018 * log10(0.35217 * (y + 3.9)) + 69.586 * y - 9.9930 * y2;
    a[7] = 2742.8 + 222.01 * log10(9.7401 * (y + 3.9)) - 4772.5 / (y + 19.773) - 6.1001 * y2;
    a[8] = -0.11857 + 0.39072 * y - 0.037813 * y2 + 0.0022265 * y3 + 0.0046931 * y4;

    double sigma = eq6(x, a);

    // kinematic limits (non diffractive - depends on particle type, table 2)
    double Lmax = 0.94 * (y + 3);
    double W_lo = 20;
    double W_hi = 45;
    double Lmin = -2.6;
    double fND_kl = 1 / (exp(W_lo * (Lmin - x)) + 1) / (exp(W_hi * (x - Lmax)) + 1);

    return rescale * sigma * fND_kl;
}

double HadronicInteractionKamae::MuonNeutrinoDiff(double Esec, double Tp) const {
    std::vector<double> b; 
    for(short i = 0; i < 8; i++)
        b.push_back(0.);
    
    double y = std::log10(Tp / TeV);
    double y2 = y * y;
    double y3 = y2 * y;
    double y4 = y3 * y;
    double z = y + 3.26; // helper variable
    double x = std::log10(Esec / GeV);

    // table 8
    if(Tp > 11.05 * GeV) {
        b[0] = 64.682 * tanh(-0.34313 * (y + 2.2)) - 5.5955 * pow_integer<2>(y + 0.44754) + 0.0050117 * pow_integer<4>(y + 9.9165);
        z = y + 2.3066;
        b[1] = -7.6016 + 3.0427e4 * exp(-1.0134e4 * pow_integer<2>(z / (1. + 41.612 * z)));
        b[2] = -1.4978 - 0.58163 * tanh(-0.36488 * (y + 1.9)) + 0.031825 * (y + 2.8097) + 0.022796 * pow_integer<2>(y - 1.8861);
        z = y + 3.8835;
        b[3] = -0.0061483 - 65.799 * exp(- 4.8239 * pow_integer<2>(z / (1. + 0.53343 * z)));
    }
    z = y + 2.95;
    b[4] = 2.8009 + z * z * (0.35341 - 0.0039779 * z * z);
    z = y + 2.2;
    b[4] +=  1.3012 * exp(-10.592 * pow_integer<2>(z - 0.19149 * z * z));
    b[5] = 1.8016 - 0.69847 * tanh(2.8627 * (y + 1.9)) - 0.015722 * (y - 45.4);
    b[6] = 1.4617 + 1.0167 * y - 0.078617 * y2 - 0.0038336 * y3 + 0.010141 * y4;
    b[7] = 3.5599 + 4.0041 * tanh(-0.41889 * (y + 2.1)) - 0.18182 * (y - 2.4209);

    double sigma = eq9(x, b);

    // kinematic limits (for all species, diff and res)
    double Wdiff = 75.;
    double Lmax = y + 3;
    double fKL = 1 / (exp(Wdiff / (x - Lmax)) + 1);    

    return sigma * fKL;
}

double HadronicInteractionKamae::MuonNeutrino1232(double Esec, double Tp) const {
    std::vector<double> c;
    for(short i = 0; i < 5; i++) 
        c.push_back(0.);
    
    double y = std::log10(Tp / TeV);
    double y2 = y * y;
    double z = y + 3.1278; // helper variable
    double x = std::log10(Esec / GeV); 

    // table 8
    c[0] = 3.6052 * exp(-60.914 * pow_integer<2>(z / (1. - 0.19497 * z))) - 0.92514 + 2.1315 / y + 0.23548 * y2;
    c[1] = 95.310 + 70.497 * y + 13.636 * y2;
    c[2] = -6.2158 - 6.2939 * tanh(21.592 * (y + 2.1)) + 0.37440 * y;
    c[3] = 2.7485 + 1.1692 * y;
    c[4] = -2.7568 - 1.8461 * y - 0.31376 * y2;

    double sigma = eq12(x, c);

    // kinematic limits (for all species, diff and res)
    double Wdiff = 75.;
    double Lmax = y + 3;
    double fKL = 1 / (exp(Wdiff / (x - Lmax)) + 1);    

    return sigma * fKL;    
}

double HadronicInteractionKamae::MuonNeutrino1600(double Esec, double Tp) const {
    std::vector<double> d;
    for(short i = 0; i < 5; i++) 
        d.push_back(0.);
    
    double y = std::log10(Tp / TeV);
    double y2 = y * y;
    double y3 = y2 * y;
    double z = y + 2.9509; // helper variable
    double x = std::log10(Esec / GeV);

    // table 8
    d[0] = 2.5489 * exp(-58.488 * pow_integer<2>(z / (1. + 1.3154 * z))) - 0.83039 - 0.34412 * y;
    d[1] = 88.173 + 65.148 * y + 12.585 * y2;
    d[2] = -7.0962 - 7.1690 * tanh(30.89 * (y + 2.1)) +  0.38032 * y;
    d[3] = -4.1440 - 3.2717 * y - 0.70537 * y2;
    d[4] = 2.2624 + 1.1806 * y - 0.0043450 * y2 - 0.04302 * y3;

    double sigma = eq12(x, d);

    // kinematic limits (for all species, diff and res)
    double Wdiff = 75.;
    double Lmax = y + 3;
    double fKL = 1 / (exp(Wdiff / (x - Lmax)) + 1);    

    return sigma * fKL;  
}

// ---------------------- Muon Neutrino --------------------------------------------------------------

double HadronicInteractionKamae::SpectrumAntiMuonNeutrino(double Esec, double Tp) const {
    if(Esec >= Tp) 
        return 0.;
    
    double s1 = AntiMuonNeutrinoND(Esec, Tp);
    double s2 = AntiMuonNeutrinoDiff(Esec, Tp);
    double s3 = AntiMuonNeutrino1232(Esec, Tp);
    double s4 = AntiMuonNeutrino1600(Esec, Tp);

    *writeOut << "numuBar" << "\t" << Tp / TeV << "\t" << Esec / GeV << "\t" << s1 << "\t" << s2 << "\t" << s3 << "\t" << s4 << "\n";
    return s1 + s2 + s3 + s4;
}

double HadronicInteractionKamae::AntiMuonNeutrinoND(double Esec, double Tp) const {
    std::vector<double> a;
    for(short i = 0; i < 9; i++)
        a.push_back(0.);
    
    double y = std::log10(Tp / TeV);
    double y2 = y * y;
    double y3 = y2 * y;
    double y4 = y3 * y;
    double z = y + 3.25; // helper variable
    double x = std::log10(Esec / GeV);

    double rescale = 1.;
    if(Tp > 15.6 * GeV) 
        rescale = 2.56 * exp(-107. * pow_integer<2>(z / (1. + 8.34* z))) - 0.385 * y - 0.125 * y2;    
    
    // table 9
    z = y + 3.3;
    a[0] = z * (-1.5243 + z * (10.107 + z * (-4.3126 + z * (0.80081 - z * 0.048724))));
    a[1] = -2.6297e-5 + 9.3858e-5 * exp(-0.32384 * (y + 3.4)) + 7.7821e-6 / (y + 4.0560) + 7.6149e-6 * y - 8.4091e-6 * y2;
    a[2] = -223.62 + 59.374 * y - 5.7356 * y2 + 1.9815 * y3 - 1.0478 * y4;
    a[3] = 0.50807 + 0.60221 * y + 0.0034120 * y2 - 0.011139 * y3;
    z = y + 3.32;
    a[4] = z * (2.6483 + z * ( 4.4585 + z * (-1.2744 + z * (-0.11659 + z * 0.0030477))));
    a[5] = 9.1101e-7 + 1.388e-6 * log10(0.016998 * (y + 4.3)) + 1.2744e-4 / pow_integer<2>(y + 4.7707);
    a[6] = -272.11 - 53.477 * log10(0.35531 * (y + 3.9)) + 56.041 / pow_integer<2>(y - 6.0876);
    a[7] = 6431.8 + 893.92 * log10(5.713e-9 *(y + 3.9)) + 2103.6 / (y + 5.6740) - 6.1125 * y2;
    a[8] = -0.11120 + 0.38144 * y - 0.040128 * y2 + 0.0047484 * y3 + 0.0054707 * y4;

    double sigma = eq6(x, a);

    // kinematic limits (non diffractive - depends on particle type, table 2)
    double Lmax = 0.98 * (y + 3);
    double W_lo = 15;
    double W_hi = 40;
    double Lmin = -2.6;
    double fND_kl = 1 / (exp(W_lo * (Lmin - x)) + 1) / (exp(W_hi * (x - Lmax)) + 1);

    return rescale * sigma * fND_kl;
}

double HadronicInteractionKamae::AntiMuonNeutrinoDiff(double Esec, double Tp) const {
    std::vector<double> b; 
    for(short i = 0; i < 8; i++)
        b.push_back(0.);
    
    double y = std::log10(Tp / TeV);
    double y2 = y * y;
    double y3 = y2 * y;
    double y4 = y3 * y;
    double z = y + 3.26; // helper variable
    double x = std::log10(Esec / GeV);

    // table 9
    if(Tp > 11.05 * GeV) {
        b[0] = 70.43 * tanh(-0.35816 * (y + 2.2)) - 6.6796 * pow_integer<2>(y + 0.52273) + 0.0065659 * pow_integer<4>(y + 9.5266);
        z = y + 2.219;
        b[1] = -8.1145 + 7686 * exp(4.4046e4 * pow_integer<2>(z / (1. + 81.105 *z)));
        b[2] = -1.3095 + 0.07127 * tanh(-0.0075463 * (y + 1.9)) + 0.067759 * (y + 5.3433) - 0.0044205 * pow_integer<2>(y - 1.8683);
        z = y + 2.8363;
        b[3] = 0.082149 - 2190.1 * exp(-533.75 * pow_integer<2>(z / (1. + 7.0976 * z)));
    }    
    z = pow_integer<2>(y + 2.95);
    b[4] = 2.754 + z * (0.33859 - 0.0034274 * z);
    z = y + 2.2;
    b[4] += 1.1679 * exp(-10.408 * pow_integer<2>(z - 0.18922 * z * z));
    b[5] = 2.1817 - 0.59584 * tanh(2.7054 * (y + 1.9)) - 0.010909 * (y - 14.9);
    b[6] = 1.4591 + 1.0275 * y - 0.074949 * y2 - 0.0060396 * y3 + 0.0097568 * y4;
    b[7] = 3.7609 + 4.2843 * tanh(-0.37148 * (y + 2.1)) - 0.16479 * (y - 2.7653);

    double sigma = eq9(x, b);

    // kinematic limits (for all species, diff and res)
    double Wdiff = 75.;
    double Lmax = y + 3;
    double fKL = 1 / (exp(Wdiff / (x - Lmax)) + 1);    

    return sigma * fKL; 
}

double HadronicInteractionKamae::AntiMuonNeutrino1232(double Esec, double Tp) const {
    std::vector<double> c;
    for(short i = 0; i < 5; i++) 
        c.push_back(0.);
    
    double y = std::log10(Tp / TeV);
    double y2 = y * y;
    double z = y + 3.1250; // helper variable
    double x = std::log10(Esec / GeV); 

    // table 9
    c[0] = 2.8262 * exp(-62.894 * pow_integer<2>(z / (1. - 0.47567 * z))) + 5.6845 + 13.409 / y - 0.097296 * y2;
    c[1] = 16.721 + 11.750 * y + 2.4637 * y2;
    c[2] = -6.0557 - 6.3378 * tanh(21.984 * (y + 2.1)) + 0.43173 * y;
    c[3] = 0.37009 + 0.27706 * y;
    c[4] = 0.047507 + 0.061570 * y + 0.0070117 * y2;

    double sigma = eq12(x, c);   

    // kinematic limits (for all species, diff and res)
    double Wdiff = 75.;
    double Lmax = y + 3;
    double fKL = 1 / (exp(Wdiff / (x - Lmax)) + 1);    

    return sigma * fKL; 
}

double HadronicInteractionKamae::AntiMuonNeutrino1600(double Esec, double Tp) const {
    std::vector<double> d;
    for(short i = 0; i < 5; i++) 
        d.push_back(0.);
    
    double y = std::log10(Tp / TeV);
    double y2 = y * y;
    double y3 = y2 * y;
    double z = y + 2.9492; // helper variable
    double x = std::log10(Esec / GeV);     

    // table 9
    d[0] = 2.24 * exp(-57.159 * pow_integer<2>(z / (1. + 1.2994 * z))) - 0.66521 - 0.27554 * y;
    d[1] = -7.065 - 4.2773 * y - 0.17648 * y2;
    d[2] = -7.0410 - 7.1977 * tanh(31.095 * (y + 2.1)) + 0.40238 * y;
    d[3] = -1.2354 - 0.87581 * y - 0.20829 * y2;
    d[4] = -0.11395 + 0.34418 * y + 0.27103 * y2 + 0.050248 * y3;

    double sigma = eq12(x, d);   

    // kinematic limits (for all species, diff and res)
    double Wdiff = 75.;
    double Lmax = y + 3;
    double fKL = 1 / (exp(Wdiff / (x - Lmax)) + 1);    

    return sigma * fKL; 
}