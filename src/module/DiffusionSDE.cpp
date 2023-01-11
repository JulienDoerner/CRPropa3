#include "crpropa/module/DiffusionSDE.h"



using namespace crpropa;

// Defining Cash-Karp coefficients
const double a[] = { 0., 0., 0., 0., 0., 0., 1. / 5., 0., 0., 0., 0.,
		0., 3. / 40., 9. / 40., 0., 0., 0., 0., 3. / 10., -9. / 10., 6. / 5.,
		0., 0., 0., -11. / 54., 5. / 2., -70. / 27., 35. / 27., 0., 0., 1631.
				/ 55296., 175. / 512., 575. / 13824., 44275. / 110592., 253.
				/ 4096., 0. };

const double b[] = { 37. / 378., 0, 250. / 621., 125. / 594., 0., 512.
		/ 1771. };

const double bs[] = { 2825. / 27648., 0., 18575. / 48384., 13525.
		/ 55296., 277. / 14336., 1. / 4. };

// turbulence dependend values for the diffusion coefficient
const std::vector<double> bB = {0.067, 0.108, 0.174, 0.28, 0.452, 0.727, 1.172, 1.887, 3.022, 3.834, 4.865, 6.172, 7.832};

const std::vector<double> kappa_para_0 = {7.62027106e+32, 2.48212330e+32, 8.33535627e+31, 2.96949725e+31,
 1.14828611e+31, 5.08992373e+30, 2.66923660e+30, 1.73875117e+30,
 1.35040055e+30, 1.25957468e+30, 1.17795162e+30, 1.11661816e+30,
 1.08112229e+30};
const std::vector<double> gamma_para = {0.563, 0.585, 0.649, 0.689, 0.73, 0.768, 0.829, 0.888, 0.913, 0.923, 0.911}; // Daten Stand: 19/02/2021

const std::vector<double> kappa_perp_0 = {2.05687400e+27, 3.76269310e+27, 7.09804916e+27, 1.44494461e+28,
 3.28214232e+28, 7.67791668e+28, 1.73937204e+29, 3.63945258e+29,
 6.27195917e+29, 7.56419514e+29, 8.54756784e+29, 9.25339363e+29,
 9.65696170e+29};
const std::vector<double> gamma_perp = {0.343, 0.378, 0.418, 0.444, 0.445, 0.481, 0.586, 0.756, 0.901, 0.926, 0.932}; // Daten Stand 19/02/2021


DiffusionSDE::DiffusionSDE(ref_ptr<MagneticField> magneticField, double tolerance,
				 double minStep, double maxStep, double epsilon) :
	minStep(0)
{
  	setMagneticField(magneticField);
  	setMaximumStep(maxStep);
  	setMinimumStep(minStep);
  	setTolerance(tolerance);
	diffusionTensor = new QLTDiffusion(epsilon);
	}

DiffusionSDE::DiffusionSDE(ref_ptr<MagneticField> magneticField, ref_ptr<DiffusionTensor> diffusionTensor, double tolerance, double minStep, double maxStep) : 
	magneticField(magneticField), diffusionTensor(diffusionTensor), minStep(0) 
{
	setMaximumStep(maxStep);
	setMinimumStep(minStep);
	setTolerance(tolerance);
}

DiffusionSDE::DiffusionSDE(ref_ptr<MagneticField> magneticField, ref_ptr<AdvectionField> advectionField, double tolerance, double minStep, double maxStep, double epsilon) :
  	minStep(0)
{
	setMagneticField(magneticField);
	setAdvectionField(advectionField);
	setMaximumStep(maxStep);
	setMinimumStep(minStep);
	setTolerance(tolerance);
	diffusionTensor = new QLTDiffusion(epsilon);
}

DiffusionSDE::DiffusionSDE(ref_ptr<MagneticField> magneticField, ref_ptr<AdvectionField> advectionField, ref_ptr<DiffusionTensor> diffusionTensor, double tolerance, double minStep, double maxStep) :
	magneticField(magneticField), advectionField(advectionField), diffusionTensor(diffusionTensor), minStep(0) 
{
	setMaximumStep(maxStep);
	setMinimumStep(minStep);
	setTolerance(tolerance);
}
	

void DiffusionSDE::process(Candidate *candidate) const {

    // save the new previous particle state

	ParticleState &current = candidate->current;
	candidate->previous = current;

	double h = clip(candidate->getNextStep(), minStep, maxStep) / c_light;
	Vector3d PosIn = current.getPosition();
	Vector3d DirIn = current.getDirection();

    // rectilinear propagation for neutral particles
    // If an advection field is provided the drift is also included
	if (current.getCharge() == 0) {
		Vector3d dir = current.getDirection();
		Vector3d Pos = current.getPosition();

		Vector3d LinProp(0.);
		if (advectionField){
			driftStep(Pos, LinProp, h);
		}

		current.setPosition(Pos + LinProp + dir*h*c_light);
		candidate->setCurrentStep(h * c_light);
		candidate->setNextStep(maxStep);
		return;
	}

	double z = candidate->getRedshift();
	double rig = current.getEnergy() / current.getCharge();
	

    // Calculate the Diffusion tensor
	double BTensor[] = {0., 0., 0., 0., 0., 0., 0., 0., 0.};

	Vector3d diffCoeff = diffusionTensor -> getDiffusionKoefficent(candidate);

	BTensor[0] = pow(2 * diffCoeff.x, 0.5); // parallel component
    BTensor[4] = pow(2 * diffCoeff.y, 0.5);
    BTensor[8] = pow(2 * diffCoeff.z, 0.5);

    // Generate random numbers
	double eta[] = {0., 0., 0.};
	for(size_t i=0; i < 3; i++) {
	  	eta[i] =  Random::instance().randNorm();
	}

	double TStep = BTensor[0] * eta[0];
	double NStep = BTensor[4] * eta[1];
	double BStep = BTensor[8] * eta[2];

	Vector3d TVec(0.);
	Vector3d NVec(0.);
	Vector3d BVec(0.);

	Vector3d DirOut = Vector3d(0.);


	double propTime = TStep * sqrt(h) / c_light;
	size_t counter = 0;
	double r=42.; //arbitrary number larger than one

	do {
		Vector3d PosOut = Vector3d(0.);
		Vector3d PosErr = Vector3d(0.);
	  	tryStep(PosIn, PosOut, PosErr, z, propTime);
	    // calculate the relative position error r and the next time step h
	  	r = PosErr.getR() / tolerance;
	  	propTime *= 0.5;
		counter += 1;

    // Check for better break condition
	} while (r > 1 && fabs(propTime) >= minStep/c_light);


	size_t stepNumber = pow(2, counter-1);
	double allowedTime = TStep * sqrt(h) / c_light / stepNumber;
	Vector3d Start = PosIn;
	Vector3d PosOut = Vector3d(0.);
	Vector3d PosErr = Vector3d(0.);
	for (size_t j=0; j<stepNumber; j++) {
		tryStep(Start, PosOut, PosErr, z, allowedTime);
		Start = PosOut;
	}

    // Normalize the tangent vector
	TVec = (PosOut-PosIn).getUnitVector();
    // Exception: If the magnetic field vanishes: Use only advection.
    // If an advection field is not provided --> rectilinear propagation.
	double tTest = TVec.getR();
	if (tTest != tTest) {
	  	Vector3d dir = current.getDirection();
		Vector3d Pos = current.getPosition();
		Vector3d LinProp(0.);
		if (advectionField){
			driftStep(Pos, LinProp, h);
			current.setPosition(Pos + LinProp);
	 		candidate->setCurrentStep(h*c_light);
	  		double newStep = 5*h*c_light;
			newStep = clip(newStep, minStep, maxStep);
	  		candidate->setNextStep(newStep);
	  		return;
		}
		current.setPosition(Pos + dir*h*c_light);
	 	candidate->setCurrentStep(h*c_light);
		double newStep = 5*h*c_light;
		newStep = clip(newStep, minStep, maxStep);
	  	candidate->setNextStep(newStep);
	  	return;
	}

    // Choose a random perpendicular vector as the Normal-vector.
    // Prevent 'nan's in the NVec-vector in the case of <TVec, NVec> = 0.
	while (NVec.getR()==0.){
	  	Vector3d RandomVector = Random::instance().randVector();
	  	NVec = TVec.cross( RandomVector );
	}
	NVec = NVec.getUnitVector();

    // Calculate the Binormal-vector
	BVec = (TVec.cross(NVec)).getUnitVector();

    // Calculate the advection step
	Vector3d LinProp(0.);
	if (advectionField){
		driftStep(PosIn, LinProp, h);
	}

    // Integration of the SDE with a Mayorama-Euler-method
	Vector3d PO = PosOut + LinProp + (NVec * NStep + BVec * BStep) * sqrt(h) ;

    // Throw error message if something went wrong with propagation.
    // Deactivate candidate.
	bool NaN = std::isnan(PO.getR());
	if (NaN == true){
		  candidate->setActive(false);
		  KISS_LOG_WARNING
			<< "\nCandidate with 'nan'-position occured: \n"
		 	<< "position = " << PO << "\n"
		  	<< "PosIn = " << PosIn << "\n"
		  	<< "TVec = " << TVec << "\n"
		  	<< "TStep = " << std::abs(TStep) << "\n"
		  	<< "NVec = " << NVec << "\n"
		  	<< "NStep = " << NStep << "\n"
		  	<< "BVec = " << BVec << "\n"
		  	<< "BStep = " << BStep << "\n"
			<< "Candidate is deactivated!\n";
		  return;
	}

	//DirOut = (PO - PosIn - LinProp).getUnitVector(); //Advection does not change the momentum vector
	// Random direction around the tangential direction accounts for the pitch angle average.
	DirOut = Random::instance().randConeVector(TVec, M_PI/2.);
	current.setPosition(PO);
	current.setDirection(DirOut);
	candidate->setCurrentStep(h * c_light);

	double nextStep;
	if (stepNumber>1){
		nextStep = h*pow(stepNumber, -2.)*c_light;
	}
	else {
		nextStep = 4 * h*c_light;
	}

	candidate->setNextStep(nextStep);

	// Debugging and Testing
	// Delete comments if additional information should be stored in candidate
	// This property "arcLength" can be interpreted as the effective arclength
	// of the propagation along a magnetic field line.

/*
	const std::string AL = "arcLength";
	if (candidate->hasProperty(AL) == false){
	  double arcLen = (TStep + NStep + BStep) * sqrt(h);
	  candidate->setProperty(AL, arcLen);
	  return;
	}
	else {
	  double arcLen = candidate->getProperty(AL);
	  arcLen += (TStep + NStep + BStep) * sqrt(h);
	  candidate->setProperty(AL, arcLen);
	}
*/

}


void DiffusionSDE::tryStep(const Vector3d &PosIn, Vector3d &POut, Vector3d &PosErr,double z, double propStep) const {

	Vector3d k[] = {Vector3d(0.),Vector3d(0.),Vector3d(0.),Vector3d(0.),Vector3d(0.),Vector3d(0.)};
	POut = PosIn;
	//calculate the sum k_i * b_i
	for (size_t i = 0; i < 6; i++) {

		Vector3d y_n = PosIn;
		for (size_t j = 0; j < i; j++)
		  y_n += k[j] * a[i * 6 + j] * propStep;

		// update k_i = direction of the regular magnetic mean field
		Vector3d BField = getMagneticFieldAtPosition(y_n, z);

		k[i] = BField.getUnitVector() * c_light;

		POut += k[i] * b[i] * propStep;
		PosErr +=  (k[i] * (b[i] - bs[i])) * propStep / kpc;

	}
}

void DiffusionSDE::driftStep(const Vector3d &pos, Vector3d &linProp, double h) const {
	Vector3d advField = getAdvectionFieldAtPosition(pos);
	linProp += advField * h;
	return;
}

void DiffusionSDE::setMinimumStep(double min) {
	if (min < 0)
		throw std::runtime_error("DiffusionSDE: minStep < 0 ");
	if (min > maxStep)
		throw std::runtime_error("DiffusionSDE: minStep > maxStep");
	minStep = min;
}

void DiffusionSDE::setMaximumStep(double max) {
	if (max < minStep)
		throw std::runtime_error("DiffusionSDE: maxStep < minStep");
	maxStep = max;
}


void DiffusionSDE::setTolerance(double tol) {
	if ((tol > 1) or (tol < 0))
		throw std::runtime_error(
				"DiffusionSDE: tolerance error not in range 0-1");
	tolerance = tol;
}

void DiffusionSDE::setMagneticField(ref_ptr<MagneticField> f) {
	magneticField = f;
}

void DiffusionSDE::setAdvectionField(ref_ptr<AdvectionField> f) {
	// std::cout << "output from setAdvectionField with field "<< f->getField(Vector3d(0.)) << "\n";
	advectionField = f;
}

void DiffusionSDE::setDiffusionTensor(ref_ptr<DiffusionTensor> t) {
	diffusionTensor = t;
}

double DiffusionSDE::getMinimumStep() const {
	return minStep;
}

double DiffusionSDE::getMaximumStep() const {
	return maxStep;
}

double DiffusionSDE::getTolerance() const {
	return tolerance;
}


ref_ptr<MagneticField> DiffusionSDE::getMagneticField() const {
	return magneticField;
}

Vector3d DiffusionSDE::getMagneticFieldAtPosition(Vector3d pos, double z) const {
	Vector3d B(0, 0, 0);
	try {
		// check if field is valid and use the field vector at the
		// position pos with the redshift z
		if (magneticField.valid())
			B = magneticField->getField(pos, z);
	}
	catch (std::exception &e) {
		KISS_LOG_ERROR 	<< "DiffusionSDE: Exception in DiffusionSDE::getMagneticFieldAtPosition.\n"
				<< e.what();
	}	
	return B;
}

ref_ptr<AdvectionField> DiffusionSDE::getAdvectionField() const {
	return advectionField;
}

Vector3d DiffusionSDE::getAdvectionFieldAtPosition(Vector3d pos) const {
	Vector3d AdvField(0.);
	try {
		// check if field is valid and use the field vector at the
		// position pos
		if (advectionField.valid())
			AdvField = advectionField->getField(pos);
	}
	catch (std::exception &e) {
		KISS_LOG_ERROR 	<< "DiffusionSDE: Exception in DiffusionSDE::getAdvectionFieldAtPosition.\n"
				<< e.what();
	}
	return AdvField;
}

std::string DiffusionSDE::getDescription() const {
	std::stringstream s;
	s << "Diffusion SDE Model with following specification: \n";
	s << "minStep: " << minStep / kpc  << " kpc, ";
	s << "maxStep: " << maxStep / kpc  << " kpc, ";
	s << "tolerance: " << tolerance << "\n";

	s << "using diffusion tensor: " << diffusionTensor -> getDescription() << "\n";
	return s.str();
}
