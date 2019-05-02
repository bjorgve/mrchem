#include <Eigen/Eigenvalues>
#include <fstream>

#include "MRCPP/MWFunctions"
#include "MRCPP/MWOperators"
#include "MRCPP/Printer"
#include "MRCPP/Timer"

#include "getkw/Getkw.hpp"
#include "getkw/Keyword.hpp"
#include "getkw/Section.hpp"

#include "initial_guess/core.h"
#include "initial_guess/gto.h"
#include "initial_guess/sad.h"
#include "utils/math_utils.h"

#include "qmfunctions/Density.h"
#include "qmfunctions/Orbital.h"

#include "SCFDriver.h"
#include "chemistry/Molecule.h"

#include "qmfunctions/density_utils.h"
#include "qmfunctions/orbital_utils.h"

#include "scf_solver/EnergyOptimizer.h"
#include "scf_solver/KAIN.h"
#include "scf_solver/LinearResponseSolver.h"
#include "scf_solver/OrbitalOptimizer.h"

#include "qmoperators/one_electron/ElectricFieldOperator.h"
#include "qmoperators/one_electron/KineticOperator.h"
#include "qmoperators/one_electron/NuclearOperator.h"
#include "qmoperators/two_electron/CoulombOperator.h"
#include "qmoperators/two_electron/ExchangeOperator.h"
#include "qmoperators/two_electron/FockOperator.h"
#include "qmoperators/two_electron/XCOperator.h"

#include "properties/DipoleMoment.h"
#include "properties/GeometryDerivatives.h"
#include "properties/Magnetizability.h"
#include "properties/Polarizability.h"

#include "qmoperators/one_electron/H_BB_dia.h"
#include "qmoperators/one_electron/H_B_dip.h"
#include "qmoperators/one_electron/H_E_dip.h"
#include "qmoperators/one_electron/H_M_pso.h"
#include "qmoperators/one_electron/NuclearGradientOperator.h"
#include "qmoperators/one_electron/X_rm3.h"

using mrcpp::Printer;
using mrcpp::Timer;
using namespace std;

namespace mrchem {

extern mrcpp::MultiResolutionAnalysis<3> *MRA; // Global MRA

SCFDriver::SCFDriver(Getkw &input) {
    max_scale = MRA->getMaxScale();
    rel_prec = input.get<double>("rel_prec");
    nuc_prec = input.get<double>("nuc_prec");
    mock_periodic = input.get<bool>("mock_periodic");

    gauge = input.getDblVec("mra.gauge_origin");
    center_of_mass = input.get<bool>("mra.center_of_mass");
    center_of_charge = input.get<bool>("mra.center_of_charge");
    periodic = input.get<bool>("mra.periodic");

    diff_kin = input.get<std::string>("derivatives.kinetic");
    diff_orb = input.get<std::string>("derivatives.h_orb");
    diff_pso = input.get<std::string>("derivatives.h_pso");

    calc_scf_energy = input.get<bool>("properties.scf_energy");
    calc_dipole_moment = input.get<bool>("properties.dipole_moment");
    calc_quadrupole_moment = input.get<bool>("properties.quadrupole_moment");
    calc_magnetizability = input.get<bool>("properties.magnetizability");
    calc_nmr_shielding = input.get<bool>("properties.nmr_shielding");
    calc_hyperfine_coupling = input.get<bool>("properties.hyperfine_coupling");
    calc_spin_spin_coupling = input.get<bool>("properties.spin_spin_coupling");
    calc_polarizability = input.get<bool>("properties.polarizability");
    calc_hyperpolarizability = input.get<bool>("properties.hyperpolarizability");
    calc_optical_rotation = input.get<bool>("properties.optical_rotation");
    calc_geometry_derivatives = input.get<bool>("properties.geometry_derivatives");

    nmr_perturbation = input.get<std::string>("nmrshielding.perturbation");
    nmr_nucleus_k = input.getIntVec("nmrshielding.nucleus_k");
    hfcc_nucleus_k = input.getIntVec("hyperfinecoupling.nucleus_k");
    sscc_nucleus_k = input.getIntVec("spinspincoupling.nucleus_k");
    sscc_nucleus_l = input.getIntVec("spinspincoupling.nucleus_l");
    pol_velocity = input.get<bool>("polarizability.velocity");
    pol_frequency = input.getDblVec("polarizability.frequency");
    optrot_velocity = input.get<bool>("opticalrotation.velocity");
    optrot_frequency = input.getDblVec("opticalrotation.frequency");
    optrot_perturbation = input.get<std::string>("opticalrotation.perturbation");

    mol_charge = input.get<int>("molecule.charge");
    mol_multiplicity = input.get<int>("molecule.multiplicity");
    mol_coords = input.getData("molecule.coords");

    wf_restricted = input.get<bool>("wavefunction.restricted");
    wf_method = input.get<std::string>("wavefunction.method");

    if (wf_method == "dft") {
        dft_spin = input.get<bool>("dft.spin");
        dft_use_gamma = input.get<bool>("dft.use_gamma");
        dft_cutoff = input.get<double>("dft.density_cutoff");
        dft_func_coefs = input.getDblVec("dft.func_coefs");
        dft_func_names = input.getData("dft.functionals");
    }

    scf_run = input.get<bool>("scf.run");
    scf_start = input.get<std::string>("scf.initial_guess");
    scf_kain = input.get<int>("scf.kain");
    scf_max_iter = input.get<int>("scf.max_iter");
    scf_rotation = input.get<int>("scf.rotation");
    scf_canonical = input.get<bool>("scf.canonical");
    scf_write_orbitals = input.get<bool>("scf.write_orbitals");
    scf_orbital_thrs = input.get<double>("scf.orbital_thrs");
    scf_property_thrs = input.get<double>("scf.property_thrs");
    scf_lambda_thrs = input.get<double>("scf.lambda_thrs");
    scf_orbital_prec = input.getDblVec("scf.orbital_prec");

    kin_free_run = input.get<bool>("kineticfree.run");
    kin_free_max_iter = input.get<int>("kineticfree.max_iter");
    kin_free_canonical = input.get<bool>("kineticfree.canonical");
    kin_free_orb_thrs = input.get<double>("kineticfree.orbital_thrs");
    kin_free_prop_thrs = input.get<double>("kineticfree.property_thrs");

    rsp_run = input.get<bool>("response.run");
    rsp_start = input.get<std::string>("response.initial_guess");
    rsp_kain = input.get<int>("response.kain");
    rsp_max_iter = input.get<int>("response.max_iter");
    rsp_canonical = input.get<bool>("response.canonical");
    rsp_write_orbitals = input.get<bool>("response.write_orbitals");
    rsp_orbital_thrs = input.get<double>("response.orbital_thrs");
    rsp_property_thrs = input.get<double>("response.property_thrs");
    rsp_directions = input.getIntVec("response.directions");
    rsp_orbital_prec = input.getDblVec("response.orbital_prec");

    ext_electric = input.get<bool>("externalfield.electric_run");
    ext_magnetic = input.get<bool>("externalfield.magnetic_run");
    if (ext_electric) {
        std::vector<double> tmp = input.getDblVec("externalfield.electric_field");
        ext_electric_field[0] = tmp[0];
        ext_electric_field[1] = tmp[1];
        ext_electric_field[2] = tmp[2];
    }
    if (ext_magnetic) {
        std::vector<double> tmp = input.getDblVec("externalfield.magnetic_field");
        ext_magnetic_field[0] = tmp[0];
        ext_magnetic_field[1] = tmp[1];
        ext_magnetic_field[2] = tmp[2];
    }

    file_start_orbitals = input.get<std::string>("files.start_orbitals");
    file_final_orbitals = input.get<std::string>("files.final_orbitals");
    file_start_x_orbs = input.get<std::string>("files.start_x_orbs");
    file_final_x_orbs = input.get<std::string>("files.final_x_orbs");
    file_start_y_orbs = input.get<std::string>("files.start_y_orbs");
    file_final_y_orbs = input.get<std::string>("files.final_y_orbs");
    file_basis_set = input.get<std::string>("files.basis_set");
    file_dens_mat_a = input.get<std::string>("files.dens_mat_a");
    file_dens_mat_b = input.get<std::string>("files.dens_mat_b");
    file_fock_mat = input.get<std::string>("files.fock_mat");
    file_energy_vec = input.get<std::string>("files.energy_vec");
    file_mo_mat_a = input.get<std::string>("files.mo_mat_a");
    file_mo_mat_b = input.get<std::string>("files.mo_mat_b");

    r_O[0] = 0.0;
    r_O[1] = 0.0;
    r_O[2] = 0.0;
}

bool SCFDriver::sanityCheck() const {
    if (scf_lambda_thrs > 0.0) {
        MSG_ERROR("Recycling of HelmholtzOperators currently disabled");
        return true;
    }
    if (wf_restricted and mol_multiplicity != 1) {
        MSG_ERROR("Restricted open-shell not implemented");
        return false;
    }
    if (calc_quadrupole_moment) {
        MSG_ERROR("Quadrupole moment not implemented");
        return false;
    }
    if (calc_hyperpolarizability) {
        MSG_ERROR("Hyperpolarizability not implemented");
        return false;
    }
    if (calc_optical_rotation) {
        MSG_ERROR("Optical rotation not implemented");
        return false;
    }
    if (calc_nmr_shielding) {
        MSG_ERROR("NMR shielding not implemented");
        return false;
    }
    if (calc_spin_spin_coupling) {
        MSG_ERROR("Spin-spin coupling not implemented");
        return false;
    }
    if (calc_hyperfine_coupling) {
        MSG_ERROR("Hyperfine coupling not implemented");
        return false;
    }
    return true;
}

void SCFDriver::setup() {
    // Setting up molecule
    molecule = std::make_shared<Molecule>(mol_coords, mol_charge, mol_multiplicity);
    molecule->printGeometry();

    // Setting up empty orbitals
    phi = molecule->getOrbitals_p();

    // Defining gauge origin
    auto COM = molecule->calcCenterOfMass();
    auto COC = molecule->calcCenterOfCharge();
    if (center_of_mass) {
        r_O[0] = COM[0];
        r_O[1] = COM[1];
        r_O[2] = COM[2];
    } else if (center_of_charge) {
        r_O[0] = COC[0];
        r_O[1] = COC[1];
        r_O[2] = COC[2];
    } else {
        r_O[0] = gauge[0];
        r_O[1] = gauge[1];
        r_O[2] = gauge[2];
    }

    // Setting up MW operators
    P = std::make_shared<mrcpp::PoissonOperator>(*MRA, rel_prec);
    PH_1 = std::make_shared<mrcpp::PHOperator<3>>(*MRA, 1); // first derivative
    PH_2 = std::make_shared<mrcpp::PHOperator<3>>(*MRA, 2); // second derivative
    ABGV_00 = std::make_shared<mrcpp::ABGVOperator<3>>(*MRA, 0.0, 0.0);
    ABGV_55 = std::make_shared<mrcpp::ABGVOperator<3>>(*MRA, 0.5, 0.5);

    // Setting up perturbation operators
    int nNucs = molecule->getNNuclei();
    h_E = new H_E_dip(r_O);
    h_B = new H_B_dip(useDerivative(diff_orb), r_O);
    h_M = new H_M_pso *[nNucs];
    for (int k = 0; k < nNucs; k++) {
        const mrcpp::Coord<3> &r_K = molecule->getNuclei()[k].getCoord();
        h_M[k] = new H_M_pso(useDerivative(diff_pso), r_K);
    }

    // Setting up properties
    if (nmr_nucleus_k[0] < 0) {
        nmr_nucleus_k.clear();
        for (int k = 0; k < nNucs; k++) nmr_nucleus_k.push_back(k);
    }
    if (sscc_nucleus_k[0] < 0) {
        sscc_nucleus_k.clear();
        for (int k = 0; k < nNucs; k++) sscc_nucleus_k.push_back(k);
    }
    if (sscc_nucleus_l[0] < 0) {
        sscc_nucleus_l.clear();
        for (int l = 0; l < nNucs; l++) sscc_nucleus_l.push_back(l);
    }
    if (hfcc_nucleus_k[0] < 0) {
        hfcc_nucleus_k.clear();
        for (int k = 0; k < nNucs; k++) hfcc_nucleus_k.push_back(k);
    }

    if (calc_magnetizability) {
        for (int d = 0; d < 3; d++) {
            if (rsp_directions[d] == 0) continue;
            rsp_calculations.push_back(h_B, 0.0, true, d, "H_B");
        }
    }
    if (calc_nmr_shielding) {
        for (auto &K : nmr_nucleus_k) {
            if (nmr_perturbation == "B") {
                for (int d = 0; d < 3; d++) {
                    if (rsp_directions[d] == 0) continue;
                    rsp_calculations.push_back(h_B, 0.0, true, d, "H_B");
                }
            } else {
                for (int d = 0; d < 3; d++) {
                    if (rsp_directions[d] == 0) continue;
                    std::string name = "H_M" + std::to_string(K);
                    rsp_calculations.push_back(h_M[K], 0.0, true, d, name);
                }
            }
        }
    }
    if (calc_polarizability) {
        for (auto &omega : pol_frequency) {
            for (int d = 0; d < 3; d++) {
                if (rsp_directions[d] == 0) continue;
                rsp_calculations.push_back(h_E, omega, false, d, "H_E");
            }
        }
    }

    // Setting up Fock operator
    auto &nuclei = molecule->getNuclei();
    T = std::make_shared<KineticOperator>((useDerivative(diff_kin)));
    // All cases need kinetic energy
    fock = std::make_shared<FockOperator>(T);

    // For non-periodic Hartree, HF and DFT we need the coulomb and nuclear part
    if ((wf_method == "hartree" or wf_method == "hf" or wf_method == "dft") and not mock_periodic) {
        J = std::make_shared<CoulombOperator>(P, phi);
        V = std::make_shared<NuclearOperator>(nuclei, nuc_prec);
        fock->getCoulombOperator() = J;
        fock->getNuclearOperator() = V;
    }
    // For periodic Hartree, HF and DFT we need the hartree potential part
    if ((wf_method == "hartree" or wf_method == "hf" or wf_method == "dft") and mock_periodic) {
        J = std::make_shared<CoulombOperator>(P, phi, nuclei);
        fock->getCoulombOperator() = J;
    }
    // For HF we need the full HF exchange
    if (wf_method == "hf") {
        K = std::make_shared<ExchangeOperator>(P, phi);
        fock->getExchangeOperator() = K;
    }
    // For DFT we need the XC operator
    double exx = 1.0;
    if (wf_method == "dft") {
        xcfun = setupFunctional(MRDFT::Gradient);
        XC = std::make_shared<XCOperator>(xcfun, phi);
        fock->getXCOperator() = XC;

        // For hybrid DFT we need a partial HF exchange
        if (xcfun->isHybrid()) {
            exx = xcfun->amountEXX();
            K = std::make_shared<ExchangeOperator>(P, phi);
            fock->getExchangeOperator() = K;
        }
    }
    // HACK we need a better way to decide whether to initialize the external potential operator
    if (ext_electric) Vext = std::make_shared<ElectricFieldOperator>(ext_electric_field, r_O);
    fock->getExtOperator() = Vext;
    fock->build(exx);
}

/** @brief choose the right derivative operator to use
 *
 * param[in] nam of the chosen derivative operator
 *
 * Returns the pointer to the correct derivative operator
 *
 */
std::shared_ptr<mrcpp::DerivativeOperator<3>> SCFDriver::useDerivative(string derivative_name) {
    if (derivative_name == "ph_1") return PH_1;
    if (derivative_name == "abgv_00") return ABGV_00;
    if (derivative_name == "abgv_55") return ABGV_55;
    MSG_FATAL("No such derivative operator");
}

void SCFDriver::clear() {
    for (int k = 0; k < molecule->getNNuclei(); k++) {
        if (h_M[k] != nullptr) delete h_M[k];
    }
    if (h_M != nullptr) delete[] h_M;
    if (h_B != nullptr) delete h_B;
    if (h_E != nullptr) delete h_E;
}

void SCFDriver::setupInitialGroundState() {
    double prec = scf_orbital_prec[0];
    if (scf_start == "gto")
        if (wf_restricted)
            *phi = initial_guess::gto::setup(prec, *molecule, file_basis_set, file_mo_mat_a);
        else
            *phi = initial_guess::gto::setup(prec, *molecule, file_basis_set, file_mo_mat_a, file_mo_mat_b);
    else if (scf_start == "core_sz")
        *phi = initial_guess::core::setup(prec, *molecule, wf_restricted, 1);
    else if (scf_start == "core_dz")
        *phi = initial_guess::core::setup(prec, *molecule, wf_restricted, 2);
    else if (scf_start == "core_tz")
        *phi = initial_guess::core::setup(prec, *molecule, wf_restricted, 3);
    else if (scf_start == "core_qz")
        *phi = initial_guess::core::setup(prec, *molecule, wf_restricted, 4);
    else if (scf_start == "sad_sz")
        *phi = initial_guess::sad::setup(prec, *molecule, wf_restricted, 1);
    else if (scf_start == "sad_dz")
        *phi = initial_guess::sad::setup(prec, *molecule, wf_restricted, 2);
    else if (scf_start == "sad_tz")
        *phi = initial_guess::sad::setup(prec, *molecule, wf_restricted, 3);
    else if (scf_start == "sad_qz")
        *phi = initial_guess::sad::setup(prec, *molecule, wf_restricted, 4);
    else if (scf_start == "mw")
        *phi = orbital::load_orbitals(file_start_orbitals);
    else
        MSG_FATAL("Invalid initial guess");
    orbital::print(*phi);
}

void SCFDriver::setupPerturbedOrbitals(const ResponseCalculation &rsp_calc) {
    molecule->initPerturbedOrbitals(rsp_calc.isDynamic());
    phi_x = molecule->getOrbitalsX_p();
    phi_y = molecule->getOrbitalsY_p();
    *phi_x = orbital::param_copy(*phi);
    *phi_y = orbital::param_copy(*phi);
    if (rsp_start == "MW") *phi_x = orbital::load_orbitals(file_start_x_orbs, rsp_calc.getFileSuffix());
    if (rsp_calc.isDynamic()) {
        if (rsp_start == "MW") *phi_y = orbital::load_orbitals(file_start_y_orbs, rsp_calc.getFileSuffix());
    }
}

void SCFDriver::setupPerturbedOperators(const ResponseCalculation &rsp_calc) {
    if (phi == nullptr) MSG_ERROR("Orbitals not initialized");
    if (phi_x == nullptr) MSG_ERROR("X orbitals not initialized");
    if (phi_y == nullptr) MSG_ERROR("Y orbitals not initialized");

    double xFac = 0.0;
    if (wf_method == "hf") {
        xFac = 1.0;
    } else if (wf_method == "dft") {
        xFac = xcfun->amountEXX();
        xcfun = setupFunctional(MRDFT::Hessian);
        dXC = std::make_shared<XCOperator>(xcfun, phi, phi_x, phi_y);
    }
    if (xFac > mrcpp::MachineZero) NOT_IMPLEMENTED_ABORT;

    dJ = std::make_shared<CoulombOperator>(P, phi, phi_x, phi_y);

    int d = rsp_calc.dir;
    RankOneTensorOperator<3> &dH = *rsp_calc.pert;
    if (rsp_calc.isDynamic()) NOT_IMPLEMENTED_ABORT;

    d_fock = std::make_shared<FockOperator>(nullptr, nullptr, dJ, dK, dXC);
    d_fock->perturbation() += dH[d];
    d_fock->build(xFac);
}

void SCFDriver::run() {
    if (not sanityCheck()) return;

    bool converged = runGroundState();
    if (converged) {
        for (int i = 0; i < rsp_calculations.size(); i++) {
            const ResponseCalculation &rsp_calc = rsp_calculations[i];
            runLinearResponse(rsp_calc);
        }
    }

    molecule->printGeometry();
    molecule->printProperties();
}

bool SCFDriver::runGroundState() {
    if (phi == nullptr) MSG_ERROR("Orbitals not initialized");
    if (fock == nullptr) MSG_ERROR("Fock operator not initialized");
    bool converged = true;

    // Setup initial guess
    setupInitialGroundState();

    auto &F = molecule->getFockMatrix();

    // Optimize orbitals
    if (scf_run) {
        OrbitalOptimizer solver;
        solver.setHistory(scf_kain);
        solver.setMaxIterations(scf_max_iter);
        solver.setRotation(scf_rotation);
        solver.setCanonical(scf_canonical);
        solver.setThreshold(scf_orbital_thrs, scf_property_thrs);
        solver.setOrbitalPrec(scf_orbital_prec[0], scf_orbital_prec[1]);
        converged = solver.optimize(*fock, *phi, F);
    } else {
        fock->setup(rel_prec);
        F = (*fock)(*phi, *phi);
        fock->clear();

        if (scf_canonical) {
            orbital::diagonalize(rel_prec, *phi, F);
        } else {
            orbital::localize(rel_prec, *phi, F);
        }
    }

    // Optimize energy
    if (kin_free_run) {
        EnergyOptimizer solver;
        solver.setMaxIterations(kin_free_max_iter);
        solver.setCanonical(kin_free_canonical);
        solver.setThreshold(kin_free_orb_thrs, kin_free_prop_thrs);
        solver.setOrbitalPrec(rel_prec, rel_prec);
        converged = solver.optimize(*fock, *phi, F);
    }

    if (scf_write_orbitals) orbital::save_orbitals(*phi, file_final_orbitals);

    // Compute requested properties
    if (converged) calcGroundStateProperties();

    return converged;
}

void SCFDriver::runLinearResponse(const ResponseCalculation &rsp_calc) {
    bool dynamic = rsp_calc.isDynamic();
    setupPerturbedOrbitals(rsp_calc);
    setupPerturbedOperators(rsp_calc);

    auto &F = molecule->getFockMatrix();

    if (d_fock->getXCOperator() != nullptr) {
        d_fock->getXCOperator()->setupDensity(rel_prec);
        d_fock->getXCOperator()->setupPotential(rel_prec);
    }

    bool converged = true;
    if (rsp_run) {
        LinearResponseSolver solver(dynamic, *fock, *phi, F);
        solver.setHistory(rsp_kain);
        solver.setMaxIterations(rsp_max_iter);
        solver.setThreshold(rsp_orbital_thrs, rsp_property_thrs);
        solver.setOrbitalPrec(rsp_orbital_prec[0], rsp_orbital_prec[1]);

        fock->setup(rsp_orbital_prec[1]);
        converged = solver.optimize(rsp_calc.freq, *d_fock, *phi_x, *phi_y);
        fock->clear();
    }
    if (rsp_write_orbitals) {
        orbital::save_orbitals(*phi_x, file_final_x_orbs, rsp_calc.getFileSuffix());
        if (dynamic) orbital::save_orbitals(*phi_y, file_final_y_orbs, rsp_calc.getFileSuffix());
    }

    // Compute requested properties
    if (converged) calcLinearResponseProperties(rsp_calc);
}

void SCFDriver::calcGroundStateProperties() {
    auto &F = molecule->getFockMatrix();
    auto &nuclei = molecule->getNuclei();

    if (calc_scf_energy) {
        fock->setup(rel_prec);
        Printer::printHeader(0, "Calculating SCF energy");
        Timer timer;
        SCFEnergy &energy = molecule->getSCFEnergy();
        energy = fock->trace(*phi, F);
        timer.stop();
        Printer::printFooter(0, timer, 2);
        fock->clear();
    }
    if (calc_dipole_moment) {
        Printer::printHeader(0, "Calculating dipole moment");
        Timer timer;
        DoubleVector &nuc = molecule->getDipoleMoment().getNuclear();
        DoubleVector &el = molecule->getDipoleMoment().getElectronic();
        H_E_dip mu(r_O);
        mu.setup(rel_prec);
        nuc = mu.trace(nuclei).real();
        el = mu.trace(*phi).real();
        mu.clear();
        timer.stop();
        Printer::printFooter(0, timer, 2);
    }
    if (calc_geometry_derivatives) {
        Printer::printHeader(0, "Calculating geometry derivatives");
        Timer timer;
        DoubleMatrix &nuc = molecule->getGeometryDerivatives().getNuclear();
        DoubleMatrix &el = molecule->getGeometryDerivatives().getElectronic();
        DoubleVector vecsum = DoubleVector::Zero(3);
        DoubleVector torque = DoubleVector::Zero(3);

        for (int k = 0; k < nuclei.size(); k++) {
            const Nucleus &nuc_k = nuclei[k];
            double Z_k = nuc_k.getCharge();
            const mrcpp::Coord<3> &R_k = nuc_k.getCoord();
            Nuclei nucs;
            nucs.push_back("H", R_k);
            NuclearGradientOperator r_rm3(nuc_k, 1.0e-2);
            r_rm3.setup(1.0e-4);
            for (int l = 0; l < nuclei.size(); l++) {
                if (l == k) continue;
                const Nucleus &nuc_l = nuclei[l];
                double Z_l = nuc_l.getCharge();
                const mrcpp::Coord<3> &R_l = nuc_l.getCoord();
                double r_x = (R_k[0] - R_l[0]);
                double r_y = (R_k[1] - R_l[1]);
                double r_z = (R_k[2] - R_l[2]);
                double R_kl = std::pow(math_utils::calc_distance(R_k, R_l), 3.0);
                ;
                nuc(k, 0) -= Z_k * Z_l * r_x / R_kl;
                nuc(k, 1) -= Z_k * Z_l * r_y / R_kl;
                nuc(k, 2) -= Z_k * Z_l * r_z / R_kl;
            }
            el.row(k) = r_rm3.trace(*phi).real();
            r_rm3.clear();
            vecsum += el.row(k);
            vecsum += nuc.row(k);
            torque[0] += R_k[1] * (el(k, 2) + nuc(k, 2)) - R_k[2] * (el(k, 1) + nuc(k, 1));
            torque[1] += R_k[2] * (el(k, 0) + nuc(k, 0)) - R_k[0] * (el(k, 2) + nuc(k, 2));
            torque[2] += R_k[0] * (el(k, 1) + nuc(k, 1)) - R_k[1] * (el(k, 0) + nuc(k, 0));
        }
        println(0, "nuclear part    ");
        println(0, nuc);
        println(0, "electronic part ");
        println(0, el);
        println(0, "Total force acting on nuclei");
        println(0, vecsum.transpose());
        println(0, "Torque acting on nuclei");
        println(0, torque.transpose());
        timer.stop();
        Printer::printFooter(0, timer, 2);
    }
    if (calc_magnetizability) {
        Printer::printHeader(0, "Calculating diamagnetic magnetizability");
        Timer timer;
        DoubleMatrix &dia = molecule->getMagnetizability().getDiamagnetic();
        H_BB_dia h(r_O);
        h.setup(rel_prec);
        dia = -h.trace(*phi).real();
        h.clear();
        timer.stop();
        Printer::printFooter(0, timer, 2);
    }
    /*
    if (calc_nmr_shielding) {
        Printer::printHeader(0, "Calculating diamagnetic NMR shielding");
        Timer timer;
        for (int k = 0; k < nmr_nucleus_k.size(); k++) {
            int K = nmr_nucleus_k[k];
            NMRShielding &nmr = molecule->getNMRShielding(K);
            MatrixXd &dia = nmr.getDiamagnetic();
            const mrcpp::Coord<3> &r_K = nmr.getNucleus().getCoord();
            H_BM_dia h(r_O, r_K);
            h.setup(rel_prec);
            dia = h.trace(*phi);
            h.clear();
        }
        timer.stop();
        Printer::printFooter(0, timer, 2);
    }
    if (calc_spin_spin_coupling) {
        Printer::printHeader(0, "Calculating diamagnetic spin-spin coupling");
        Timer timer;
        for (int k = 0; k < sscc_nucleus_k.size(); k++) {
            int K = sscc_nucleus_k[k];
            for (int l = 0; l < sscc_nucleus_l.size(); l++) {
                int L = sscc_nucleus_l[l];
                if (K == L) continue;
                SpinSpinCoupling &sscc = molecule->getSpinSpinCoupling(K, L);
                MatrixXd &dia = sscc.getDiamagnetic();
                const mrcpp::Coord<3> &r_K = sscc.getNucleusK().getCoord();
                const mrcpp::Coord<3> &r_L = sscc.getNucleusL().getCoord();
                NOT_IMPLEMENTED_ABORT;
            }
        }
        timer.stop();
        Printer::printFooter(0, timer, 2);
    }
    if (calc_hyperfine_coupling) {
        Printer::printHeader(0, "Calculating HyperFine Coupling Constant");
        Timer timer;

        for (int k = 0; k < hfcc_nucleus_k.size(); k++) {
            int K = hfcc_nucleus_k[k];
            HyperFineCoupling &hfc = molecule->getHyperFineCoupling(K);
            const Nuclei &nucs = molecule->getNuclei();
            const Nucleus &nuc = nucs[K];
            const mrcpp::Coord<3> &r_K = nuc.getCoord();
            NOT_IMPLEMENTED_ABORT;
        }
        timer.stop();
        Printer::printFooter(0, timer, 2);
    }
    */
}

void SCFDriver::calcLinearResponseProperties(const ResponseCalculation &rsp_calc) {
    int j = rsp_calc.dir;
    double freq = rsp_calc.freq;

    if (calc_magnetizability and rsp_calc.pert == h_B) {
        Printer::printHeader(0, "Calculating paramagnetic magnetizability");
        Timer timer;
        DoubleMatrix &para = molecule->getMagnetizability().getParamagnetic();
        h_B->setup(rel_prec);
        para.row(j) = -h_B->trace(*phi, *phi_x, *phi_y).real();
        h_B->clear();
        timer.stop();
        Printer::printFooter(0, timer, 2);
    }

    if (calc_polarizability and rsp_calc.pert == h_E) {
        Printer::printHeader(0, "Calculating polarizability");
        Timer timer;
        DoubleMatrix &pol = molecule->getPolarizability(freq).getTensor();
        h_E->setup(rel_prec);
        pol.row(j) = -h_E->trace(*phi, *phi_x, *phi_y).real();
        h_E->clear();
        timer.stop();
        Printer::printFooter(0, timer, 2);
    }
    /*
    if (calc_nmr_shielding) {
        if (nmr_perturbation == "B" and rsp_calc.pert == h_B) {
            Timer timer;
            Printer::printHeader(0, "Calculating paramagnetic NMR shielding ");
            for (int k = 0; k < nmr_nucleus_k.size(); k++) {
                int K = nmr_nucleus_k[k];
                MatrixXd &para = molecule->getNMRShielding(K).getParamagnetic();
                h_M[K]->setup(rel_prec);
                para.row(j) = -h_M[K]->trace(*phi, *phi_x, *phi_y);
                h_M[K]->clear();
            }
            timer.stop();
            Printer::printFooter(0, timer, 2);
        }
        if (nmr_perturbation == "M") {
            for (int k = 0; k < nmr_nucleus_k.size(); k++) {
                int K = nmr_nucleus_k[k];
                if (rsp_calc.pert == h_M[K]) {
                    Timer timer;
                    Printer::printHeader(0, "Calculating paramagnetic NMR shielding");

                    MatrixXd &para = molecule->getNMRShielding(K).getParamagnetic();
                    h_B->setup(rel_prec);
                    para.col(j) = -h_B->trace(*phi, *phi_x, *phi_y);
                    h_B->clear();

                    timer.stop();
                    Printer::printFooter(0, timer, 2);
                }
            }
        }
    }
    if (calc_spin_spin_coupling) {
        Printer::printHeader(0, "Calculating paramagnetic spin-spin coupling");
        Timer timer;
        NOT_IMPLEMENTED_ABORT;
        timer.stop();
        Printer::printFooter(0, timer, 2);
    }
    */
}

/** @brief Build initial density grid for DFT
 *
 * This will refine the density grid in XCFunctional around each nuclear site
 * of the molecule. Uses the refinement algorithm for Gaussians in MRCPP by
 * placing a narrow Gaussian function on each atom with exponent set to the
 * square of the nuclear charge. This grid will be adaptively refined during
 * the SCF procedure.
 */
void SCFDriver::setupInitialGrid(mrdft::XCFunctional &func, const Molecule &mol) {
    Printer::printHeader(0, "Initialize DFT grid");
    println(0, " Nr  Element                        nNodes          nPoints");
    Printer::printSeparator(0, '-');
    Timer timer;
    const Nuclei &nucs = mol.getNuclei();
    for (int k = 0; k < nucs.size(); k++) {
        func.buildGrid(nucs[k].getCharge(), nucs[k].getCoord());
        printout(0, std::setw(3) << k);
        printout(0, std::setw(7) << nucs[k].getElement().getSymbol());
        printout(0, std::setw(32) << func.getNNodes());
        printout(0, std::setw(17) << func.getNPoints() << "\n");
    }
    timer.stop();
    Printer::printFooter(0, timer, 2);
}

/** @brief helper routine to set up the correct parameters in the functional before using it
 *
 * param[in] order the requested order of the derivative (order=1 for SCF)
 *
 */

std::shared_ptr<mrdft::XCFunctional> SCFDriver::setupFunctional(int order) {
    auto fun = std::make_shared<mrdft::XCFunctional>(*MRA, dft_spin);
    for (int i = 0; i < dft_func_names.size(); i++) {
        double f_coef = dft_func_coefs[i];
        const std::string &f_name = dft_func_names[i];
        fun->setFunctional(f_name, f_coef);
    }
    fun->setUseGamma(dft_use_gamma);
    fun->setDensityCutoff(dft_cutoff);
    fun->evalSetup(order);
    // setupInitialGrid(*fun, *molecule);
    return fun;
}

} // namespace mrchem
