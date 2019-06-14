#include "MRCPP/MWOperators"
#include "MRCPP/Printer"
#include "MRCPP/Timer"

#include "HartreePotential.h"
#include "chemistry/chemistry_utils.h"
#include "parallel.h"
#include "qmfunctions/Orbital.h"
#include "qmfunctions/density_utils.h"
#include "qmfunctions/orbital_utils.h"
#include "qmfunctions/qmfunction_utils.h"

using mrcpp::PoissonOperator;
using mrcpp::Printer;
using mrcpp::Timer;

namespace mrchem {

/** @brief constructor
 *
 * @param P: MW Poisson operator
 * @param Phi: orbitals defining the operator
 *
 * Density is not spin-separated since this operator requires only total density.
 * This operator will always point to the same OrbitalVector, but the orbitals within
 * the vector can change throughout the calculation. The density and (*this)
 * QMPotential is uninitialized at this point and will be computed at setup.
 */

HartreePotential::HartreePotential(std::shared_ptr<mrcpp::PoissonOperator> P,
                                   std::shared_ptr<OrbitalVector> Phi,
                                   const Nuclei &nucs,
                                   double prec)
        : CoulombPotential(P, Phi) {
    this->nuclei = nucs;
    this->local = false;
    this->nuc_prec = prec;
}

void HartreePotential::setupGlobalDensity(double prec) {
    if (hasDensity()) return;
    if (this->orbitals == nullptr) MSG_ERROR("Orbitals not initialized");

    OrbitalVector &Phi = *this->orbitals;
    Density &rho = this->density;

    Timer timer;
    Density rho_el(false);
    density::compute(prec, rho_el, Phi, DENSITY::Total);
    auto rho_nuc = chemistry::compute_nuclear_density(this->apply_prec, this->nuclei, 1.0 / this->nuc_prec);
    qmfunction::add(rho, 1.0, rho_el, -1.0, rho_nuc, -1.0);
    if (std::norm(rho.integrate()) > this->apply_prec) MSG_FATAL("Non-zero net charge on unit cell");

    timer.stop();
    double t = timer.getWallTime();
    int n = rho.getNNodes(NUMBER::Total);
    Printer::printTree(0, "Hartree density", n, t);
}

} // namespace mrchem
