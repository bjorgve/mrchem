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

#include "qmoperators/one_electron/PositionOperator.h"
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
                                   const double &rc)

        : CoulombPotential(P, Phi)
        , nuclei(nucs)
        , rc(rc) {}

void HartreePotential::setupGlobalDensity(double prec) {
    if (hasDensity()) return;
    if (this->orbitals == nullptr) MSG_ERROR("Orbitals not initialized");

    OrbitalVector &Phi = *this->orbitals;
    Density &rho = this->density;

    Timer timer;
    // Density rho_el(false);
    density::compute(prec, this->rho_el, Phi, DensityType::Total);
    auto period = (*MRA).getWorldBox().getScalingFactor(0);
    this->b_smeared =
        chemistry::compute_nuclear_density_smeared(this->apply_prec, this->nuclei, this->rc, period * 2.0);
    this->b_smeared.real().rescale(0.0);
        //chemistry::hack_density(this->apply_prec, this->nuclei, this->rc, period * 2.0, Phi);
    this->b_corr = chemistry::calc_bcorr(this->apply_prec, this->nuclei, this->rc, period * 2.0, Phi);
    // this->b_corr.real().rescale(0.0);

    Density rho_tree(false);
    if (not rho_tree.hasReal()) rho_tree.alloc(NUMBER::Real);
    qmfunction::add(rho_tree, 1.0, this->b_smeared, 1.0, this->b_corr, -1.0);

    qmfunction::add(rho, 1.0, this->rho_el, 1.0, rho_tree, -1.0);
    println(0, "rho_el.int " << this->rho_el.integrate())
    println(0, "rho_tree.int " << rho_tree.integrate())
    // if (std::norm(rho.integrate()) > this->apply_prec) MSG_ABORT("Non-zero net charge on unit cell");
    auto dip_oper = PositionOperator({0.0, 0.0, 0.0});
    dip_oper.setup(prec);
    auto mu = dip_oper.trace(rho);
    println(0, "mu look " << mu[0] << " " << mu[1] << " " << mu[2])
    dip_oper.clear();
    timer.stop();
    // double t = timer.getWallTime();
    // int n = rho.getNNodes(NUMBER::Total);
    // Printer::printTree(0, "Hartree density", n, t);
}

} // namespace mrchem
