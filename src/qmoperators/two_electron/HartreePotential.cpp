#include "MRCPP/MWOperators"
#include "MRCPP/MWFunctions"
#include "MRCPP/Printer"
#include "MRCPP/Timer"

#include "HartreePotential.h"
#include "chemistry/chemistry_utils.h"
#include "parallel.h"
#include "qmfunctions/Orbital.h"
#include "qmfunctions/density_utils.h"
#include "qmfunctions/orbital_utils.h"
#include "qmfunctions/qmfunction_utils.h"
#include "utils/print_utils.h"
#include "utils/math_utils.h"

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

    Nuclei nucs = this->nuclei;

    Timer timer;
    // Density rho_el(false);
    density::compute(prec, this->rho_el, Phi, DensityType::Total);
    auto period = (*MRA).getWorldBox().getScalingFactor(0);
    auto rc = 0.5;
    auto b_smear = [nucs, rc](const mrcpp::Coord<3> &r) -> double {
        auto g_rc = 0.0;
        for (auto &nuc : nucs) {
            auto R = math_utils::calc_distance(r, nuc.getCoord());
            auto g_i = 0.0;
            if (R <= rc and R >= 0) {
                g_i = 21.0 * std::pow((R - rc), 3.0) * (6.0 * R * R + 3.0 * R * rc + rc * rc) /
                      (5.0 * mrcpp::pi * std::pow(rc, 8.0));
            }
            g_rc += g_i * nuc.getCharge();
        }
        return g_rc;
    };

    Density rho_tree(false);
    if (not rho_tree.hasReal()) rho_tree.alloc(NUMBER::Real);

    mrcpp::project<3>(1.0e-6, rho_tree.real(), b_smear);
    // qmfunction::add(rho_tree, 1.0, this->b_smeared, 1.0, this->b_corr, -1.0);

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

void HartreePotential::setupGlobalPotential(double prec) {
    if (this->poisson == nullptr) MSG_ERROR("Poisson operator not initialized");

    PoissonOperator &P = *this->poisson;
    QMFunction &V = *this;
    QMFunction &rho = this->density;
    Nuclei &nucs = this->nuclei;
    println(0, "nucs.size( ) " << nucs.size())

    mrcpp::FunctionTree<3> V_tree(*MRA);
    mrcpp::FunctionTree<3> corr_tree(*MRA);


    if (V.hasReal()) MSG_ERROR("Potential not properly cleared");
    if (V.hasImag()) MSG_ERROR("Potential not properly cleared");

    // Adjust precision by system size
    double abs_prec = prec / rho.norm();
    bool need_to_apply = not(V.isShared()) or mpi::share_master();

    Timer timer;
    V.alloc(NUMBER::Real);

    static int count = 0;
    count++;
    println(0, "count  " << count)
    println(0, "Hartree version")

    if (need_to_apply) {
        if (this->isFarField()) {
            mrcpp::apply_near_field(abs_prec, V.real(), P, rho.real());
        } else {
            println(0, "here on hartree?")
            mrcpp::apply(abs_prec, V_tree, P, rho.real());
        }
    }
    auto rc = 0.5;
    // if (count == 1) {
    auto V_smear = [nucs, rc](const mrcpp::Coord<3> &r) -> double {
    // auto V_smear = [nucs, rc](const mrcpp::Coord<3> &r) -> double {

        auto v_g = 0.0;
        for (auto &nuc : nucs) {
            auto R = math_utils::calc_distance(r, nuc.getCoord());
            auto v_i = 0.0;
            if (R <= rc and R >= 0) {
                v_i = -(9.0 * std::pow(R, 7.0) - 30.0 * std::pow(R, 6.0) * rc + 28.0 * std::pow(R, 5.0) * rc * rc -
                        14.0 * R * R * std::pow(rc, 5.0) + 12.0 * std::pow(rc, 7.0)) /
                      (5.0 * std::pow(rc, 8.0));
            } else {
                v_i = -1.0 / R;
            }
            v_g += v_i * nuc.getCharge();
        }
        return v_g;
    };
    mrcpp::project<3>(prec, corr_tree, V_smear);
    // }

    println(0, "V_tree " << V_tree.getSquareNorm());
    mrcpp::add(1.0e-6, V.real(), 1.0, V_tree, -1.0, corr_tree);
    mpi::share_function(V, 0, 22445, mpi::comm_share);
    print_utils::qmfunction(2, "Coulomb potential", V, timer);
    println(0, "V.real " << V.real().getSquareNorm());



    auto b_smear = [nucs, rc](const mrcpp::Coord<3> &r) -> double {
        auto g_rc = 0.0;
        for (auto &nuc : nucs) {
            auto R = math_utils::calc_distance(r, nuc.getCoord());
            auto g_i = 0.0;
            if (R <= rc and R >= 0) {
                g_i = 21.0 * std::pow((R - rc), 3.0) * (6.0 * R * R + 3.0 * R * rc + rc * rc) /
                      (5.0 * mrcpp::pi * std::pow(rc, 8.0));
            }
            g_rc += g_i * nuc.getCharge();
        }
        return g_rc;
    };

    mrcpp::FunctionTree<3> b_smear_tree(*MRA);
    mrcpp::FunctionTree<3> V_smear_tree(*MRA);

    // mrcpp::PoissonOperator P(*MRA);

    mrcpp::project<3>(1.0e-6, b_smear_tree, b_smear);

    mrcpp::apply<3>(1.0e-6, V_smear_tree, P, b_smear_tree);
    println(0, "V_smear_tree.norm " << V_smear_tree.getSquareNorm())
    println(0, "corr_tree.norm " << corr_tree.getSquareNorm())


}

} // namespace mrchem
