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

#include "qmoperators/one_electron/H_E_dip.h"

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


    //// BCOR PART
    auto dip_oper = H_E_dip({0.0, 0.0, 0.0});
    dip_oper.setup(prec);
    DoubleVector nuc_dip = -dip_oper.trace(nucs).real();
    DoubleVector dip_el = dip_oper.trace(Phi).real();
    DoubleVector tot_dip = nuc_dip + dip_el;

    auto new_charge = tot_dip[2]/8.0;

    mrcpp::Coord<3> pos_coord{0.0, 0.0, 4.};
    mrcpp::Coord<3> neg_coord{0.0, 0.0, -4.};
    std::vector<double> charges{-new_charge, new_charge};
    std::vector<mrcpp::Coord<3>> coords{pos_coord, neg_coord};

    dip_oper.clear();

    auto b_corr_smear = [charges, coords, rc](const mrcpp::Coord<3> &r) -> double {
        auto g_rc = 0.0;
        for (auto i = 0; i < charges.size(); i++) {
            auto R = math_utils::calc_distance(r, coords[i]);
            auto g_i = 0.0;
            if (R <= rc and R >= 0) {
                g_i = 21.0 * std::pow((R - rc), 3.0) * (6.0 * R * R + 3.0 * R * rc + rc * rc) /
                      (5.0 * mrcpp::pi * std::pow(rc, 8.0));
            }
            g_rc += g_i * charges[i];
        }
        return g_rc;
    };

    Density rho_corr_tree(false);
    if (not rho_corr_tree.hasReal()) rho_corr_tree.alloc(NUMBER::Real);
    mrcpp::project<3>(1.0e-6, rho_corr_tree.real(), b_corr_smear);
    //// BCOR PART

    Density rho_tree(false);
    if (not rho_tree.hasReal()) rho_tree.alloc(NUMBER::Real);

    mrcpp::project<3>(1.0e-6, rho_tree.real(), b_smear);
    // qmfunction::add(rho_tree, 1.0, this->b_smeared, 1.0, this->b_corr, -1.0);
    Density add_tree(false);
    if (not add_tree.hasReal()) add_tree.alloc(NUMBER::Real);
    qmfunction::add(add_tree, 1.0, rho_tree, 0.0, rho_corr_tree, -1.0);

    qmfunction::add(rho, 1.0, this->rho_el, 1.0, add_tree, -1.0);
    println(0, "rho_el.int " << this->rho_el.integrate())
    println(0, "rho_tree.int " << rho_tree.integrate())
    // if (std::norm(rho.integrate()) > this->apply_prec) MSG_ABORT("Non-zero net charge on unit cell");
    auto dip_oper_2 = PositionOperator({0.0, 0.0, 0.0});
    dip_oper_2.setup(prec);
    auto mu = dip_oper_2.trace(rho);
    println(0, "mu look " << mu[0] << " " << mu[1] << " " << mu[2])
    dip_oper_2.clear();
    timer.stop();
    // double t = timer.getWallTime();
    // int n = rho.getNNodes(NUMBER::Total);
    // Printer::printTree(0, "Hartree density", n, t);
}

void HartreePotential::setupGlobalPotential(double prec) {
    if (this->poisson == nullptr) MSG_ERROR("Poisson operator not initialized");

    PoissonOperator &P = *this->poisson;
    OrbitalVector &Phi = *this->orbitals;
    QMFunction &V = *this;
    QMFunction &rho = this->density;
    Nuclei &nucs = this->nuclei;
    println(0, "nucs.size( ) " << nucs.size())

    mrcpp::FunctionTree<3> V_tree(*MRA);
    mrcpp::FunctionTree<3> V_smear_tree(*MRA);
    mrcpp::FunctionTree<3> V_corr_tree(*MRA);
    mrcpp::FunctionTree<3> add_tree(*MRA);


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
    mrcpp::project<3>(prec, V_smear_tree, V_smear);
    // }

    //// BCOR PART
    auto dip_oper = H_E_dip({0.0, 0.0, 0.0});
    dip_oper.setup(prec);
    DoubleVector nuc_dip = -dip_oper.trace(nucs).real();
    DoubleVector dip_el = dip_oper.trace(Phi).real();
    DoubleVector tot_dip = nuc_dip + dip_el;

    auto new_charge = tot_dip[2]/8.0;

    mrcpp::Coord<3> pos_coord{0.0, 0.0, 4.};
    mrcpp::Coord<3> neg_coord{0.0, 0.0, -4.};
    std::vector<double> charges{-new_charge, new_charge};
    std::vector<mrcpp::Coord<3>> coords{pos_coord, neg_coord};
    std::vector<mrcpp::Coord<3>> pos_coords{};
    std::vector<mrcpp::Coord<3>> neg_coords{};
    pos_coords.push_back(pos_coord);
    neg_coords.push_back(neg_coord);

    // auto period = (*MRA).getWorldBox().getScalingFactor(0);
    // for (auto x  = -5; x < 6; x ++) {
    //     for (auto y = -5; y < 6; y++) {
    //         for (auto z = -5; z < 6; z++) {
    //             mrcpp::Coord<3> tmp_pos{pos_coord[0] + x*period*2.0, pos_coord[1] + y*period*2.0, pos_coord[2] + z*period*2.0 };
    //             mrcpp::Coord<3> tmp_neg{neg_coord[0] + x*period*2.0, neg_coord[1] + y*period*2.0, neg_coord[2] + z*period*2.0 };
    //             pos_coords.push_back(tmp_pos);
    //             neg_coords.push_back(tmp_neg);
    //         }
    //     }
    // }

    dip_oper.clear();

    auto V_corr_smear = [charges, neg_coords, pos_coords, rc](const mrcpp::Coord<3> &r) -> double {
        auto v_g = 0.0;
        for (auto i = 0; i < pos_coords.size(); i++) {
            auto R = math_utils::calc_distance(r, pos_coords[i]);
            auto v_i = 0.0;
            if (R <= rc and R >= 0) {
                v_i = -(9.0 * std::pow(R, 7.0) - 30.0 * std::pow(R, 6.0) * rc + 28.0 * std::pow(R, 5.0) * rc * rc -
                        14.0 * R * R * std::pow(rc, 5.0) + 12.0 * std::pow(rc, 7.0)) /
                      (5.0 * std::pow(rc, 8.0));
            } else {
                v_i = -1.0 / R;
            }
            v_g += v_i * charges[1];
        }

        for (auto i = 0; i < neg_coords.size(); i++) {
            auto R = math_utils::calc_distance(r, neg_coords[i]);
            auto v_i = 0.0;
            if (R <= rc and R >= 0) {
                v_i = -(9.0 * std::pow(R, 7.0) - 30.0 * std::pow(R, 6.0) * rc + 28.0 * std::pow(R, 5.0) * rc * rc -
                        14.0 * R * R * std::pow(rc, 5.0) + 12.0 * std::pow(rc, 7.0)) /
                      (5.0 * std::pow(rc, 8.0));
            } else {
                v_i = -1.0 / R;
            }
            v_g += v_i * charges[0];
        }
        return v_g;
    };

    mrcpp::project<3>(1.0e-6, V_corr_tree, V_corr_smear);
    //// BCOR PART
    println(0, "V_tree " << V_tree.getSquareNorm());
    mrcpp::add(1.0e-6, add_tree, 1.0, V_smear_tree, 0.0, V_corr_tree);
    mrcpp::add(1.0e-6, V.real(), 1.0, V_tree, -1.0, add_tree);
    mpi::share_function(V, 0, 22445, mpi::comm_share);
    print_utils::qmfunction(2, "Coulomb potential", V, timer);
    println(0, "V.real " << V.real().getSquareNorm());



    // auto b_smear = [nucs, rc](const mrcpp::Coord<3> &r) -> double {
    //     auto g_rc = 0.0;
    //     for (auto &nuc : nucs) {
    //         auto R = math_utils::calc_distance(r, nuc.getCoord());
    //         auto g_i = 0.0;
    //         if (R <= rc and R >= 0) {
    //             g_i = 21.0 * std::pow((R - rc), 3.0) * (6.0 * R * R + 3.0 * R * rc + rc * rc) /
    //                   (5.0 * mrcpp::pi * std::pow(rc, 8.0));
    //         }
    //         g_rc += g_i * nuc.getCharge();
    //     }
    //     return g_rc;
    // };
    //
    // mrcpp::FunctionTree<3> b_smear_tree(*MRA);
    //
    // // mrcpp::PoissonOperator P(*MRA);
    //
    // mrcpp::project<3>(1.0e-6, b_smear_tree, b_smear);
    //
    // mrcpp::apply<3>(1.0e-6, V_smear_tree, P, b_smear_tree);
    // println(0, "V_smear_tree.norm " << V_smear_tree.getSquareNorm())
    // println(0, "corr_tree.norm " << corr_tree.getSquareNorm())


}

} // namespace mrchem
