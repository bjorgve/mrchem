#include "MRCPP/Plotter"
#include "MRCPP/Printer"
#include "MRCPP/Timer"

#include "SmearedNuclearPotential.h"
#include "chemistry/Nucleus.h"
#include "chemistry/chemistry_utils.h"
#include "parallel.h"
#include "qmfunctions/qmfunction_utils.h"
#include "utils/math_utils.h"
#include "utils/periodic_utils.h"
#include "utils/print_utils.h"

using mrcpp::Printer;
using mrcpp::Timer;

namespace mrchem {

/** @brief projects an analytic expression for a smoothed nuclear potential
 *
 * @param[in] nucs collection of nuclei that defines the potential
 * @param[in] prec precision used both in smoothing and projection
 *
 * We create two different analytic functions for the nuclear potential:
 *
 * 1) this->func: The total potential from all nuclei of the system.
 *                This is needed later for analytic calculations.
 *
 * 2) loc_func: Temporary function that contains only some of the
 *              nuclei, which are distributed among the available
 *              MPIs. This is used only for the projection below.
 */
SmearedNuclearPotential::SmearedNuclearPotential(Nuclei nucs,
                                                 double proj_prec,
                                                 double smooth_prec,
                                                 double rc,
                                                 bool mpi_share)
        : NuclearPotential(nucs, proj_prec, smooth_prec, mpi_share) {
    if (proj_prec < 0.0) MSG_ABORT("Negative projection precision");
    if (smooth_prec < 0.0) smooth_prec = proj_prec;

    int pprec = Printer::getPrecision();
    int w0 = Printer::getWidth() - 1;
    int w1 = 5;
    int w2 = 8;
    int w3 = 2 * w0 / 9;
    int w4 = w0 - w1 - w2 - 3 * w3;

    std::stringstream o_head;
    o_head << std::setw(w1) << "N";
    o_head << std::setw(w2) << "Atom";
    o_head << std::string(w4, ' ');
    o_head << std::setw(w3) << "Charge";
    o_head << std::setw(w3) << "Precision";
    o_head << std::setw(w3) << "Smoothing";

    mrcpp::print::header(2, "Projecting nuclear potential");
    println(2, o_head.str());
    mrcpp::print::separator(2, '-');

    auto period = (*MRA).getWorldBox().getScalingFactor(0);
    nucs = periodic::periodify_nuclei(nucs, period * 2.0);

    NuclearFunction loc_func;
    double c = 0.00435 * smooth_prec;
    for (int k = 0; k < nucs.size(); k++) {
        const Nucleus &nuc = nucs[k];
        double Z = nuc.getCharge();
        double Z_5 = std::pow(Z, 5.0);
        double smooth = std::pow(c / Z_5, 1.0 / 3.0);

        // All projection must be done on grand master in order to be exact
        int proj_rank = (mpi::numerically_exact) ? 0 : k % mpi::orb_size;

        this->func.push_back(nuc, smooth);
        if (mpi::orb_rank == proj_rank) loc_func.push_back(nuc, smooth);

        std::stringstream o_row;
        o_row << std::setw(w1) << k;
        o_row << std::setw(w2) << nuc.getElement().getSymbol();
        o_row << std::string(w4, ' ');
        o_row << std::setw(w3) << std::setprecision(pprec) << std::scientific << Z;
        o_row << std::setw(w3) << std::setprecision(pprec) << std::scientific << smooth_prec;
        o_row << std::setw(w3) << std::setprecision(pprec) << std::scientific << smooth;
        println(2, o_row.str());
    }

    auto V_smear = [nucs, rc](const mrcpp::Coord<3> &r) -> double {
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

    Timer t_tot;

    // Scale precision by system size
    int Z_tot = chemistry::get_total_charge(nucs);
    double abs_prec = proj_prec / Z_tot;

    QMFunction V_loc(false);
    QMFunction V_smear_loc(false);

    Timer t_loc;
    qmfunction::project(V_loc, loc_func, NUMBER::Real, abs_prec);
    t_loc.stop();

    qmfunction::project(V_smear_loc, V_smear, NUMBER::Real, abs_prec);

    if (not this->hasReal()) this->alloc(NUMBER::Real);
    qmfunction::add(*this, 1.0, V_loc, -1.0, V_smear_loc, -1);
    t_tot.stop();
    mrcpp::print::separator(2, '-');
    print_utils::qmfunction(2, "Local potential", V_loc, t_loc);
    mrcpp::print::footer(2, t_tot, 2);
}

} // namespace mrchem
