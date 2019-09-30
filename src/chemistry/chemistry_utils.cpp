/*
 * MRChem, a numerical real-space code for molecular electronic structure
 * calculations within the self-consistent field (SCF) approximations of quantum
 * chemistry (Hartree-Fock and Density Functional Theory).
 * Copyright (C) 2019 Stig Rune Jensen, Luca Frediani, Peter Wind and contributors.
 *
 * This file is part of MRChem.
 *
 * MRChem is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * MRChem is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with MRChem.  If not, see <https://www.gnu.org/licenses/>.
 *
 * For information on the complete list of contributors to MRChem, see:
 * <https://mrchem.readthedocs.io/>
 */

#include "MRCPP/Gaussians"
#include "MRCPP/Printer"
#include "MRCPP/Plotter"

#include "Nucleus.h"
#include "chemistry_utils.h"
#include "qmfunctions/Density.h"
#include "qmfunctions/density_utils.h"
#include "utils/math_utils.h"

namespace mrchem {

/** @brief computes the repulsion self energy of a set of nuclei
 *
 * @param[in] nucs the set of nuclei
 *
 */
double chemistry::compute_nuclear_repulsion(const Nuclei &nucs) {
    int nNucs = nucs.size();
    double E_nuc = 0.0;
    for (int i = 0; i < nNucs; i++) {
        const Nucleus &nuc_i = nucs[i];
        const double Z_i = nuc_i.getCharge();
        const mrcpp::Coord<3> &R_i = nuc_i.getCoord();
        for (int j = i + 1; j < nNucs; j++) {
            const Nucleus &nuc_j = nucs[j];
            const double Z_j = nuc_j.getCharge();
            const mrcpp::Coord<3> &R_j = nuc_j.getCoord();
            double R_ij = math_utils::calc_distance(R_i, R_j);
            E_nuc += (Z_i * Z_j) / R_ij;
        }
    }
    return E_nuc;
}

/** @brief Returns the sum of atomic charges*/
double chemistry::get_total_charge(const Nuclei &nucs) {
    double charge = 0;
    for (const auto &nuc : nucs) charge += nuc.getCharge();
    return charge;
}

/** @breif computes the nuclear density as a sum of narrow Gaussians */
Density chemistry::compute_nuclear_density(double prec, const Nuclei &nucs, double alpha) {
    auto beta = std::pow(alpha / MATHCONST::pi, 3.0 / 2.0);
    auto gauss = mrcpp::GaussExp<3>();
    for (auto i = 0; i < nucs.size(); i++) {
        const auto &nuc_i = nucs[i];
        const auto Z_i = nuc_i.getCharge();
        const auto &R_i = nuc_i.getCoord();
        auto gauss_f = mrcpp::GaussFunc<3>(alpha, beta * Z_i, R_i);
        gauss.append(gauss_f);
    }
    auto scaling_factor = (*MRA).getWorldBox().getScalingFactor();
    gauss.makePeriodic(scaling_factor);
    Density rho(false);


    int nPts = 10000;                                // Number of points
    double a[3] = {-20.0, 0.0, 0.0};                 // Start point of plot
    double b[3] = { 20.0, 0.0, 0.0};                 // End point of plot

    mrcpp::Plotter<3> plot;                         // Plotter of 3D functions
    plot.setNPoints(nPts);                          // Set number of points
    plot.setRange(a, b);                            // Set plot range

    density::compute(prec, rho, gauss);
    plot.linePlot(rho.real(), "func_tree_periodic_gauss");



    return rho;
}
/** @breif computes the nuclear density as a sum of narrow Gaussians */
double chemistry::compute_nuclear_self_repulsion(const Nuclei &nucs, double alpha) {
    auto beta = std::pow(alpha / MATHCONST::pi, 3.0 / 2.0);

    auto self_rep = 0.0;
    for (auto i = 0; i < nucs.size(); i++) {
        const auto &nuc_i = nucs[i];
        const auto Z_i = nuc_i.getCharge();
        const auto &R_i = nuc_i.getCoord();
        auto gauss_f = mrcpp::GaussFunc<3>(alpha, beta, R_i);
        self_rep += Z_i * Z_i * gauss_f.calcCoulombEnergy(gauss_f);
    }
    return self_rep;
}
} // namespace mrchem
