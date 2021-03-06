/*
 * MRChem, a numerical real-space code for molecular electronic structure
 * calculations within the self-consistent field (SCF) approximations of quantum
 * chemistry (Hartree-Fock and Density Functional Theory).
 * Copyright (C) 2020 Stig Rune Jensen, Luca Frediani, Peter Wind and contributors.
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

#pragma once

#include "mrchem.h"

#include "utils/math_utils.h"
#include "utils/print_utils.h"

namespace mrchem {

// clang-format off
class Magnetizability final {
public:
    DoubleMatrix getTensor() const { return getDiamagnetic() + getParamagnetic(); }
    DoubleMatrix &getDiamagnetic() { return this->dia_tensor; }
    DoubleMatrix &getParamagnetic() { return this->para_tensor; }
    const DoubleMatrix &getDiamagnetic() const { return this->dia_tensor; }
    const DoubleMatrix &getParamagnetic() const { return this->para_tensor; }

    void print() const {
        auto iso_au_d = getDiamagnetic().trace() / 3.0;
        auto iso_au_p = getParamagnetic().trace() / 3.0;
        auto iso_au_t = iso_au_d + iso_au_p;

        // SI units (J/T^2 10^{-30})
        auto iso_si_t = iso_au_t * PHYSCONST::JT_m2;
        auto iso_si_d = iso_au_d * PHYSCONST::JT_m2;
        auto iso_si_p = iso_au_p * PHYSCONST::JT_m2;

        mrcpp::print::header(0, "Magnetizability");
        print_utils::matrix(0, "Diamagnetic", getDiamagnetic());
        print_utils::scalar(0, "Isotropic average", iso_au_d, "(au)");
        print_utils::scalar(0, "                 ", iso_si_d, "(SI)");
        mrcpp::print::separator(0, '-');
        print_utils::matrix(0, "Paramagnetic", getParamagnetic());
        print_utils::scalar(0, "Isotropic average", iso_au_p, "(au)");
        print_utils::scalar(0, "                 ", iso_si_p, "(SI)");
        mrcpp::print::separator(0, '-');
        print_utils::matrix(0, "Total tensor", getTensor());
        print_utils::scalar(0, "Isotropic average", iso_au_t, "(au)");
        print_utils::scalar(0, "                 ", iso_si_t, "(SI)");
        mrcpp::print::separator(0, '=', 2);
    }

private:
    DoubleMatrix dia_tensor{math_utils::init_nan(3,3)};
    DoubleMatrix para_tensor{math_utils::init_nan(3,3)};
};
// clang-format on

} // namespace mrchem
