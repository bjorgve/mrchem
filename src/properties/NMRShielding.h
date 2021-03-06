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
class NMRShielding final {
public:
    NMRShielding(int k, const Nucleus &n) : K(k), nuc(n) {}

    int getK() const { return this->K; }
    const Nucleus &getNucleus() const { return this->nuc; }

    DoubleMatrix getTensor() const { return getDiamagnetic() + getParamagnetic(); }
    DoubleMatrix &getDiamagnetic() { return this->dia_tensor; }
    DoubleMatrix &getParamagnetic() { return this->para_tensor; }
    const DoubleMatrix &getDiamagnetic() const { return this->dia_tensor; }
    const DoubleMatrix &getParamagnetic() const { return this->para_tensor; }

    void print() const {
        auto iso_ppm_d = getDiamagnetic().trace() / 3.0;
        auto iso_ppm_p = getParamagnetic().trace() / 3.0;
        auto iso_ppm_t = iso_ppm_d + iso_ppm_p;

        mrcpp::print::header(0, "NMR shielding");
        print_utils::scalar(0, "Nucleus K", getK(), getNucleus().getElement().getSymbol(), 0);
        mrcpp::print::separator(0, '-');
        print_utils::matrix(0, "Diamagnetic", getDiamagnetic());
        print_utils::scalar(0, "Isotropic average", iso_ppm_d, "(ppm)");
        mrcpp::print::separator(0, '-');
        print_utils::matrix(0, "Paramagnetic", getParamagnetic(), -1);
        print_utils::scalar(0, "Isotropic average", iso_ppm_p, "(ppm)");
        mrcpp::print::separator(0, '-');
        print_utils::matrix(0, "Total tensor", getTensor());
        print_utils::scalar(0, "Isotropic average", iso_ppm_t, "(ppm)");
        mrcpp::print::separator(0, '=', 2);
    }

private:
    const int K;
    const Nucleus nuc;
    DoubleMatrix dia_tensor{math_utils::init_nan(3,3)};
    DoubleMatrix para_tensor{math_utils::init_nan(3,3)};
};
// clang-format on

} // namespace mrchem
