#pragma once

#include "chemistry/chemistry_utils.h"
#include "CoulombPotential.h"
#include "CoulombPotentialD1.h"
#include "CoulombPotentialD2.h"
#include "HartreePotential.h"
#include "qmoperators/RankZeroTensorOperator.h"
#include "qmfunctions/qmfunction_utils.h"

/** @class CoulombOperator
 *
 * @brief Operator containing a single CoulombPotential
 *
 * This class is a simple TensorOperator realization of @class CoulombPotential.
 *
 */

namespace mrchem {

class CoulombOperator final : public RankZeroTensorOperator {
public:
    CoulombOperator(mrcpp::PoissonOperator *P)
            : potential(new CoulombPotential(P)) {
        RankZeroTensorOperator &J = (*this);
        J = *this->potential;
    }
    CoulombOperator(mrcpp::PoissonOperator *P, OrbitalVector *Phi)
            : potential(new CoulombPotentialD1(P, Phi)) {
        RankZeroTensorOperator &J = (*this);
        J = *this->potential;
    }
    CoulombOperator(mrcpp::PoissonOperator *P, OrbitalVector *Phi, OrbitalVector *X, OrbitalVector *Y)
            : potential(new CoulombPotentialD2(P, Phi, X, Y)) {
        RankZeroTensorOperator &J = (*this);
        J = *this->potential;
    }
    CoulombOperator(mrcpp::PoissonOperator *P, OrbitalVector *Phi, const Nuclei &nucs)
            : potential(new HartreePotential(P, Phi, nucs)) {
        RankZeroTensorOperator &J = (*this);
        J = *this->potential;
    }
    ~CoulombOperator() {
        if (this->potential != nullptr) delete this->potential;
    }

    const Nuclei &getNuclei() const { return this->potential->getNuclei(); }

    Density &getDensity() {
        if (potential != nullptr) return this->potential->getDensity();
        MSG_FATAL("Coulomb operator not properly initialized");
    }

    ComplexDouble trace(OrbitalVector &Phi) { return 0.5 * RankZeroTensorOperator::trace(Phi); }
    ComplexDouble trace(const Nuclei &nucs) {
        QMFunction &V = *this->potential;
        Density rho_nuc(false);
        chemistry::compute_nuclear_density(this->potential->prec(), rho_nuc, nucs, 1.0e6);
        return 0.5 * qmfunction::dot(V, rho_nuc);
    }

private:
    CoulombPotential *potential;
};

} //namespace mrchem
