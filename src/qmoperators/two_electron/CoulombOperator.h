#pragma once

#include "MRCPP/Printer"

#include "CoulombPotential.h"
#include "CoulombPotentialD1.h"
#include "CoulombPotentialD2.h"
#include "HartreePotential.h"
#include "chemistry/chemistry_utils.h"
#include "qmfunctions/qmfunction_utils.h"
#include "qmoperators/RankZeroTensorOperator.h"

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
    CoulombOperator(std::shared_ptr<mrcpp::PoissonOperator> P) {
        potential = std::make_shared<CoulombPotential>(P);

        // Invoke operator= to assign *this operator
        RankZeroTensorOperator &J = (*this);
        J = potential;
    }
    CoulombOperator(std::shared_ptr<mrcpp::PoissonOperator> P, std::shared_ptr<OrbitalVector> Phi) {
        potential = std::make_shared<CoulombPotentialD1>(P, Phi);

        // Invoke operator= to assign *this operator
        RankZeroTensorOperator &J = (*this);
        J = potential;
    }
    CoulombOperator(std::shared_ptr<mrcpp::PoissonOperator> P,
                    std::shared_ptr<OrbitalVector> Phi,
                    std::shared_ptr<OrbitalVector> X,
                    std::shared_ptr<OrbitalVector> Y) {
        potential = std::make_shared<CoulombPotentialD2>(P, Phi, X, Y);

        // Invoke operator= to assign *this operator
        RankZeroTensorOperator &J = (*this);
        J = potential;
    }

    CoulombOperator(std::shared_ptr<mrcpp::PoissonOperator> P, std::shared_ptr<OrbitalVector> Phi, const Nuclei &nucs) {

        println(0, "Couloumb nuclei size prior to init " << nucs.size()) potential =
            std::make_shared<HartreePotential>(P, Phi, nucs);
        println(0, "Couloumb nuclei size " << (*potential).getNuclei().size()) RankZeroTensorOperator &J = (*this);
        J = potential;
    }
    const Nuclei &getNuclei() const { return this->potential->getNuclei(); }

    ~CoulombOperator() override = default;

    auto &getPoisson() { return this->potential->getPoisson(); }
    auto &getDensity() { return this->potential->getDensity(); }

    ComplexDouble trace(OrbitalVector &Phi) { return 0.5 * RankZeroTensorOperator::trace(Phi); }
    ComplexDouble trace(const Nuclei &nucs) {
        QMFunction &V = *this->potential;
        double nuc_prec = 10e-4; // this->potential->getNucPrec();
        Density rho_nuc(false);
        rho_nuc = chemistry::compute_nuclear_density(this->potential->prec(), nucs, 1.0 / nuc_prec);
        println(0, "Colomb trace, rho_nuc.integrate() " << rho_nuc.integrate())
                println(0, "Colomb trace, rho_nuc.norm() " << rho_nuc.norm())
                    println(0, "Colomb trace, V.integrate() " << V.integrate())
                        println(0, "Colomb trace, V.norm() " << V.norm())
            // println(0, qmfunction::dot(V, rho_nuc).real().integrate())
            return 0.5 *
            qmfunction::dot(V, rho_nuc);
    }

private:
    std::shared_ptr<CoulombPotential> potential{nullptr};
};

} // namespace mrchem
