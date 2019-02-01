#pragma once

#include "chemistry/Nucleus.h"
#include "qmoperators/two_electron/CoulombPotential.h"
/** @class HartrePotential
 *
 * @brief Coulomb potential defined by a particular electron density
 *
 * The Coulomb potential is computed by application of the Poisson operator
 * an electron density. There are two ways of defining the density:
 *
 *  1) Use getDensity() prior to setup() and build the density as you like.
 *  2) Provide a default set of orbitals in the constructor that is used to
 *     compute the density on-the-fly in setup().
 *
 * If a set of orbitals has NOT been given in the constructor, the density
 * MUST be explicitly computed prior to setup(). The density will be computed
 * on-the-fly in setup() ONLY if it is not already available. After setup() the
 * operator will be fixed until clear(), which deletes both the density and the
 * potential.
 */

namespace mrchem {

class HartreePotential final : public CoulombPotential {
public:
    HartreePotential(mrcpp::PoissonOperator *P, OrbitalVector *Phi, const Nuclei &nucs, double prec);
    virtual ~HartreePotential() = default;

private:
    OrbitalVector *orbitals; ///< Unperturbed orbitals defining the ground-state electron density

    void setupDensity(double prec);
};

} //namespace mrchem
