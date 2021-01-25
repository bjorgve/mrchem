#pragma once

#include "CoulombPotential.h"
#include "chemistry/Nucleus.h"

namespace mrchem {

class HartreePotential final : public CoulombPotential {
public:
    HartreePotential(std::shared_ptr<mrcpp::PoissonOperator> P,
                     std::shared_ptr<OrbitalVector> Phi,
                     const Nuclei &nucs,
                     const double &rc);
    ~HartreePotential() override = default;

    Density &getBSmear() { return this->b_smeared; }
    Density &getBCorr() { return this->b_corr; }
    Density &getEl() { return this->rho_el; }
    double getRc() { return this->rc; }

private:
    Nuclei nuclei;
    Density b_smeared{false};
    Density b_corr{false};
    Density rho_el{false}; ///< Ground-state electron density
    const double rc;

    void setupGlobalDensity(double prec);
};

} // namespace mrchem
