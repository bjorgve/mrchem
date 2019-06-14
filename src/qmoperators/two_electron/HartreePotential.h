#pragma once

#include "CoulombPotential.h"
#include "chemistry/Nucleus.h"

namespace mrchem {

class HartreePotential final : public CoulombPotential {
public:
    HartreePotential(std::shared_ptr<mrcpp::PoissonOperator> P,
                     std::shared_ptr<OrbitalVector> Phi,
                     const Nuclei &nucs,
                     double prec);
    ~HartreePotential() override = default;

private:
    void setupGlobalDensity(double prec);
};

} // namespace mrchem
