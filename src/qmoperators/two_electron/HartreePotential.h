#pragma once

#include "CoulombPotential.h"
#include "chemistry/Nucleus.h"

namespace mrchem {

class HartreePotential final : public CoulombPotential {
public:
    HartreePotential(std::shared_ptr<mrcpp::PoissonOperator> P,
                     std::shared_ptr<OrbitalVector> Phi,
                     const Nuclei &nucs);
    ~HartreePotential() override = default;

private:
    Nuclei nuclei;

    void setupGlobalDensity(double prec);
};

} // namespace mrchem
