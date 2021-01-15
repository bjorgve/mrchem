#pragma once

#include "NuclearPotential.h"
#include "analyticfunctions/NuclearFunction.h"
#include "qmoperators/RankZeroTensorOperator.h"

namespace mrchem {

class SmearedNuclearPotential final : public NuclearPotential {
public:
    SmearedNuclearPotential(Nuclei nucs,
                            double proj_prec,
                            double smooth_prec = -1.0,
                            double rc = 1.0,
                            bool mpi_share = false,
                            std::shared_ptr<OrbitalVector> Phi = nullptr);
    ~SmearedNuclearPotential() override { free(NUMBER::Total); }

    void setup(double prec) override;
    // void setup(double  prec ) { setApplyPrec(prec); }
    void clear() override { clearApplyPrec(); }

    std::shared_ptr<QMPotential> getCorrection() override { return corr; }
private:
    std::shared_ptr<QMPotential> corr;
    NuclearFunction func;
    Nuclei nucs;
    double proj_prec;
    double smooth_prec;
    double rc;
    std::shared_ptr<OrbitalVector> orbitals; ///< Unperturbed orbitals defining the ground-state electron density
};

} // namespace mrchem
