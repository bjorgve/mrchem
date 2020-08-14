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
                            bool mpi_share = false);
    ~SmearedNuclearPotential() override { free(NUMBER::Total); }

    void setup(double prec) override { setApplyPrec(prec); }
    void clear() override { clearApplyPrec(); }

private:
    NuclearFunction func;
};

} // namespace mrchem
