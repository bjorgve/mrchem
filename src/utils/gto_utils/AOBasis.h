#pragma once

#include "AOContraction.h"
#include <vector>

namespace mrchem {
namespace gto_utils {

class AOBasis final {
public:
    AOBasis();
    AOBasis(const AOBasis &bas);
    AOBasis &operator=(const AOBasis &bas);
    ~AOBasis();

    void append(const AOContraction &ctr);

    mrcpp::GaussExp<3> getAO(int n, const mrcpp::Coord<3> &center) const;
    mrcpp::GaussExp<3> getBasis(const mrcpp::Coord<3> &center) const;
    mrcpp::GaussExp<3> getNormBasis(const mrcpp::Coord<3> &center) const;

    AOContraction &getContraction(int n) { return *this->ctrs[n]; }
    const AOContraction &getContraction(int n) const { return *this->ctrs[n]; }

    int size() const { return this->ctrs.size(); }
    int getNFunc() const { return this->nFunc; }

    // This should print shell by shell
    friend std::ostream &operator<<(std::ostream &o, const AOBasis &b) {
        o << "    nFunc " << b.nFunc << std::endl;
        for (auto ctr : b.ctrs) o << "    " << *ctr;
        return o;
    }

private:
    int nPrim; ///< Total number of primitives in set
    int nFunc; ///< Total number of functions (all l-components included)
    std::vector<AOContraction *> ctrs;
};

} // namespace gto_utils
} // namespace mrchem
