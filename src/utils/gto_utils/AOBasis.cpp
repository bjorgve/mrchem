#include "MRCPP/Printer"

#include "AOBasis.h"
#include "AOContraction.h"

using mrcpp::GaussExp;

namespace mrchem {
namespace gto_utils {

AOBasis::AOBasis() {
    this->nPrim = 0;
    this->nFunc = 0;
}

AOBasis::AOBasis(const AOBasis &bas) {
    this->nPrim = 0;
    this->nFunc = 0;
    for (int i = 0; i < bas.size(); i++) append(bas.getContraction(i));
}

AOBasis &AOBasis::operator=(const AOBasis &bas) {
    if (this != &bas) {
        this->nPrim = 0;
        this->nFunc = 0;
        for (int i = 0; i < bas.size(); i++) append(bas.getContraction(i));
    }
    return *this;
}

AOBasis::~AOBasis() {
    for (unsigned int i = 0; i < ctrs.size(); i++) {
        if (this->ctrs[i] != nullptr) delete this->ctrs[i];
    }
}

void AOBasis::append(const AOContraction &ctr) {
    this->ctrs.push_back(new AOContraction(ctr));
    this->nPrim++;
    this->nFunc += ctr.getNComp();
}

GaussExp<3> AOBasis::getAO(int n, const mrcpp::Coord<3> &center) const {
    assert(n >= 0 and n < nFunc);
    int m = 0;
    for (auto i : ctrs) {
        const AOContraction &ctr = *i;
        for (int j = 0; j < ctr.getNComp(); j++) {
            if (m == n) return ctr.getNormContraction(j, center);
            m++;
        }
    }
    MSG_ABORT("Something is terribly wrong");
}

GaussExp<3> AOBasis::getBasis(const mrcpp::Coord<3> &center) const {
    NOT_IMPLEMENTED_ABORT;
    GaussExp<3> abas;
    for (auto i : this->ctrs) {
        const AOContraction &ctr = *i;
        for (int m = 0; m < ctr.getNComp(); m++) abas.append(ctr.getContraction(m, center));
    }
    return abas;
}

GaussExp<3> AOBasis::getNormBasis(const mrcpp::Coord<3> &center) const {
    NOT_IMPLEMENTED_ABORT;
    GaussExp<3> abas;
    for (auto i : this->ctrs) {
        const AOContraction &ctr = *i;
        for (int m = 0; m < ctr.getNComp(); m++) abas.append(ctr.getNormContraction(m, center));
    }
    return abas;
}

} // namespace gto_utils
} // namespace mrchem
