#include "MRCPP/Printer"
#include "MRCPP/Timer"

#include "CoulombOperator.h"
#include "ExchangeOperator.h"
#include "FockOperator.h"
#include "XCOperator.h"
#include "chemistry/chemistry_utils.h"
#include "properties/SCFEnergy.h"
#include "qmfunctions/Orbital.h"
#include "qmfunctions/orbital_utils.h"
#include "qmoperators/one_electron/ElectricFieldOperator.h"
#include "qmoperators/one_electron/KineticOperator.h"
#include "qmoperators/one_electron/NuclearOperator.h"
#include "utils/math_utils.h"

using mrcpp::Printer;
using mrcpp::Timer;

using KineticOperator_p = std::shared_ptr<mrchem::KineticOperator>;
using NuclearOperator_p = std::shared_ptr<mrchem::NuclearOperator>;
using CoulombOperator_p = std::shared_ptr<mrchem::CoulombOperator>;
using ExchangeOperator_p = std::shared_ptr<mrchem::ExchangeOperator>;
using XCOperator_p = std::shared_ptr<mrchem::XCOperator>;
using ElectricFieldOperator_p = std::shared_ptr<mrchem::ElectricFieldOperator>;

namespace mrchem {
extern mrcpp::MultiResolutionAnalysis<3> *MRA; // Global MRA

/** @brief constructor
 *
 * @param t:   Kinetic operator
 * @param v:   Nuclear potential operator
 * @param j:   Coulomb operator
 * @param k:   HF exchange operator
 * @param xc:  Exchange-Correlation operator
 * @param ext: External Field operator
 *
 * Each of the arguments can be NULL, so this operators includes both core Hamiltonian,
 * the Hartree(-Fock) method and (pure/hybrid) Density Functional Theory.
 */
FockOperator::FockOperator(KineticOperator_p t,
                           NuclearOperator_p v,
                           CoulombOperator_p j,
                           ExchangeOperator_p k,
                           XCOperator_p xc,
                           ElectricFieldOperator_p ext)
        : kin(t)
        , nuc(v)
        , coul(j)
        , ex(k)
        , xc(xc)
        , ext(ext) {}

/** @brief build the Fock operator once all contributions are in place
 *
 */
void FockOperator::build(double exx) {
    this->exact_exchange = exx;
    this->T = RankZeroTensorOperator();
    if (this->kin != nullptr) this->T += (*this->kin);

    this->V = RankZeroTensorOperator();
    if (this->nuc != nullptr) this->V += (*this->nuc);
    if (this->coul != nullptr) this->V += (*this->coul);
    if (this->ex != nullptr) this->V -= this->exact_exchange * (*this->ex);
    if (this->xc != nullptr) this->V += (*this->xc);
    if (this->ext != nullptr) this->V += (*this->ext);

    RankZeroTensorOperator &F = (*this);
    F = this->kinetic() + this->potential();
}

/** @brief prepare operator for application
 *
 * @param prec: apply precision
 *
 * This will call the setup function of all underlying operators, and in particular
 * it will compute the internal exchange if there is an ExchangeOperator.
 */
void FockOperator::setup(double prec) {
    Timer timer;
    Printer::printHeader(0, "Setting up Fock operator");
    Printer::printDouble(0, "Precision", prec, 5);
    Printer::printSeparator(0, '-');
    this->kinetic().setup(prec);
    this->potential().setup(prec);
    this->perturbation().setup(prec);
    if (this->ex != nullptr) this->ex->setupInternal(prec);
    timer.stop();
    Printer::printFooter(0, timer, 2);
}

/** @brief clear operator after application
 *
 * This will call the clear function of all underlying operators, and bring them back
 * to the state after construction. The operator can now be reused after another setup.
 */
void FockOperator::clear() {
    this->kinetic().clear();
    this->potential().clear();
    this->perturbation().clear();
}

/** @brief rotate orbitals of two-electron operators
 *
 * @param U: unitary transformation matrix
 *
 * This function should be used in case the orbitals are rotated *after* the FockOperator
 * has been setup. In particular the ExchangeOperator needs to rotate the precomputed
 * internal exchange potentials.
 */
void FockOperator::rotate(const ComplexMatrix &U) {
    if (this->ex != nullptr) this->ex->rotate(U);
}

/** @brief compute the SCF energy
 *
 * @param Phi: orbitals
 * @param F: Fock matrix
 *
 * This function will compute the total energy for a given OrbitalVector and
 * the corresponding Fock matrix. Tracing the kinetic energy operator is avoided
 * by tracing the Fock matrix and subtracting all other contributions.
 */
SCFEnergy FockOperator::trace(OrbitalVector &Phi, const ComplexMatrix &F) {
    auto E_nuc = 0.0; // Nuclear repulsion
    auto E_el = 0.0;  // Electronic energy
    auto E_orb = 0.0; // Orbital energy
    auto E_kin = 0.0; // Kinetic energy
    auto E_en = 0.0;  // Nuclear-electronic interaction
    auto E_ee = 0.0;  // Electronic repulsion
    auto E_x = 0.0;   // Exact Exchange
    auto E_xc = 0.0;  // Exchange and Correlation
    auto E_xc2 = 0.0; // Trace of the XC operator
    auto E_ext = 0.0; // External field contribution to the electronic energy
    auto E_nex = 0.0; // External field contribution to the nuclear energy

    // Nuclear part
    if (this->nuc != nullptr) {
        Nuclei &nucs = this->nuc->getNuclei();
        E_nuc = chemistry::compute_nuclear_repulsion(nucs);
        if (this->ext != nullptr) {
            E_nex = this->ext->trace(nucs).real();
            E_nuc += E_nex;
        }
    }

    // Orbital energies
    for (auto i = 0; i < Phi.size(); i++) {
        auto occ = static_cast<double>(Phi[i].occ());
        E_orb += occ * F(i, i).real();
    }

    // Electronic part
    if (this->nuc != nullptr) E_en = this->nuc->trace(Phi).real();
    if (this->coul != nullptr) E_ee = this->coul->trace(Phi).real();
    if (this->ex != nullptr) E_x = -this->exact_exchange * this->ex->trace(Phi).real();
    if (this->xc != nullptr) E_xc = this->xc->getEnergy();
    if (this->xc != nullptr) E_xc2 = this->xc->trace(Phi).real();
    if (this->ext != nullptr) E_ext = this->ext->trace(Phi).real();

    println(0, "E_en " << E_en);
    println(0, "E_ee " << E_ee);
    println(0, "E_x " << E_xc);
    println(0, "E_xc2 " << E_xc2);
    println(0, "E_ext " << E_ext);
    auto E_eex = E_ee + E_x;
    println(0, "E_eex " << E_eex);
    auto E_orbxc2 = E_orb - E_xc2;
    E_kin = E_orbxc2 - 2.0 * E_eex - E_en - E_ext;
    println(0, "E_kin " << E_kin);
    E_el = E_orbxc2 - E_eex + E_xc;
    println(0, "E_el " << E_el);
    return SCFEnergy{E_nuc, E_el, E_orb, E_kin, E_en, E_ee, E_xc, E_x, E_nex, E_ext};
}

ComplexMatrix FockOperator::operator()(OrbitalVector &bra, OrbitalVector &ket) {
    Timer t_tot;
    Printer::printHeader(0, "Calculating Fock matrix");

    auto t = this->getKineticOperator();
    auto v = this->potential();

    Timer t_kin;
    ComplexMatrix T = ComplexMatrix::Zero(bra.size(), ket.size());
    if (t != nullptr) T += (*t)(bra, ket);
    t_kin.stop();
    Printer::printDouble(0, "Kinetic part", t_kin.getWallTime());

    Timer t_pot;
    ComplexMatrix V = ComplexMatrix::Zero(bra.size(), ket.size());
    if (v.size() > 0) V += v(bra, ket);
    t_pot.stop();
    Printer::printDouble(0, "Potential part", t_pot.getWallTime());

    t_tot.stop();
    Printer::printFooter(0, t_tot, 2);
    return T + V;
}

ComplexMatrix FockOperator::dagger(OrbitalVector &bra, OrbitalVector &ket) {
    NOT_IMPLEMENTED_ABORT;
}

} // namespace mrchem
