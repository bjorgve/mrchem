#include "MRCPP/Printer"
#include "MRCPP/Timer"

#include "CoulombOperator.h"
#include "ExchangeOperator.h"
#include "FarFieldOperator.h"
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
using FarFieldOperator_p = std::shared_ptr<mrchem::FarFieldOperator>;
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
                           FarFieldOperator_p f,
                           ExchangeOperator_p k,
                           XCOperator_p xc,
                           ElectricFieldOperator_p ext)
        : kin(t)
        , nuc(v)
        , coul(j)
        , ff(f)
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
    if (this->ff != nullptr) this->V += (*this->ff);
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
    Timer t_tot;
    auto plevel = Printer::getPrintLevel();
    mrcpp::print::header(2, "Building Fock operator");
    mrcpp::print::value(2, "Precision", prec, "(rel)", 5);
    mrcpp::print::separator(2, '-');
    this->kinetic().setup(prec);
    this->potential().setup(prec);
    this->perturbation().setup(prec);
    t_tot.stop();
    mrcpp::print::footer(2, t_tot, 2);
    if (plevel == 1) mrcpp::print::time(1, "Building Fock operator", t_tot);
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
SCFEnergy FockOperator::trace(OrbitalVector &Phi, const Nuclei &nucs) {
    Timer t_tot;
    auto plevel = Printer::getPrintLevel();
    mrcpp::print::header(2, "Computing molecular energy");

    double E_kin = 0.0;  // Kinetic energy
    double E_nn = 0.0;   // Nuclear repulsion
    double E_en = 0.0;   // Nuclear-electronic interaction
    double E_ee = 0.0;   // Electronic repulsion
    double E_x = 0.0;    // Exact Exchange
    double E_xc = 0.0;   // Exchange and Correlation
    double E_eext = 0.0; // External field contribution to the electronic energy
    double E_next = 0.0; // External field contribution to the nuclear energy
    double E_eeff = 0.0; // Far field contribution to the electronic energy
    double E_nff = 0.0;  // Far field contribution to the nuclear energy
    double E_en_hack = 0.0;  // Far field contribution to the nuclear energy

    //// PERIODIC ENERGIES ////
    double E_se = 0.0;
    double E_ncor = 0.0;

    // Nuclear part
    if (this->nuc != nullptr) E_nn = chemistry::compute_nuclear_repulsion(nucs);
    if (this->ext != nullptr) E_next = -this->ext->trace(nucs).real();
    if (this->ff != nullptr) E_nff = -0.5 * this->ff->trace(nucs).real();

    // Electronic part
    if (this->kin != nullptr) E_kin = this->kin->trace(Phi).real();
    if (this->nuc != nullptr) E_en = this->nuc->trace(Phi).real();
    if (this->nuc == nullptr) println(0, "nullpointer fucker")
    if (this->nuc != nullptr) println(0, "ikke nullpointer fucker")
    println(0, "her da?")
    if (this->coul != nullptr) E_ee = 0.5 * this->coul->trace(Phi).real();
    if (this->ex != nullptr) E_x = -this->exact_exchange * this->ex->trace(Phi).real();
    if (this->xc != nullptr) E_xc = this->xc->getEnergy();
    if (this->ext != nullptr) E_eext = this->ext->trace(Phi).real();
    if (this->ff != nullptr) E_eeff = this->ff->trace(Phi).real();




    if (this->coul->getPotential()->getRc() > 0.0) {
        // RankZeroTensorOperator hack;
        // hack += this->nuc->getCorrection();
        // // if (this->nuc != nullptr) E_en_hack = 0.5*hack.trace(Phi).real();
        // // if (this->nuc != nullptr) E_en_hack = hack.trace(Phi).real();
        // E_nn = 0.0;
        // E_se = 0.5 * mrcpp::dot<3>(this->coul->getPotential()->getBSmear().real(), this->coul->getPotential()->real());
        // auto rc = this->coul->getPotential()->getRc();
        // // rc = rc*0.1;
        // if (rc <= 0.0) MSG_ABORT("RC has to be positive");
        // auto Ig = 10976.0 / (17875.0 * rc);
        // E_en -= 2.0*E_en_hack;
        // for (auto &nuc : nucs) { E_ncor += 0.5 * nuc.getCharge() * nuc.getCharge() * (Ig - 12.0 / (5.0 * rc)); }
        // for (auto &nuc : nucs) { E_ncor -= 0.5 * nuc.getCharge() * nuc.getCharge() * (Ig + 12.0 / (5.0 * rc)); }

        // /// Equation 16:
        //
        // auto eq_16 = 0.0;
        // for (auto &nuc : nucs) eq_16 += nuc.getCharge() * nuc.getCharge() * Ig;
        // println(0, "IG part " << eq_16)
        // println(0, "E_en " << E_en)
        // eq_16 += E_en;
        // eq_16 += hack.trace(Phi).real();
        // println(0, "Equation 16 " << eq_16)
        // println(0, "SUM " << E_ee + E_se + E_ncor)
        //
        // // Equation 11
        // auto eq_11 = -this->coul->trace(nucs).real() - mrcpp::dot(this->coul->getPotential()->getBSmear().real(), this->coul->getPotential()->real());
        // println(0, "11 a " << this->coul->trace(nucs).real())
        // println(0, "11 b " << mrcpp::dot(this->coul->getPotential()->getBSmear().real(), this->coul->getPotential()->real()))
        // println(0, "Equation 11 " << eq_11)


        // auto E_pp = 0.0;
        // for (auto &nuc : nucs) E_pp -= 0.5*nuc.getCharge() * nuc.getCharge() * (Ig + 12.0 / (5.0 * rc));
        //
        // auto E_pm = -this->coul->trace(nucs).real() - mrcpp::dot(this->coul->getPotential()->getBSmear().real(), this->coul->getPotential()->real());        // auto dip_oper = H_E_dip({0.0, 0.0, 0.0});
        //
        // auto E_mm_a = 0.5*mrcpp::dot(this->coul->getPotential()->getBSmear().real(), this->coul->getPotential()->real());
        // auto E_mm_b = 0.5*this->coul->trace(Phi).real();
        // auto E_mm = E_mm_a + E_mm_b;
        // println(0, "E_pp " << E_pp)
        // println(0, "E_pm " << E_pm)
        //
        // println(0, "E_mm_a " << E_mm_a)
        // println(0, "E_mm_b " << E_mm_b)
        // println(0, "E_mm " << E_mm)
        //
        // println(0, "SUM " << E_pp + E_pm + E_mm)
        // dip_oper.setup(1.0e-6);
        // DoubleVector nuc_dip = -dip_oper.trace(nucs).real();
        // DoubleVector dip_el = dip_oper.trace(Phi).real();
        // DoubleVector tot_dip = nuc_dip + dip_el;
        //
        // auto new_charge = tot_dip[2]/8.0;
        //
        // auto rc_tmp = rc;
        //
        // mrcpp::Coord<3> pos_coord{0.0, 0.0, 4.};
        // mrcpp::Coord<3> neg_coord{0.0, 0.0, -4.};
        // Ig = 10976.0 / (17875.0 * rc_tmp);
        //
        // dip_oper.clear();
        // println(0, 0.5 * new_charge * new_charge* Ig)    //- 12.0 / (5.0 * rc_tmp)));
        // E_ncor += 0.5 * new_charge * new_charge* (Ig - 12.0 / (5.0 * rc_tmp));

    }

    ///////////
    std::cout << std::setprecision(16) << std::fixed;

    // auto test = -0.5*mrcpp::dot(this->coul->getPotential()->getEl().real(), this->coul->getPotential()->real());
    // println(0, "test " << test)
    // auto test_2 = -0.5*this->coul->trace(Phi).real();
    // // println(0, "test 2 " << test_2)
    // auto E_mm = 0.5*mrcpp::dot(this->coul->getPotential()->getDensity().real(), this->coul->getPotential()->real());
    // // println(0, "E_mm" << E_mm)
    // auto E_mp = mrcpp::dot(this->coul->getPotential()->getDensity().real(), this->nuc->getPotential()->real());
    // println(0, "E_mp " <<  E_mp)
    // // println(0, "E_mm + E_mp " << E_mm + E_mm)
    // auto E_pp_dir = mrcpp::dot(this->coul->getPotential()->getBCorr().real(), this->nuc->getPotential()->real());
    // E_pp_dir += chemistry::compute_nuclear_repulsion(nucs);
    // E_pp_dir += mrcpp::dot(this->coul->getPotential()->getBCorr().real(), this->nuc->getCorrection()->real());
    // println(0, "E_pp_dir " << E_pp_dir)
    // //// E_pp
    // RankZeroTensorOperator hack;
    // hack += this->nuc->getCorrection();
    // auto rc = this->coul->getPotential()->getRc();
    // auto Ig = 10976.0 / (17875.0 * rc);

    // auto eq_1 = 0.0;
    // for (auto &nuc : nucs) { eq_1 -= 0.5 * nuc.getCharge() * nuc.getCharge() * (Ig + 12.0 / (5.0 * rc)); }


    // auto eq_2_a = -0.5*hack.trace(nucs).real();
    // auto eq_2_b = 0.5*mrcpp::dot(this->coul->getPotential()->getBSmear().real(), this->nuc->getCorrection()->real());
    //
    // auto eq_3 = 0.5*mrcpp::dot(this->coul->getPotential()->getBCorr().real(), this->nuc->getCorrection()->real());

    ///
    // auto eq_2_a = -0.5*hack.trace(nucs).real();
    // auto eq_2_b = -0.5*mrcpp::dot(this->coul->getPotential()->getBSmear().real(), this->nuc->getCorrection()->real());
    // auto eq_3 = -0.5*mrcpp::dot(this->coul->getPotential()->getBCorr().real(), this->nuc->getNoCorr()->real());
    // // println(0, "show stig " << -0.5*mrcpp::dot(this->coul->getPotential()->getBCorr().real(), this->nuc->getPotential()->real()))
    // println(0, "eq_2_a " << eq_2_a)
    // println(0, "eq_2_b " << eq_2_b)
    // println(0, "eq_2_a + eq_2_b " << eq_2_a + eq_2_b)
    // println(0, "eq_3 " << eq_3)
    // // auto eq_3 = -0.5*mrcpp::dot(this->coul->getPotential()->getBCorr().real(), this->nuc->getPotential()->real());
    // auto eq_4 = 0.5*mrcpp::dot(this->coul->getPotential()->getBCorr().real(), this->nuc->getCorrection()->real());
    // println(0, "eq_4 " << eq_4)
    // // auto E_pp = eq_1 + eq_2_a + eq_2_b + eq_3 + eq_4;
    // auto E_pp = eq_2_a + eq_2_b + eq_3 + eq_4;
    // println(0, "E_pp " << E_pp )
    //
    // println(0, "Da sum " << E_mm + E_mp + E_pp)
    // println(0, "Da sum v2" << E_mm + E_mp + E_pp + E_kin + E_xc)


    //
    // println(0, "new_e_pp " << eq_2+eq_3+eq_4+eq_1)
    // println(0, "SUM new_e_pp " << eq_2+eq_3+eq_4+eq_1+E_mm+E_mp)
    // println(0, "SUM_2 new_e_pp " << eq_2+eq_3+eq_4+eq_1+E_mm+E_mp+E_kin+E_xc)

    /// Investigate E_pm
    // auto e_pm_2 = mrcpp::dot(this->coul->getPotential()->getBCorr().real(), this->coul->getPotential()->real());
    // auto e_pm_1_a = -this->coul->trace(nucs).real();
    // auto e_pm_1_b = -mrcpp::dot(this->coul->getPotential()->getBSmear().real(), this->coul->getPotential()->real());
    // println(0, "e_pm_2 " << e_pm_2)
    // println(0, "e_pm_1_a "<< e_pm_1_a)
    // println(0, "e_pm_1_b "<< e_pm_1_b)
    // println(0, "sum " << e_pm_2 + e_pm_1_a + e_pm_1_b)
    //
    // // auto E_pp = eq_1 + eq_2_a + eq_2_b + eq_3;
    // // println(0, "E_pp " << E_pp)
    // //
    // // println(0, "eq_1 " << eq_1)
    // // println(0, "eq_2_a " << eq_2_a)
    // // println(0, "eq_2_b " << eq_2_b)
    // // println(0, "eq_3 " << eq_3)
    // // println(0, "SUM " << E_mm + E_mp + E_pp)
    //
    // //////////
    // auto E_pp = 0.0;
    // for (auto &nuc : nucs) E_pp -= 0.5*nuc.getCharge() * nuc.getCharge() * (Ig + 12.0 / (5.0 * rc));
    //
    // auto E_pm = -this->coul->trace(nucs).real() - mrcpp::dot(this->coul->getPotential()->getBSmear().real(), this->coul->getPotential()->real());        // auto dip_oper = H_E_dip({0.0, 0.0, 0.0});
    //                   mrcpp::dot(this->coul->getPotential()->getBSmear().real(), this->coul->getPotential()->real())
    // auto E_mm_a = 0.5*mrcpp::dot(this->coul->getPotential()->getBSmear().real(), this->coul->getPotential()->real());
    // auto E_mm_b = 0.5*this->coul->trace(Phi).real();
    // auto E_mm = E_mm_a + E_mm_b;
    // auto stig_1 = -this->coul->trace(nucs).real();
    // auto stig_2 = -hack.trace(nucs).real();
    // auto stig_3 = 0.5*this->coul->trace(Phi).real();
    // auto stig_4 = -0.5*mrcpp::dot(this->coul->getPotential()->getBSmear().real(), this->coul->getPotential()->real());
    // auto stig_5 = -0.5*mrcpp::dot(this->coul->getPotential()->getBCorr().real(), this->coul->getPotential()->real());
    // auto stig_6 = 0.5*hack.trace(Phi).real();
    // auto stig_7 = -0.5*mrcpp::dot(this->coul->getPotential()->getBSmear().real(), this->nuc->getCorrection()->real());
    //
    // auto E_pm = stig_1 + stig_4;
    // auto E_mm =
    // auto stig_8 = 0.0;
    // for (auto &nuc : nucs) { stig_8 -= 0.5 * nuc.getCharge() * nuc.getCharge() * (Ig + 12.0 / (5.0 * rc)); }
    // println(0, "stig_1 " << stig_1)
    // println(0, "stig_2 " << stig_2)
    // println(0, "stig_3 " << stig_3)
    // println(0, "stig_4 " << stig_4)
    // println(0, "stig_5 " << stig_5)
    // println(0, "stig_6 " << stig_6)
    // println(0, "stig_7 " << stig_7)
    // println(0, "stig_8 " << stig_8)
    //
    // auto stig = stig_1 + stig_2 + stig_3 + stig_4 + stig_5 + stig_6 + stig_7 + stig_8;
    // println(0, "stig " << stig)
    // println(0, "tot enrgy " << stig + E_xc + E_kin)
    // println(0, "stig minus 16 " << stig - 0.5*stig_1 - 0.5*stig_4)

    mrcpp::print::footer(2, t_tot, 2);
    if (plevel == 1) mrcpp::print::time(1, "Computing molecular energy", t_tot);

    return SCFEnergy{E_kin, E_nn, E_en, E_ee, E_x, E_xc, E_next, E_eext, E_nff, E_eeff, E_se, E_ncor};
}

ComplexMatrix FockOperator::operator()(OrbitalVector &bra, OrbitalVector &ket) {
    Timer t_tot;
    auto plevel = Printer::getPrintLevel();
    mrcpp::print::header(2, "Computing Fock matrix");

    auto t = this->getKineticOperator();
    auto v = this->potential();

    ComplexMatrix T = ComplexMatrix::Zero(bra.size(), ket.size());
    if (t != nullptr) T += (*t)(bra, ket);

    ComplexMatrix V = ComplexMatrix::Zero(bra.size(), ket.size());
    if (v.size() > 0) V += v(bra, ket);

    mrcpp::print::footer(2, t_tot, 2);
    if (plevel == 1) mrcpp::print::time(1, "Computing Fock matrix", t_tot);
    return T + V;
}

ComplexMatrix FockOperator::dagger(OrbitalVector &bra, OrbitalVector &ket) {
    NOT_IMPLEMENTED_ABORT;
}

} // namespace mrchem
