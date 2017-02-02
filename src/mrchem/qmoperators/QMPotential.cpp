#include "QMPotential.h"
#include "MWAdder.h"
#include "MWMultiplier.h"
#include "GridGenerator.h"
#include "Orbital.h"
#include "Timer.h"

extern MultiResolutionAnalysis<3> *MRA; // Global MRA

using namespace std;

QMPotential::QMPotential()
    : ComplexFunction<3>(0, 0),
      QMOperator(MRA->getMaxScale()) {
}

QMPotential::~QMPotential() {
    if (this->hasReal()) MSG_ERROR("Potential not cleared");
    if (this->hasImag()) MSG_ERROR("Potential not cleared");
}

void QMPotential::setup(double prec) {
    QMOperator::setup(prec);
    if (this->hasReal()) MSG_ERROR("Potential not cleared");
    if (this->hasImag()) MSG_ERROR("Potential not cleared");
}

void QMPotential::clear() {
    clearReal(true);
    clearImag(true);
    QMOperator::clear();
}

Orbital* QMPotential::operator() (Orbital &phi) {
    Orbital *Vphi = new Orbital(phi);

    Timer timer;
    calcRealPart(*Vphi, phi);
    calcImagPart(*Vphi, phi, false);
    timer.stop();

    int n = Vphi->getNNodes();
    double t = timer.getWallTime();
    TelePrompter::printTree(1, "Applied QM potential", n, t);

    return Vphi;
}

Orbital* QMPotential::adjoint(Orbital &phi) {
    Orbital *Vphi = new Orbital(phi);

    Timer timer;
    calcRealPart(*Vphi, phi);
    calcImagPart(*Vphi, phi, true);
    timer.stop();

    int n = Vphi->getNNodes();
    double t = timer.getWallTime();
    TelePrompter::printTree(1, "Applied QM adjoint potential", n, t);

    return Vphi;
}

void QMPotential::calcRealPart(Orbital &Vphi, Orbital &phi) {
    if (this->apply_prec < 0.0) MSG_ERROR("Uninitialized operator");

    MWAdder<3> add(this->apply_prec, this->max_scale);
    MWMultiplier<3> mult(this->apply_prec, this->max_scale);
    GridGenerator<3> grid(this->max_scale);

    QMPotential &V = *this;
    FunctionTreeVector<3> vec;

    if (Vphi.hasReal()) MSG_ERROR("Orbital not empty");
    if (V.hasReal() and phi.hasReal()) {
        FunctionTree<3> *tree = new FunctionTree<3>(*MRA);
        grid(*tree, phi.real());
        mult(*tree, 1.0, V.real(), phi.real(), 0);
        vec.push_back(1.0, tree);
    }
    if (V.hasImag() and phi.hasImag()) {
        FunctionTree<3> *tree = new FunctionTree<3>(*MRA);
        grid(*tree, phi.imag());
        mult(*tree, 1.0, V.imag(), phi.imag(), 0);
        vec.push_back(-1.0, tree);
    }
    if (vec.size() == 1) {
        Vphi.setReal(vec[0]);
        vec.clear(false);
    }
    if (vec.size() == 2) {
        Vphi.allocReal();
        grid(Vphi.real(), vec);
        add(Vphi.real(), vec, 0);
        vec.clear(true);
    }
}

void QMPotential::calcImagPart(Orbital &Vphi, Orbital &phi, bool adjoint) {
    if (this->apply_prec < 0.0) MSG_ERROR("Uninitialized operator");

    MWAdder<3> add(this->apply_prec, this->max_scale);
    MWMultiplier<3> mult(this->apply_prec, this->max_scale);
    GridGenerator<3> grid(this->max_scale);

    QMPotential &V = *this;
    FunctionTreeVector<3> vec;

    if (Vphi.hasImag()) MSG_ERROR("Orbital not empty");
    if (V.hasReal() and phi.hasImag()) {
        FunctionTree<3> *tree = new FunctionTree<3>(*MRA);
        grid(*tree, phi.imag());
        mult(*tree, 1.0, V.real(), phi.imag(), 0);
        vec.push_back(1.0, tree);
    }
    if (V.hasImag() and phi.hasReal()) {
        FunctionTree<3> *tree = new FunctionTree<3>(*MRA);
        grid(*tree, phi.real());
        mult(*tree, 1.0, V.imag(), phi.real(), 0);
        if (adjoint) {
            vec.push_back(-1.0, tree);
        } else {
            vec.push_back(1.0, tree);
        }
    }
    if (vec.size() == 1) {
        Vphi.setImag(vec[0]);
        vec.clear(false);
    }
    if (vec.size() == 2) {
        Vphi.allocImag();
        grid(Vphi.imag(), vec);
        add(Vphi.imag(), vec, 0);
        vec.clear(true);
    }
}
