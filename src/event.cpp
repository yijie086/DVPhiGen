#include "event.h"
#include <iostream>
#include <iomanip>
#include <stdexcept>

DVphiEvent::DVphiEvent(double EB_in, double Q2_in, double W_in, double t_in)
    : EB(EB_in),
      Q2(Q2_in),
      W(W_in),
      t(t_in),
      xB(ComputeXB(Q2_in, W_in)),
      epsilon(ComputeEpsilon(EB_in, Q2_in, W_in)),
      gammaFlux(ComputeGammaFlux(EB_in, Q2_in, W_in)),
      jacobian(ComputeJacobian(ComputeXB(Q2_in, W_in), Q2_in, W_in)),
      hasAngles(false),
      phi_e(0.0),
      phi(0.0),
      cosThetaH(0.0),
      phiH(0.0)
{
    if (EB <= 0.0) throw std::runtime_error("EB must be positive.");
    if (Q2 <= 0.0) throw std::runtime_error("Q2 must be positive.");
    if (W <= 0.0)  throw std::runtime_error("W must be positive.");
    if (xB <= 0.0 || xB >= 1.0) throw std::runtime_error("Computed xB out of range.");
}

void DVphiEvent::SetAngles(double phi_e_in, double phi_in, double cosThetaH_in, double phiH_in) {
    if (cosThetaH_in < -1.0 || cosThetaH_in > 1.0) {
        throw std::runtime_error("cosThetaH must be in [-1,1].");
    }
    phi_e = phi_e_in;
    phi = phi_in;
    cosThetaH = cosThetaH_in;
    phiH = phiH_in;
    hasAngles = true;
}

double DVphiEvent::ComputeXB(double Q2, double W) {
    const double denom = W * W - Mp * Mp + Q2;
    if (denom <= 0.0) {
        throw std::runtime_error("Invalid denominator for xB.");
    }
    return Q2 / denom;
}

double DVphiEvent::ComputeEpsilon(double EB, double Q2, double W) {
    double y = (W * W + Q2 - Mp * Mp) / (2.0 * Mp * EB);
    if (y <= 0.0 || y >= 1.0) {
        throw std::runtime_error("Invalid y for epsilon calculation.");
    }
    double epsilon_num = 1.0 - y - (Q2 / (4.0 * EB * EB));
    double epsilon_den = 1.0 - y + (y * y / 2.0) + (Q2 / (4.0 * EB * EB));
    if (epsilon_den <= 0.0) {
        throw std::runtime_error("Invalid denominator for epsilon.");
    }
    return epsilon_num / epsilon_den;
}

double DVphiEvent::ComputeGammaFlux(double EB, double Q2, double W) {
    double nu = (W * W + Q2 - Mp * Mp) / (2.0 * Mp);
    double xB = ComputeXB(Q2, W);
    double epsilon = ComputeEpsilon(EB, Q2, W);
    if (nu <= 0.0) {
        throw std::runtime_error("Invalid nu for gamma flux.");
    }
    return (alpha / (8.0 * m_pi )) * (Q2 / EB*EB*Mp*Mp) * ((1-xB)/xB*xB*xB) * (1.0 / (1.0 - epsilon));
}

double DVphiEvent::ComputeJacobian(double xB, double Q2, double W) {
    double dW_dxb = Q2 / (2.0*W*xB*xB);
    return std::abs(dW_dxb);
}

void DVphiEvent::Print(const std::string& name) const {
    std::cout << "==== " << name << " ====\n";
    std::cout << std::fixed << std::setprecision(6);
    std::cout << "EB         = " << EB << " GeV\n";
    std::cout << "Q2         = " << Q2 << " GeV^2\n";
    std::cout << "W          = " << W  << " GeV\n";
    std::cout << "t          = " << t  << " GeV^2\n";
    std::cout << "xB         = " << xB << "\n";
    std::cout << "epsilon     = " << epsilon << "\n";
    std::cout << "gammaFlux   = " << gammaFlux << " GeV^-2\n";

    if (hasAngles) {
        std::cout << "phi_e      = " << phi_e << " rad\n";
        std::cout << "phi        = " << phi << " rad\n";
        std::cout << "cosThetaH  = " << cosThetaH << "\n";
        std::cout << "phiH       = " << phiH << " rad\n";
    } else {
        std::cout << "angles      = [not set]\n";
    }
    std::cout << std::endl;
}