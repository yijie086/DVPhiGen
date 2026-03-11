#ifndef EVENT_H
#define EVENT_H

#include <string>

struct DVphiEvent {
    static constexpr double Mp = 0.9382720813;
    static constexpr double alpha = 1.0 / 137.035999084;
    static constexpr double m_pi = 3.14159265358979323846;

    double EB;
    double Q2;
    double W;
    double t;
    double xB;

    double epsilon;   // virtual photon polarization
    double gammaFlux; // Hand convention virtual photon flux
    double jacobian;    // Jacobian for transformation from (xB, Q2, t) to (W, Q2, t)

    bool hasAngles;

    double phi_e;      // scattered electron azimuth in lab [rad]
    double phi;        // hadron plane azimuth relative to lepton plane [rad]
    double cosThetaH;  // cos(theta_H)
    double phiH;       // helicity azimuth [rad]

    DVphiEvent(double EB_in, double Q2_in, double W_in, double t_in);

    void SetAngles(double phi_e_in, double phi_in, double cosThetaH_in, double phiH_in);

    static double ComputeXB(double Q2, double W);
    static double ComputeEpsilon(double EB, double Q2, double W);
    static double ComputeGammaFlux(double EB, double Q2, double W);
    static double ComputeJacobian(double xB, double Q2, double W);

    void Print(const std::string& name = "DVphiEvent") const;
};

#endif