#include "kinematics.h"
#include <iostream>
#include <iomanip>
#include <stdexcept>
#include <cmath>
#include <algorithm>

static constexpr double PI   = 3.14159265358979323846;
static constexpr double Me   = 0.000511;
static constexpr double Mp   = 0.9382720813;
static constexpr double Mphi = 1.019461;
static constexpr double MK   = 0.493677;

// ==================== Vec3 ====================

Vec3::Vec3() : x(0.0), y(0.0), z(0.0) {}
Vec3::Vec3(double x_, double y_, double z_) : x(x_), y(y_), z(z_) {}

double Vec3::Mag2() const { return x * x + y * y + z * z; }
double Vec3::Mag() const { return std::sqrt(Mag2()); }

Vec3 Vec3::Unit() const {
    double m = Mag();
    if (m <= 0.0) return Vec3(0.0, 0.0, 0.0);
    return (*this) / m;
}

Vec3 Vec3::operator+(const Vec3& other) const {
    return Vec3(x + other.x, y + other.y, z + other.z);
}

Vec3 Vec3::operator-(const Vec3& other) const {
    return Vec3(x - other.x, y - other.y, z - other.z);
}

Vec3 Vec3::operator*(double a) const {
    return Vec3(a * x, a * y, a * z);
}

Vec3 Vec3::operator/(double a) const {
    return Vec3(x / a, y / a, z / a);
}

double Vec3::Dot(const Vec3& other) const {
    return x * other.x + y * other.y + z * other.z;
}

Vec3 Vec3::Cross(const Vec3& other) const {
    return Vec3(
        y * other.z - z * other.y,
        z * other.x - x * other.z,
        x * other.y - y * other.x
    );
}

// ==================== FourVector ====================

FourVector::FourVector() : E(0.0), px(0.0), py(0.0), pz(0.0) {}

FourVector::FourVector(double E_, double px_, double py_, double pz_)
    : E(E_), px(px_), py(py_), pz(pz_) {}

double FourVector::M2() const {
    return E * E - px * px - py * py - pz * pz;
}

double FourVector::M() const {
    double m2 = M2();
    return (m2 >= 0.0) ? std::sqrt(m2) : -std::sqrt(-m2);
}

Vec3 FourVector::Vect() const {
    return Vec3(px, py, pz);
}

FourVector FourVector::operator+(const FourVector& other) const {
    return FourVector(E + other.E, px + other.px, py + other.py, pz + other.pz);
}

FourVector FourVector::operator-(const FourVector& other) const {
    return FourVector(E - other.E, px - other.px, py - other.py, pz - other.pz);
}

void FourVector::Boost(const Vec3& beta) {
    double b2 = beta.Mag2();
    if (b2 >= 1.0) {
        throw std::runtime_error("Boost beta^2 >= 1.");
    }

    double gamma = 1.0 / std::sqrt(1.0 - b2);
    double bp = beta.x * px + beta.y * py + beta.z * pz;
    double gamma2 = (b2 > 0.0) ? (gamma - 1.0) / b2 : 0.0;

    double new_px = px + gamma2 * bp * beta.x + gamma * beta.x * E;
    double new_py = py + gamma2 * bp * beta.y + gamma * beta.y * E;
    double new_pz = pz + gamma2 * bp * beta.z + gamma * beta.z * E;
    double new_E  = gamma * (E + bp);

    px = new_px;
    py = new_py;
    pz = new_pz;
    E  = new_E;
}

// ==================== FinalState ====================

void FinalState::Print() const {
    auto printP4 = [](const std::string& name, const FourVector& p) {
        std::cout << std::fixed << std::setprecision(6);
        std::cout << name
                  << " : E="  << p.E
                  << " px="   << p.px
                  << " py="   << p.py
                  << " pz="   << p.pz
                  << " M="    << p.M()
                  << "\n";
    };

    printP4("e' ", electron);
    printP4("p' ", proton);
    printP4("K+ ", kp);
    printP4("K- ", km);
    std::cout << std::endl;
}

// ==================== KinematicsCalculator ====================

double KinematicsCalculator::Lambda(double x, double y, double z) {
    return x * x + y * y + z * z - 2.0 * x * y - 2.0 * x * z - 2.0 * y * z;
}

FinalState KinematicsCalculator::BuildEvent(const DVphiEvent& evt) {
    if (!evt.hasAngles) {
        throw std::runtime_error("Event angles are not set.");
    }

    // ---------------- initial state in lab ----------------
    double pz_beam = std::sqrt(std::max(0.0, evt.EB * evt.EB - Me * Me));
    FourVector k_in(evt.EB, 0.0, 0.0, pz_beam);
    FourVector p_in(Mp, 0.0, 0.0, 0.0);

    // ---------------- scattered electron from EB,Q2,W ----------------
    double nu = (evt.W * evt.W + evt.Q2 - Mp * Mp) / (2.0 * Mp);
    double Eep = evt.EB - nu;
    if (Eep <= Me) {
        throw std::runtime_error("Scattered electron energy <= me.");
    }
    
    double cosThetaE = 1.0 - evt.Q2 / (2.0 * evt.EB * Eep);
    if (cosThetaE < -1.000001 || cosThetaE > 1.000001) {
        throw std::runtime_error("Unphysical combination in BuildEvent.");
    }
    cosThetaE = std::clamp(cosThetaE, -1.0, 1.0);

    double sinThetaE = std::sqrt(std::max(0.0, 1.0 - cosThetaE * cosThetaE));
    double pep = std::sqrt(std::max(0.0, Eep * Eep - Me * Me));

    // scattered electron with random lab azimuth phi_e
    FourVector k_out(
        Eep,
        pep * sinThetaE * std::cos(evt.phi_e),
        pep * sinThetaE * std::sin(evt.phi_e),
        pep * cosThetaE
    );

    FourVector q = k_in - k_out;

    // ---------------- gamma* p CM ----------------
    FourVector total_gp = q + p_in;
    Vec3 beta_cm = total_gp.Vect() / total_gp.E;

    FourVector q_cm = q;
    FourVector p_cm = p_in;
    FourVector k_in_cm = k_in;

    q_cm.Boost(beta_cm * (-1.0));
    p_cm.Boost(beta_cm * (-1.0));
    k_in_cm.Boost(beta_cm * (-1.0));

    double s = evt.W * evt.W;

    double pq_cm = std::sqrt(std::max(0.0, Lambda(s, -evt.Q2, Mp * Mp))) / (2.0 * evt.W);
    double pphi_cm = std::sqrt(std::max(0.0, Lambda(s, Mphi * Mphi, Mp * Mp))) / (2.0 * evt.W);

    double Eq_cm   = (s - Mp * Mp - evt.Q2) / (2.0 * evt.W);
    double Ephi_cm = (s + Mphi * Mphi - Mp * Mp) / (2.0 * evt.W);
    double Ep_cm   = (s + Mp * Mp - Mphi * Mphi) / (2.0 * evt.W);

    // t = Mphi^2 - Q^2 - 2 Eq Ephi + 2 pq pphi cos(theta*)
    double numerator = evt.t - Mphi * Mphi + evt.Q2 + 2.0 * Eq_cm * Ephi_cm;
    double denominator = 2.0 * pq_cm * pphi_cm;
    if (std::abs(denominator) < 1e-12) {
        throw std::runtime_error("Denominator too small in cos(theta*) calculation.");
    }

    double cosThetaStar = numerator / denominator;
    cosThetaStar = std::clamp(cosThetaStar, -1.0, 1.0);
    double sinThetaStar = std::sqrt(std::max(0.0, 1.0 - cosThetaStar * cosThetaStar));

    // ---------------- hadron production basis in gamma*p CM ----------------
    // z-axis along q
    Vec3 zAxis = q_cm.Vect().Unit();
    if (zAxis.Mag() < 1e-12) {
        throw std::runtime_error("q_cm direction is ill-defined.");
    }

    // x-axis lies in the lepton plane
    Vec3 kinVec = k_in_cm.Vect().Unit();

    // remove z component to get transverse direction in the lepton plane
    Vec3 xAxis = kinVec - zAxis * kinVec.Dot(zAxis);
    if (xAxis.Mag() < 1e-10) {
        xAxis = Vec3(1.0, 0.0, 0.0);
    } else {
        xAxis = xAxis.Unit();
    }

    Vec3 yAxis = zAxis.Cross(xAxis);
    if (yAxis.Mag() < 1e-10) {
        yAxis = Vec3(0.0, 1.0, 0.0);
    } else {
        yAxis = yAxis.Unit();
    }

    // phi meson direction in gamma*p CM
    Vec3 phi_dir =
        xAxis * (sinThetaStar * std::cos(evt.phi)) +
        yAxis * (sinThetaStar * std::sin(evt.phi)) +
        zAxis * cosThetaStar;

    phi_dir = phi_dir.Unit();

    FourVector phi_cm(
        Ephi_cm,
        pphi_cm * phi_dir.x,
        pphi_cm * phi_dir.y,
        pphi_cm * phi_dir.z
    );

    FourVector proton_cm(
        Ep_cm,
        -phi_cm.px,
        -phi_cm.py,
        -phi_cm.pz
    );

    // ---------------- phi -> K+ K- in phi rest ----------------
    double EK_star = Mphi / 2.0;
    double pK_star = std::sqrt(std::max(0.0, EK_star * EK_star - MK * MK));
    double sinThetaH = std::sqrt(std::max(0.0, 1.0 - evt.cosThetaH * evt.cosThetaH));

    // helicity frame basis
    Vec3 zH = phi_cm.Vect().Unit();
    if (zH.Mag() < 1e-12) {
        throw std::runtime_error("phi_cm direction is ill-defined.");
    }

    Vec3 zQ = q_cm.Vect().Unit();

    Vec3 yH = zQ.Cross(zH);
    if (yH.Mag() < 1e-10) {
        yH = Vec3(0.0, 1.0, 0.0);
    } else {
        yH = yH.Unit();
    }

    Vec3 xH = yH.Cross(zH);
    if (xH.Mag() < 1e-10) {
        xH = Vec3(1.0, 0.0, 0.0);
    } else {
        xH = xH.Unit();
    }

    Vec3 kplus_dir =
        xH * (sinThetaH * std::cos(evt.phiH)) +
        yH * (sinThetaH * std::sin(evt.phiH)) +
        zH * evt.cosThetaH;

    kplus_dir = kplus_dir.Unit();

    FourVector kp_phiRest(
        EK_star,
        pK_star * kplus_dir.x,
        pK_star * kplus_dir.y,
        pK_star * kplus_dir.z
    );

    FourVector km_phiRest(
        EK_star,
        -kp_phiRest.px,
        -kp_phiRest.py,
        -kp_phiRest.pz
    );

    // ---------------- boost kaons from phi rest -> gamma*p CM ----------------
    Vec3 beta_phi_cm = phi_cm.Vect() / phi_cm.E;

    FourVector kp_cm = kp_phiRest;
    FourVector km_cm = km_phiRest;

    kp_cm.Boost(beta_phi_cm);
    km_cm.Boost(beta_phi_cm);

    // ---------------- boost hadrons from gamma*p CM -> lab ----------------
    FourVector proton_lab = proton_cm;
    FourVector kp_lab = kp_cm;
    FourVector km_lab = km_cm;

    proton_lab.Boost(beta_cm);
    kp_lab.Boost(beta_cm);
    km_lab.Boost(beta_cm);

    FinalState fs;
    fs.electron = k_out;
    fs.proton   = proton_lab;
    fs.kp       = kp_lab;
    fs.km       = km_lab;

    return fs;
}