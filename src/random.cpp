#include "random.h"
#include <cmath>
#include <algorithm>
#include <stdexcept>
#include "event.h"

static constexpr double PI   = 3.14159265358979323846;
static constexpr double Mp   = 0.9382720813;
static constexpr double Mphi = 1.019461;

RandomGenerator::RandomGenerator(unsigned int seed)
    : state(seed),
      Q2_min(1.0),
      Q2_max(8.0),
      W_min(2.0),
      W_max(4.0),
      wmax(1.0),
      EB(10.6)
{
    ValidateRanges();
    UpdateWmax(EB);
}

unsigned int RandomGenerator::NextUInt() {
    state = 1664525u * state + 1013904223u;
    return state;
}

double RandomGenerator::Uniform(double a, double b) {
    double x = static_cast<double>(NextUInt()) / 4294967295.0;
    return a + (b - a) * x;
}

void RandomGenerator::SetRanges(double Q2_min_in, double Q2_max_in,
                                double W_min_in,  double W_max_in, double EB_in) {
    Q2_min = Q2_min_in;
    Q2_max = Q2_max_in;
    W_min  = W_min_in;
    W_max  = W_max_in;
    EB = EB_in;

    ValidateRanges();
    UpdateWmax(EB);
}

void RandomGenerator::ValidateRanges() const {
    if (Q2_min <= 0.0) {
        throw std::runtime_error("Q2_min must be > 0.");
    }
    if (Q2_max <= Q2_min) {
        throw std::runtime_error("Q2_max must be > Q2_min.");
    }
    if (W_min <= 0.0) {
        throw std::runtime_error("W_min must be > 0.");
    }
    if (W_max <= W_min) {
        throw std::runtime_error("W_max must be > W_min.");
    }

    // also require upper corner to be physically allowed for phi production
    double s_test = W_max * W_max;
    double lambda_out = Lambda(s_test, Mphi * Mphi, Mp * Mp);
    if (lambda_out < 0.0) {
        throw std::runtime_error("W range does not allow phi production.");
    }
}

double RandomGenerator::GenerateQ2() {
    return Uniform(Q2_min, Q2_max);
}

double RandomGenerator::GenerateW() {
    return Uniform(W_min, W_max);
}

double RandomGenerator::GeneratePhiE() {
    return Uniform(0.0, 2.0 * PI);
}

double RandomGenerator::GeneratePhi() {
    return Uniform(0.0, 2.0 * PI);
}

double RandomGenerator::WeightCosThetaH(double EB, double Q2, double W, double t, double cosThetaH) const { 
    (void)t;
    (void)cosThetaH;

    double w = 1.0;
    double epsilon = DVphiEvent::ComputeEpsilon(EB, Q2, W);
    double cR = 0.4;
    double R = cR*Q2 / (Mphi*Mphi);
    double r0004 = epsilon * R / (1.0 + epsilon * R);
    w = (3.0/4.0)*((1-r0004) + (3.0*r0004-1.0)*cosThetaH*cosThetaH);
    if (!std::isfinite(w) || w < 0.0) return 0.0;
    return w;
}

/*
// Original weighted cos(theta_H) generator. Restore this block to recover
// the helicity-angle distribution from WeightCosThetaH().
double RandomGenerator::GenerateCosThetaH(double EB, double Q2, double W, double t) {
    double w0   = WeightCosThetaH(EB, Q2, W, t, 0.0);
    double w1   = WeightCosThetaH(EB, Q2, W, t, 1.0);
    double wmax = 2 * std::max(w0, w1);

    if (wmax <= 0.0 || !std::isfinite(wmax)) {
        throw std::runtime_error("Invalid max weight for cosThetaH.");
    }

    while (true) {
        double cosThetaH = Uniform(-1.0, 1.0);
        double w = WeightCosThetaH(EB, Q2, W, t, cosThetaH);

        if (w > wmax) {
            throw std::runtime_error("WeightCosThetaH exceeds estimated maximum.");
        }

        if (Uniform(0.0, wmax) < w) {
            return cosThetaH;
        }
    }
}
*/

double RandomGenerator::GenerateCosThetaH(double EB, double Q2, double W, double t) {
    (void)EB;
    (void)Q2;
    (void)W;
    (void)t;
    return Uniform(-1.0, 1.0);
}

double RandomGenerator::GeneratePhiH() {
    return Uniform(0.0, 2.0 * PI);
}

double RandomGenerator::Lambda(double x, double y, double z) {
    return x*x + y*y + z*z - 2.0*x*y - 2.0*x*z - 2.0*y*z;
}

void RandomGenerator::GetTRange(double Q2, double W, double& t_low, double& t_high) const {
    double s = W * W;

    double lambda_in  = Lambda(s, -Q2, Mp * Mp);
    double lambda_out = Lambda(s, Mphi * Mphi, Mp * Mp);

    if (lambda_in < 0.0 || lambda_out < 0.0) {
        throw std::runtime_error("Unphysical Q2/W combination in GetTRange.");
    }

    double pq   = std::sqrt(lambda_in)  / (2.0 * W);
    double pphi = std::sqrt(lambda_out) / (2.0 * W);

    double Eq   = (s - Mp * Mp - Q2) / (2.0 * W);
    double Ephi = (s + Mphi * Mphi - Mp * Mp) / (2.0 * W);

    // cos(theta*) = +1
    double t1 = Mphi * Mphi - Q2 - 2.0 * Eq * Ephi + 2.0 * pq * pphi;

    // cos(theta*) = -1
    double t2 = Mphi * Mphi - Q2 - 2.0 * Eq * Ephi - 2.0 * pq * pphi;

    t_low  = std::min(t1, t2);
    t_high = std::max(t1, t2);
}

double RandomGenerator::GenerateT(double Q2, double W) {
    double t_low, t_high;
    GetTRange(Q2, W, t_low, t_high);
    return Uniform(t_low, t_high);
}

double RandomGenerator::WeightQ2Wt(double EB, double Q2, double W, double t) const {
    double Wth= Mp + Mphi; // threshold W for phi production
    double nuT = 3.0;
    double alpha1 = 400;
    double alpha2 = 1.0;
    double alpha3 = 0.32;
    double cR = 0.4;
    double cT = alpha1 * std::pow( 1 - (Wth*Wth)/(W*W), alpha2) * std::pow(W, alpha3);
    double sigmaT = cT / std::pow((1+Q2/(Mphi*Mphi)), nuT);
    double R = cR*Q2 / (Mphi*Mphi);
    double sigmaL = R * sigmaT;
    double epsilon = DVphiEvent::ComputeEpsilon(EB, Q2, W);
    double gammaFlux = DVphiEvent::ComputeGammaFlux(EB, Q2, W);
    double sigma = gammaFlux * (sigmaT + epsilon * sigmaL);
    double jacobian = DVphiEvent::ComputeJacobian(DVphiEvent::ComputeXB(Q2, W), Q2, W);
    sigma *= jacobian;

    double B0 = 2.2;
    double alphaprime = 0.24;
    double slope_b = B0 + 4.0 * alphaprime * std::log(W);

    sigma = sigma * std::exp(slope_b * t);

    // current target distribution:
    //   f(Q2,W,t) ~ exp(b t)
    // since t < 0, this suppresses large |t|
    double w = sigma;

    if (!std::isfinite(w) || w < 0.0) return 0.0;
    return w;
}

void RandomGenerator::UpdateWmax(double EB) {
    // According to your rule:
    // wmax = 2 * weight(minQ2, maxW, t_min)
    //
    // Here "t_min" is interpreted as the physical t closest to 0
    // at (Q2_min, W_max), i.e. t_high after sorting.

    double t_low, t_high;
    GetTRange(Q2_min, W_max, t_low, t_high);

    double t_min_near_zero = t_high;
    double w_ref = WeightQ2Wt(EB, Q2_min, W_max, t_min_near_zero);

    if (w_ref <= 0.0 || !std::isfinite(w_ref)) {
        throw std::runtime_error("Invalid reference weight in UpdateWmax.");
    }

    wmax = 2.0 * w_ref;
}

/*
// Original weighted Q2/W/t generator. Restore this block to recover
// accept-reject sampling with WeightQ2Wt().
void RandomGenerator::GenerateKinematics(double EB, double& Q2, double& W, double& t) {
    while (true) {
        double q2_try = GenerateQ2();
        double w_try  = GenerateW();
        double t_try  = GenerateT(q2_try, w_try);

        double weight = WeightQ2Wt(EB, q2_try, w_try, t_try);

        if (weight > wmax) {
            throw std::runtime_error("Weight exceeds wmax. Increase safety factor or check model.");
        }

        double r = Uniform(0.0, wmax);
        if (r < weight) {
            Q2 = q2_try;
            W  = w_try;
            t  = t_try;
            return;
        }
    }
}
*/

void RandomGenerator::GenerateKinematics(double EB, double& Q2, double& W, double& t) {
    (void)EB;

    Q2 = GenerateQ2();
    W  = GenerateW();
    t  = GenerateT(Q2, W);
}
