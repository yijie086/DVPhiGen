#include "kinematics.h"

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

static constexpr double Me = 0.000511;
static constexpr double Mp = 0.9382720813;
static constexpr double PI = 3.14159265358979323846;

struct Particle {
    int index = 0;
    int charge = 0;
    int type = 0;
    int pid = 0;
    int parent = 0;
    int daughter = 0;
    FourVector p4;
    double mass = 0.0;
};

static double NormalizePhi(double phi) {
    while (phi < 0.0) phi += 2.0 * PI;
    while (phi >= 2.0 * PI) phi -= 2.0 * PI;
    return phi;
}

static double Clamp(double x, double lo, double hi) {
    return std::max(lo, std::min(hi, x));
}

static double Lambda(double x, double y, double z) {
    return x * x + y * y + z * z - 2.0 * x * y - 2.0 * x * z - 2.0 * y * z;
}

static void GetTRange(double Q2, double W, double Mphi, double& tLow, double& tHigh) {
    double s = W * W;
    double lambdaIn = Lambda(s, -Q2, Mp * Mp);
    double lambdaOut = Lambda(s, Mphi * Mphi, Mp * Mp);

    if (lambdaIn < 0.0 || lambdaOut < 0.0) {
        throw std::runtime_error("Unphysical Q2/W combination while computing t range.");
    }

    double pq = std::sqrt(lambdaIn) / (2.0 * W);
    double pphi = std::sqrt(lambdaOut) / (2.0 * W);
    double Eq = (s - Mp * Mp - Q2) / (2.0 * W);
    double Ephi = (s + Mphi * Mphi - Mp * Mp) / (2.0 * W);

    double t1 = Mphi * Mphi - Q2 - 2.0 * Eq * Ephi + 2.0 * pq * pphi;
    double t2 = Mphi * Mphi - Q2 - 2.0 * Eq * Ephi - 2.0 * pq * pphi;
    tLow = std::min(t1, t2);
    tHigh = std::max(t1, t2);
}

static Particle ParseParticleLine(const std::string& line) {
    Particle p;
    double px = 0.0;
    double py = 0.0;
    double pz = 0.0;
    double e = 0.0;
    std::istringstream iss(line);
    iss >> p.index >> p.charge >> p.type >> p.pid >> p.parent >> p.daughter
        >> px >> py >> pz >> e >> p.mass;
    if (!iss) {
        throw std::runtime_error("Could not parse LUND particle line: " + line);
    }
    p.p4 = FourVector(e, px, py, pz);
    return p;
}

static double Momentum(const FourVector& p) {
    return std::sqrt(p.px * p.px + p.py * p.py + p.pz * p.pz);
}

static double ThetaDeg(const FourVector& p) {
    double mag = Momentum(p);
    if (mag <= 0.0) return 0.0;
    return std::acos(Clamp(p.pz / mag, -1.0, 1.0)) * 180.0 / PI;
}

static double PhiDeg(const FourVector& p) {
    return NormalizePhi(std::atan2(p.py, p.px)) * 180.0 / PI;
}

static Vec3 UnitOr(const Vec3& v, const Vec3& fallback) {
    if (v.Mag() < 1e-12) return fallback;
    return v.Unit();
}

int main(int argc, char* argv[]) {
    if (argc < 4) {
        std::cerr << "Usage: " << argv[0]
                  << " input.lund particle.csv event.csv [beam_energy]\n";
        return 2;
    }

    const std::string inputPath = argv[1];
    const std::string particleCsvPath = argv[2];
    const std::string eventCsvPath = argv[3];
    const double fallbackEB = (argc >= 5) ? std::stod(argv[4]) : 10.6;

    std::ifstream fin(inputPath);
    if (!fin.is_open()) {
        throw std::runtime_error("Could not open input LUND file: " + inputPath);
    }

    std::ofstream particleCsv(particleCsvPath);
    std::ofstream eventCsv(eventCsvPath);
    particleCsv << "event,index,pid,p,theta_deg,phi_deg,px,py,pz,E,mass\n";
    eventCsv << "event,EB,xB,Q2,W,t,minus_t,t_forward,tprime,phi,cosThetaH,phiH,epsilon,Mkk\n";

    std::string header;
    int eventNumber = 0;
    while (std::getline(fin, header)) {
        if (header.empty()) continue;

        std::istringstream hss(header);
        std::vector<double> hvals;
        double hv = 0.0;
        while (hss >> hv) hvals.push_back(hv);
        if (hvals.empty()) continue;

        int nParticles = static_cast<int>(hvals[0]);
        double EB = fallbackEB;
        if (hvals.size() > 6 && hvals[6] > 0.0) EB = hvals[6];

        std::vector<Particle> particles;
        particles.reserve(nParticles);
        for (int i = 0; i < nParticles; ++i) {
            std::string line;
            if (!std::getline(fin, line)) {
                throw std::runtime_error("Unexpected EOF inside event.");
            }
            particles.push_back(ParseParticleLine(line));
        }

        ++eventNumber;
        std::map<int, FourVector> byPid;
        for (const Particle& p : particles) {
            byPid[p.pid] = p.p4;
            particleCsv << eventNumber << ","
                        << p.index << ","
                        << p.pid << ","
                        << Momentum(p.p4) << ","
                        << ThetaDeg(p.p4) << ","
                        << PhiDeg(p.p4) << ","
                        << p.p4.px << ","
                        << p.p4.py << ","
                        << p.p4.pz << ","
                        << p.p4.E << ","
                        << p.mass << "\n";
        }

        if (!byPid.count(11) || !byPid.count(2212) ||
            !byPid.count(321) || !byPid.count(-321)) {
            continue;
        }

        FourVector kOut = byPid[11];
        FourVector pOut = byPid[2212];
        (void)pOut;
        FourVector kp = byPid[321];
        FourVector km = byPid[-321];

        double pzBeam = std::sqrt(std::max(0.0, EB * EB - Me * Me));
        FourVector kIn(EB, 0.0, 0.0, pzBeam);
        FourVector pIn(Mp, 0.0, 0.0, 0.0);
        FourVector q = kIn - kOut;
        FourVector gp = q + pIn;
        FourVector phiMeson = kp + km;

        double Q2 = -q.M2();
        double W = gp.M();
        double t = (q - phiMeson).M2();
        double minusT = -t;
        double xB = DVphiEvent::ComputeXB(Q2, W);
        double epsilon = DVphiEvent::ComputeEpsilon(EB, Q2, W);
        double Mkk = phiMeson.M();
        double tLow = 0.0;
        double tHigh = 0.0;
        GetTRange(Q2, W, Mkk, tLow, tHigh);
        double tForward = tHigh;
        double tprime = tForward - t;

        Vec3 betaCm = gp.Vect() / gp.E;
        FourVector qCm = q;
        FourVector kInCm = kIn;
        FourVector phiCm = phiMeson;
        FourVector kpCm = kp;

        qCm.Boost(betaCm * (-1.0));
        kInCm.Boost(betaCm * (-1.0));
        phiCm.Boost(betaCm * (-1.0));
        kpCm.Boost(betaCm * (-1.0));

        Vec3 zAxis = UnitOr(qCm.Vect(), Vec3(0.0, 0.0, 1.0));
        Vec3 kinVec = UnitOr(kInCm.Vect(), Vec3(0.0, 0.0, 1.0));
        Vec3 xAxis = kinVec - zAxis * kinVec.Dot(zAxis);
        xAxis = UnitOr(xAxis, Vec3(1.0, 0.0, 0.0));
        Vec3 yAxis = UnitOr(zAxis.Cross(xAxis), Vec3(0.0, 1.0, 0.0));

        Vec3 phiDir = UnitOr(phiCm.Vect(), Vec3(0.0, 0.0, 1.0));
        double eventPhi = NormalizePhi(std::atan2(phiDir.Dot(yAxis), phiDir.Dot(xAxis)));

        Vec3 zH = phiDir;
        Vec3 zQ = zAxis;
        Vec3 yH = UnitOr(zQ.Cross(zH), Vec3(0.0, 1.0, 0.0));
        Vec3 xH = UnitOr(yH.Cross(zH), Vec3(1.0, 0.0, 0.0));

        Vec3 betaPhiCm = phiCm.Vect() / phiCm.E;
        FourVector kpPhiRest = kpCm;
        kpPhiRest.Boost(betaPhiCm * (-1.0));
        Vec3 kpDir = UnitOr(kpPhiRest.Vect(), Vec3(0.0, 0.0, 1.0));

        double cosThetaH = Clamp(kpDir.Dot(zH), -1.0, 1.0);
        double phiH = NormalizePhi(std::atan2(kpDir.Dot(yH), kpDir.Dot(xH)));

        eventCsv << eventNumber << ","
                 << EB << ","
                 << xB << ","
                 << Q2 << ","
                 << W << ","
                 << t << ","
                 << minusT << ","
                 << tForward << ","
                 << tprime << ","
                 << eventPhi << ","
                 << cosThetaH << ","
                 << phiH << ","
                 << epsilon << ","
                 << Mkk << "\n";
    }

    return 0;
}
