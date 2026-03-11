#include "event.h"
#include "random.h"
#include "kinematics.h"

#include <iostream>
#include <fstream>
#include <iomanip>
#include <stdexcept>
#include <string>
#include <cmath>

// ---------------- masses ----------------
static constexpr double Me = 0.000511;
static constexpr double Mp = 0.9382720813;
static constexpr double MK = 0.493677;

// ---------------- write one particle line ----------------
// A simple LUND-style particle record:
// index charge type pid parent daughter px py pz E mass vx vy vz
void WriteLundParticle(std::ofstream& fout,
                       int index,
                       int charge,
                       int type,
                       int pid,
                       int parent,
                       int daughter,
                       const FourVector& p,
                       double mass,
                       double vx = 0.0,
                       double vy = 0.0,
                       double vz = 0.0) {
    fout << std::setw(4)  << index
         << std::setw(4)  << charge
         << std::setw(4)  << type
         << std::setw(8)  << pid
         << std::setw(4)  << parent
         << std::setw(4)  << daughter
         << std::setw(16) << std::setprecision(8) << std::fixed << p.px
         << std::setw(16) << std::setprecision(8) << std::fixed << p.py
         << std::setw(16) << std::setprecision(8) << std::fixed << p.pz
         << std::setw(16) << std::setprecision(8) << std::fixed << p.E
         << std::setw(16) << std::setprecision(8) << std::fixed << mass
         << std::setw(12) << std::setprecision(4) << std::fixed << vx
         << std::setw(12) << std::setprecision(4) << std::fixed << vy
         << std::setw(12) << std::setprecision(4) << std::fixed << vz
         << "\n";
}

// ---------------- write one event ----------------
// First line is a simple event header.
// Then 4 particle lines: e', p', K+, K-
void WriteLundEvent(std::ofstream& fout,
                    const DVphiEvent& evt,
                    const FinalState& fs) {
    // Simple header:
    // NPART, eventID, beamE, Q2, W, xB, t, phi, cosThetaH, phiH
    fout << std::setw(4)  << 4
        << std::setw(4)  << 1
        << std::setw(4)  << 1
        << std::setw(4)  << 0
        << std::setw(4)  << 0
        << std::setw(6)  << 11
        << std::setw(14) << std::fixed << std::setprecision(6) << evt.EB
        << std::setw(8)  << 2212
        << std::setw(6)  << 0
        << std::setw(14) << std::fixed << std::setprecision(6) << 1.0
        << "\n";

    // particle lines
    // electron: pid=11, charge=-1
    WriteLundParticle(fout, 1, -1, 1,   11, 0, 0, fs.electron, Me);

    // proton: pid=2212, charge=+1
    WriteLundParticle(fout, 2, +1, 1, 2212, 0, 0, fs.proton,   Mp);

    // K+: pid=321, charge=+1
    WriteLundParticle(fout, 3, +1, 1,  321, 0, 0, fs.kp,       MK);

    // K-: pid=-321, charge=-1
    WriteLundParticle(fout, 4, -1, 1, -321, 0, 0, fs.km,       MK);
}

int main() {
    try {
        const double EB = 10.6;
        const int NEvents = 100000;

        RandomGenerator rng(20260310);

        rng.SetRanges(1.0, 8.0,   // Q2_min, Q2_max
              2.0, 4.0,  // W_min,  W_max
              EB);

        std::ofstream fout("dvphi.lund");
        if (!fout.is_open()) {
            throw std::runtime_error("Cannot open output file dvphi.lund");
        }

        int nAccepted = 0;
        int iEvent = 0;

        while (nAccepted < NEvents) {
            ++iEvent;

            try {
                double Q2, W, t;
                rng.GenerateKinematics(EB, Q2, W, t);
                double phi_e = rng.GeneratePhiE();
                double phi  = rng.GeneratePhi();
                double cth  = rng.GenerateCosThetaH(EB, Q2, W, t);
                double phiH = rng.GeneratePhiH();

                DVphiEvent evt(EB, Q2, W, t);
                evt.SetAngles(phi_e, phi, cth, phiH);

                FinalState fs = KinematicsCalculator::BuildEvent(evt);

                WriteLundEvent(fout, evt, fs);
                ++nAccepted;
            }
            catch (const std::exception&) {
                // skip unphysical event and continue
                continue;
            }
        }

        fout.close();

        std::cout << "Generated " << nAccepted
                  << " events into dvphi.lund\n";
    }
    catch (const std::exception& ex) {
        std::cerr << "Error: " << ex.what() << std::endl;
        return 1;
    }

    return 0;
}