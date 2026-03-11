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
    WriteLundParticle(fout, 1, -1, 1,   11, 0, 0, fs.electron, Me);
    WriteLundParticle(fout, 2, +1, 1, 2212, 0, 0, fs.proton,   Mp);
    WriteLundParticle(fout, 3, +1, 1,  321, 0, 0, fs.kp,       MK);
    WriteLundParticle(fout, 4, -1, 1, -321, 0, 0, fs.km,       MK);
}

void PrintHelp(const char* progName) {
    std::cout
        << "Usage:\n"
        << "  " << progName << " [EB] [Q2min] [Q2max] [Wmin] [Wmax] [NEvents] [outfile]\n\n"
        << "Arguments (all optional):\n"
        << "  EB       Beam energy in GeV        (default: 10.6)\n"
        << "  Q2min    Minimum Q^2 in GeV^2      (default: 1.0)\n"
        << "  Q2max    Maximum Q^2 in GeV^2      (default: 8.0)\n"
        << "  Wmin     Minimum W in GeV          (default: 2.0)\n"
        << "  Wmax     Maximum W in GeV          (default: 4.0)\n"
        << "  NEvents  Number of accepted events (default: 10000000)\n"
        << "  outfile  Output LUND file name     (default: dvphi.lund)\n\n"
        << "Examples:\n"
        << "  " << progName << "\n"
        << "  " << progName << " 10.6\n"
        << "  " << progName << " 10.6 1.0 8.0 2.0 4.0 500000 dvphi.lund\n"
        << "  " << progName << " --help\n\n"
        << "Notes:\n"
        << "  - If an argument is omitted, the default value is used.\n"
        << "  - Progress is printed every 10% of accepted events.\n";
}

int main(int argc, char* argv[]) {
    try {
        // ---------- defaults ----------
        double EB     = 10.6;
        double Q2min  = 1.0;
        double Q2max  = 8.0;
        double Wmin   = 2.0;
        double Wmax   = 4.0;
        int NEvents   = 100000;
        std::string outFile = "dvphi.lund";

        // ---------- help ----------
        if (argc > 1) {
            std::string arg1 = argv[1];
            if (arg1 == "-h" || arg1 == "--help") {
                PrintHelp(argv[0]);
                return 0;
            }
        }

        // ---------- positional args ----------
        if (argc > 1) EB     = std::stod(argv[1]);
        if (argc > 2) Q2min  = std::stod(argv[2]);
        if (argc > 3) Q2max  = std::stod(argv[3]);
        if (argc > 4) Wmin   = std::stod(argv[4]);
        if (argc > 5) Wmax   = std::stod(argv[5]);
        if (argc > 6) NEvents = std::stoi(argv[6]);
        if (argc > 7) outFile = argv[7];

        if (argc > 8) {
            std::cerr << "Error: too many arguments.\n\n";
            PrintHelp(argv[0]);
            return 1;
        }

        // ---------- sanity checks ----------
        if (EB <= 0.0) {
            throw std::runtime_error("EB must be > 0.");
        }
        if (Q2min <= 0.0 || Q2max <= 0.0 || Q2min >= Q2max) {
            throw std::runtime_error("Require 0 < Q2min < Q2max.");
        }
        if (Wmin <= 0.0 || Wmax <= 0.0 || Wmin >= Wmax) {
            throw std::runtime_error("Require 0 < Wmin < Wmax.");
        }
        if (NEvents <= 0) {
            throw std::runtime_error("NEvents must be > 0.");
        }

        // ---------- print config ----------
        std::cout << "========== Generator configuration ==========\n";
        std::cout << "EB      = " << EB << " GeV\n";
        std::cout << "Q2 range= [" << Q2min << ", " << Q2max << "] GeV^2\n";
        std::cout << "W  range= [" << Wmin  << ", " << Wmax  << "] GeV\n";
        std::cout << "NEvents = " << NEvents << "\n";
        std::cout << "Output  = " << outFile << "\n";
        std::cout << "============================================\n";

        RandomGenerator rng(20260310);
        rng.SetRanges(Q2min, Q2max, Wmin, Wmax, EB);

        std::ofstream fout(outFile);
        if (!fout.is_open()) {
            throw std::runtime_error("Cannot open output file dvphi.lund");
        }

        int nAccepted = 0;
        int nTried = 0;

        // progress every 10%
        int reportStep = NEvents / 10;
        if (reportStep < 1) reportStep = 1;
        int nextReport = reportStep;

        while (nAccepted < NEvents) {
            ++nTried;

            try {
                double Q2, W, t;
                rng.GenerateKinematics(EB, Q2, W, t);

                double phi_e = rng.GeneratePhiE();
                double phi   = rng.GeneratePhi();
                double cth   = rng.GenerateCosThetaH(EB, Q2, W, t);
                double phiH  = rng.GeneratePhiH();

                DVphiEvent evt(EB, Q2, W, t);
                evt.SetAngles(phi_e, phi, cth, phiH);

                FinalState fs = KinematicsCalculator::BuildEvent(evt);

                WriteLundEvent(fout, evt, fs);
                ++nAccepted;

                if (nAccepted >= nextReport || nAccepted == NEvents) {
                    double percent = 100.0 * static_cast<double>(nAccepted) / static_cast<double>(NEvents);
                    double accRate = 100.0 * static_cast<double>(nAccepted) / static_cast<double>(nTried);

                    std::cout << "Progress: "
                              << std::fixed << std::setprecision(1)
                              << percent << "%  "
                              << "(" << nAccepted << "/" << NEvents << " accepted, "
                              << "tried " << nTried << ", "
                              << "acceptance " << std::setprecision(3) << accRate << "%)"
                              << std::endl;

                    nextReport += reportStep;
                }
            }
            catch (const std::exception&) {
                // skip unphysical event and continue
                continue;
            }
        }

        fout.close();

        std::cout << "\nGenerated " << nAccepted
                  << " accepted events into dvphi.lund\n";
        std::cout << "Total tried events: " << nTried << "\n";
        std::cout << "Final acceptance: "
                  << std::fixed << std::setprecision(3)
                  << 100.0 * static_cast<double>(nAccepted) / static_cast<double>(nTried)
                  << "%\n";
    }
    catch (const std::exception& ex) {
        std::cerr << "Error: " << ex.what() << "\n";
        return 1;
    }

    return 0;
}