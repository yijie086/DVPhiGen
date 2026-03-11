#ifndef RANDOM_H
#define RANDOM_H

class RandomGenerator {
public:
    RandomGenerator(unsigned int seed = 12345);

    // basic random
    double Uniform(double a, double b);

    // set generation ranges
    void SetRanges(double Q2_min_in, double Q2_max_in,
                   double W_min_in,  double W_max_in, double EB_in);

    // proposal distributions
    double GenerateQ2();
    double GenerateW();
    double GenerateT(double Q2, double W);

    // angular variables
    double GeneratePhiE();
    double GeneratePhi();
    double GeneratePhiH();

    
    double GenerateCosThetaH(double EB, double Q2, double W, double t);
    double WeightCosThetaH(double EB, double Q2, double W, double t, double cosThetaH) const;

    // target weight in (Q2, W, t)
    double WeightQ2Wt(double EB, double Q2, double W, double t) const;

    // generate one accepted kinematic point
    void GenerateKinematics(double EB, double& Q2, double& W, double& t);

private:
    unsigned int state;

    double Q2_min;
    double Q2_max;
    double W_min;
    double W_max;
    double wmax;
    double EB;

    unsigned int NextUInt();

    // helpers
    static double Lambda(double x, double y, double z);

    void ValidateRanges() const;
    void UpdateWmax(double EB);

    void GetTRange(double Q2, double W, double& t_low, double& t_high) const;
};

#endif