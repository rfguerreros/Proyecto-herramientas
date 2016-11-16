#include <iostream>
#include <cstdlib>
#include <fstream>
#include <cmath>

const double C = 1.0; //micro F/cm*cm 
const double gK = 36; // mmho/cm*cm
const double gNa = 120; // mmho/cm*cm
const double gl = 0.3; // mmho/cm*cm
const double VK = 12; // mV
const double VNa = -115; // mV
const double Vl = 10.6; // mV

double I(double t);
double fV(double V, double n, double m, double h, double t, double (*I) (double));
double fn(double n, double V, double t);
double fm(double m, double V, double t);
double fh(double h, double V, double t);


int main(void)
{
  return 0;
}

double I(double t)
{
  double i;
  i = 1;
  return i;
}

double fV(double V, double n, double m, double h, double t, double (*I) (double))
{
  double dV;
  dV = -gK * std::pow(n,4)*(V-VK)- gNa * std::pow(m,3)*(V-VNa) -gK * (V-Vl) + I(t);  
  return dV;
}

double fn(double n, double V, double t)
{
  double dn;
  double an, bn;
  an = 0.01 * (V + 10.0) / (exp(1.0 + V / 10.0) - 1.0);
  bn = 0.125 * exp(V / 80.0);
  dn = an * (1.0 - n) - bn * n;
  return dn;
}

double fm(double m, double V, double t)
{
  double dm;
  double am, bm;
  am = 0.01 * (V + 25.0) / (exp(2.5 + V / 10.0) - 1.0);
  bm = 4.0 * exp(V / 18.0);
  dm = am * (1.0 - m) - bm * m;
  return dm;
}

double fh(double h, double V, double t)
{
  double dh;
  double ah, bh;
  ah = 0.07 * exp(V / 20.0);
  bh = 1.0 / (exp(3.0 + V / 10.0) + 1.0);
  dh = ah * (1.0 - h) - bh * h;
  return dh;
}
