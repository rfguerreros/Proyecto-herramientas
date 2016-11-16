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


double fV(double V, double n, double m, double h, double t);
double fn(double n, double V, double t);
double fm(double n, double V, double t);
double fh(double n, double V, double t);


int main(void)
{
  return 0;
}
