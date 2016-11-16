#include <iostream>
#include <cstdlib>
#include <fstream>
#include <cmath>

const double C = 1.0; //micro F/cm*cm 
const double gK = 36.0; // mmho/cm*cm
const double gNa = 120.0; // mmho/cm*cm
const double gl = 0.3; // mmho/cm*cm
const double VK = 12.0; // mV
const double VNa = -115.0; // mV
const double Vl = 10.6; // mV

double fV(double V, double n, double m, double h, double t);
double fn(double n, double V, double t);
double fm(double m, double V, double t);
double fh(double h, double V, double t);
double rk4(double (*fV)(double ,double ,double, double, double),
	   double (*fn)(double ,double ,double),
	   double (*fm)(double ,double ,double),
	   double (*fh)(double ,double ,double),
	   double ti, double vi, double ni, double mi, 
	   double hi, double tf, double & vf, double & nf,
	   double & mf, double & hf);

int main(void)
{
  double vi = 0, ni = 0, mi = 0, hi = 0, 
    vf = 0, nf = 0, mf = 0, hf = 0, 
    ti = 0, tf = 0, tmax = 50;
  
  std::ofstream imp("impulso.data");

  double dt = 0.001;
  
  for (tf = ti + dt; tf < tmax; tf = ti + dt){
    rk4(fV, fn, fm, fh, ti, vi, ni, mi, hi, tf, vf, nf, mf, hf);
    
    imp << tf << "\t" << vf << std::endl;
    
    vi = vf;
    ni = nf;
    mi = mf;
    hi = hf;
    ti = tf;
  }
  
  return 0;
}


double fV(double V, double n, double m, double h, double t)
{
  double dV;
  double i = 0;
  dV = (-gK * std::pow(n,4)*(V - VK)- gNa * std::pow(m,3)* h *(V-VNa) - gl * (V-Vl) + i) / C;  
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

double rk4(double (*fV)(double ,double ,double, double, double),
	   double (*fn)(double ,double ,double),
	   double (*fm)(double ,double ,double),
	   double (*fh)(double ,double ,double),
	   double ti, double vi, double ni, double mi, 
	   double hi, double tf, double & vf, double & nf,
	   double & mf, double & hf)
{
  double delta = tf - ti;
  double k1v, k2v, k3v, k4v;
  double k1n, k2n, k3n, k4n;
  double k1m, k2m, k3m, k4m;
  double k1h, k2h, k3h, k4h;
  
  k1v = delta * fV(vi, ni, mi, hi, ti);
  k1n = delta * fn(ni, vi, ti);
  k1m = delta * fm(mi, vi, ti);
  k1h = delta * fh(hi, vi, ti);
  
  k2v = delta * fV(vi + k1v/2.0, ni + k1n/2.0,
		   mi + k1m/2.0, hi + k1h/2.0, ti + delta/2.0);
  k2n = delta * fn(ni + k1n/2.0, vi + k1v/2.0, ti + delta/2.0);
  k2m = delta * fm(mi + k1m/2.0, vi + k1v/2.0, ti + delta/2.0);
  k2h = delta * fh(hi + k1h/2.0, vi + k1v/2.0, ti + delta/2.0);

  k3v = delta * fV(vi + k2v/2.0, ni + k2n/2.0,
		   mi + k2m/2.0, hi + k2h/2.0, ti + delta/2.0);
  k3n = delta * fn(ni + k2n/2.0, vi + k2v/2.0, ti + delta/2.0);
  k3m = delta * fm(mi + k2m/2.0, vi + k2v/2.0, ti + delta/2.0);
  k3h = delta * fh(hi + k2h/2.0, vi + k2v/2.0, ti + delta/2.0);
  
  k4v = delta * fV(vi + k3v, ni + k3n, mi + k3m, hi + k3h, ti + delta);
  k4n = delta * fn(ni + k3n, vi + k3v, ti + delta);
  k4m = delta * fm(mi + k3m, vi + k3v, ti + delta);
  k4h = delta * fh(hi + k3h, vi + k3v, ti + delta);

  vf = vi + (k1v + 2.0 * (k2v + k3v) + k4v)/6.0;
  nf = ni + (k1n + 2.0 * (k2n + k3n) + k4n)/6.0;
  mf = mi + (k1m + 2.0 * (k2m + k3m) + k4m)/6.0;
  hf = hi + (k1m + 2.0 * (k2h + k3h) + k4h)/6.0;
  
}
