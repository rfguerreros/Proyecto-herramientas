#include <iostream>
#include <cmath>
#include <cmath>
#include <fstream>

// Estructura de membrana
struct Membrane
{
  // Constantes a utilizar
  double C_m; // [uF/cm^2] Capacitancia
  double g_Na, g_K, g_L; // [(mohm*cm^2)^-1] Conductancia
  double V_Na, V_K, V_L; // [mV] Potenciales de reversa
  //Vectores con las funciones a integrar
  double f[4]; //f[0] = V, f[1] = m, f[2] = h, f[3] = n
  double ff[4]; //Vector final
  //Declaración de funciones
  void set_constants(void);
  double a_m (double V); 
  double b_m (double V);
  double a_h (double V);
  double b_h (double V);
  double a_n (double V);
  double b_n (double V);
  double m (double V);
  double h (double V);
  double n (double V);
  void set_initial_conditions(void); 
  double I_Na (double V, double m, double h); // [uA/cm^2] Corriente asociada a canales de Sodio
  double I_K (double V, double n); // [uA/cm^2] Corriente asociada a canales de Potasio
  double I_L (double V); // [uA/cm^2] Corriente asociada a filtraciones
  double I_ext (double t); // [uA/cm^2] Corriente externa proporcionada DATO DE ENTRADA
  void Df(double t, double x[], double dx[]); // Ecuaciones acopladas
  void time_advance(double ti, double tf); // Metodo de solucion numérico
};


int main(void)
{
  //Declaracion de archivos de salida
  std::ofstream sta("station.data");
  std::ofstream imp("impulso.data");
  std::ofstream can("canales.data");
  //Datos de tiempo [ms]
  double ti = 0, tf, dt = 0.001;
  double vi, vf = 100.0, dV = 0.01;
  // Declaración de membrana
  Membrane membrane;  
  membrane.set_constants();
  // Estudio de canales de activacion y su dependencia del voltaje
  for (vi = -120.0; vi <= vf;vi+=dV)
    {
      sta << vi << "\t" << membrane.m(vi) << "\t" << membrane.h(vi) << "\t" << membrane.n(vi) << std::endl;
    }
  //Condiciones iniciales
  membrane.set_initial_conditions();
  //Iteracion temporal
  for (tf = ti + dt; tf < 100.0; tf += dt)
    {
      membrane.time_advance(ti, tf);
      double V = membrane.f[0];
      double m = membrane.f[1];
      double h = membrane.f[2];
      double n = membrane.f[3];
      imp << tf << "\t" << V << "\t" << m << "\t" << h << "\t" << n << std::endl;
      can << tf << "\t" << membrane.I_Na(V,m,h) << "\t" << membrane.I_K(V,n) << "\t" << membrane.I_L(V) << "\t" << membrane.I_ext(ti) << std::endl;
      ti = tf;
    }
  
  return 0;
}

void Membrane::set_constants(void)
{
  C_m = 1.0;
  g_Na = 120.0;
  g_K = 36.0;
  g_L = 0.3;
  V_Na = 50.0;
  V_K = -77.0;
  V_L = -54.387;
}

double Membrane::a_m(double V)
{
  return 0.1 * (V + 40.0) / (1.0 - std::exp(-(V + 40.0) / 10.0));
}

double Membrane::b_m(double V)
{
  return 4.0 * std::exp(-(V + 65.0) / 18.0);
}

double Membrane::a_h(double V)
{
  return 0.07 * std::exp(-(V + 65.0)/20.0);
}

double Membrane::b_h(double V)
{
  return 1.0 / (1.0 + std::exp(-(V+35.0)/10.0));
}

double Membrane::a_n(double V)
{
  return 0.01 * (V + 55.0) / (1.0 - std::exp(-(V + 55.0)/10.0));
}

double Membrane::b_n(double V)
{
  return 0.125 * exp(-(V + 65.0) / 80.0);
}

double Membrane::I_Na (double V, double m, double h)
{
  double a = (g_Na * std::pow(m,3) * h * (V - V_Na)); 
  return a;
}

double Membrane::I_K (double V, double n)
{
  double a = (g_K * std::pow(n,4) * (V - V_K));
  return a;
}

double Membrane::I_L (double V)
{
  double a = g_L * (V - V_L);
  return a;
}

double Membrane::I_ext (double t)
{
  double i;

  if((t >= 10.0 && t <= 30.0))
    {
      i = 0.0;
    }
  else if (t >= 30 && t <= 100)
    {
      i = 0.0;
    }
  else
    {
      i = 0;
    }
  return i;
}

void Membrane::set_initial_conditions(void)
{
  f[0] = -65.0;
  f[1] = a_m(f[0])/(a_m(f[0])+b_m(f[0]));
  f[2] = a_h(f[0])/(a_h(f[0])+b_h(f[0]));
  f[3] = a_n(f[0])/(a_n(f[0])+b_n(f[0]));
}

void Membrane::Df(double t, double x[], double dx[])
{
  dx[0] = (I_ext(t) - I_Na(x[0],x[1],x[2]) - I_K(x[0],x[3]) - I_L(x[0]))/C_m;
  dx[1] = a_m(x[0]) * (1.0 - x[1]) - b_m(x[0]) * x[1];
  dx[2] = a_h(x[0]) * (1.0 - x[2]) - b_h(x[0]) * x[2];
  dx[3] = a_n(x[0]) * (1.0 - x[3]) - b_n(x[0]) * x[3];
}

void Membrane::time_advance(double ti, double tf)
{
  double x[4], dx[4], k1[4], k2[4], k3[4], k4[4];
  double h = tf - ti;
  
  Df(ti, f, dx);
  for (int ii = 0; ii < 4 ; ++ii)
    {
      k1[ii] = h * dx[ii] ;
      x[ii] = f[ii] + k1[ii]/2.0;  
    };

  Df(ti + h/2.0, x, dx);

  for (int ii = 0; ii < 4 ; ++ii)
    {
      k2[ii] = h * dx[ii] ;
      x[ii] = f[ii] + k2[ii]/2.0;  
    };

  Df(ti + h/2.0, x, dx);

  for (int ii = 0; ii < 4 ; ++ii)
    {
      k3[ii] = h * dx[ii] ;
      x[ii] = f[ii] + k3[ii];  
    };

  Df(ti + h/2.0, x, dx);

  for (int ii = 0; ii < 4 ; ++ii)
    {
      k4[ii] = h * dx[ii] ;
      ff[ii] = f[ii] + (k1[ii] + 2.0 * k2[ii] + 2.0 * k3[ii]+ k4[ii]) /6.0;  
    };
   
  for (int ii = 0; ii < 4 ; ++ii)
    {
      f[ii] = ff[ii];
    };
}

double Membrane::m (double V)
{
  return a_m(V)/(a_m(V)+b_m(V));
}

double Membrane::h (double V)
{
  return a_h(V)/(a_h(V)+b_h(V));
}

double Membrane::n (double V)
{
  return a_n(V)/(a_n(V)+b_n(V));
}
