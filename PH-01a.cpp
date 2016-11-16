#include <iostream>
#include <math.h>

const double LV = 300; // longitud del intervalo del dominio de la funcion a graficar
const int N = 1000; // número de pasos.
const double dV = LV/N;

void set_gnuplot_n();
void calculo_n();
void set_gnuplot_m();
void calculo_m();
void set_gnuplot_h();
void calculo_h();
void set_gnuplot_m_h();
void calculo_m_h();

int main()
{
  set_gnuplot_n();
  calculo_n();
  set_gnuplot_m();
  calculo_m();
  set_gnuplot_h();
  calculo_h();
  set_gnuplot_m_h();
  calculo_m_h();
  return 0;
}



void set_gnuplot_n(void)
{  
  std::cout << "set terminal svg" << std::endl; 
  std::cout << "set output 'n-punto-a.svg'" << std::endl;
  std::cout << "set title 'Steady state for function {/Helvetica-Italic n}'" << std::endl;
  std::cout << "set xlabel 'Potential (mV)'" << std::endl;
  std::cout << "set ylabel 'n'" << std::endl;
  std::cout << "plot '-' w dots" << std::endl;
}

void calculo_n(void)
{
  double n, V;
  for(int ii = 0; ii < N; ++ii){
    V = - 200 + ii*dV;
    n = ( 0.01*(V+10) ) / ( 0.01*(V+10) + 0.125 * ( exp(1+(V/10)) - 1 ) * exp(V/80) );
    std::cout << V << " " << n << std::endl;
  }
  std::cout << "e" << std::endl;
}



void set_gnuplot_m(void)
{  
  std::cout << "set terminal svg" << std::endl; 
  std::cout << "set output 'm-punto-a.svg'" << std::endl;
  std::cout << "set title 'Steady state for function {/Helvetica-Italic m}'" << std::endl;
  std::cout << "set xlabel 'Potential (mV)'" << std::endl;
  std::cout << "set ylabel 'm'" << std::endl;
  std::cout << "plot '-' w dots" << std::endl;
}

void calculo_m(void)
{
  double m, V;
  for(int ii = 0; ii < N; ++ii){
    V = - 150 + ii*2*dV/3;
    m = ( 0.01*(V+25) ) / ( 0.01*(V+25) + 4 * ( exp(2.5+(V/10)) - 1 ) * exp(V/18) );
    std::cout << V << " " << m << std::endl;
  }
  std::cout << "e" << std::endl;
}



void set_gnuplot_h(void)
{  
  std::cout << "set terminal svg" << std::endl; 
  std::cout << "set output 'h-punto-a.svg'" << std::endl;
  std::cout << "set title 'Steady state for function {/Helvetica-Italic h}'" << std::endl;
  std::cout << "set xlabel 'Potential (mV)'" << std::endl;
  std::cout << "set ylabel 'h'" << std::endl;
  std::cout << "plot '-' w dots" << std::endl;
}

void calculo_h(void)
{
  double h, V;
  for(int ii = 0; ii < N; ++ii){
    V = - 100 + ii*2*dV/3;
    h = ( 0.07*exp(V/20) ) / (  0.07*exp(V/20) + ( 1/( exp(3+(V/10)) + 1 ) ) );
    std::cout << V << " " << h << std::endl;
  }
  std::cout << "e" << std::endl;
}



void set_gnuplot_m_h(void)
{  
  std::cout << "set terminal svg" << std::endl; 
  std::cout << "set output 'm-h-punto-a.svg'" << std::endl;
  std::cout << "set title 'Steady state for function {/Helvetica-Italic m^3 * h}'" << std::endl;
  std::cout << "set xlabel 'Potential (mV)'" << std::endl;
  std::cout << "set ylabel 'm^3 * h'" << std::endl;
  std::cout << "plot '-' w dots" << std::endl;
}

void calculo_m_h(void)
{
  double h, m, coef, V;
  for(int ii = 0; ii < N; ++ii){
    V = - 150 + ii*dV;
    h = ( 0.07*exp(V/20) ) / (  0.07*exp(V/20) + ( 1/( exp(3+(V/10)) + 1 ) ) );
    m = ( 0.01*(V+25) ) / ( 0.01*(V+25) + 4 * ( exp(2.5+(V/10)) - 1 ) * exp(V/18) );
    coef = m*m*m*h;
    std::cout << V << " " << coef << std::endl;
  }
  std::cout << "e" << std::endl;
}
