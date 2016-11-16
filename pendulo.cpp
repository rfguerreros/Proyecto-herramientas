#include<iostream>
#include<cmath>
#include<fstream>

#define pi 3.14159265358979323846
#define g 9.80665



using namespace std;

double L = 10;
double l = 8;
double k = 1;
double m = 1;

double fx(double t,double x,double y, double vx,double vy){
  double dx;
  dx = vx;
  return dx;
}

double fy(double t,double x,double y, double vx,double vy){
  double dy;
  dy = vy;
  return dy;
}

double fvx(double t,double x,double y, double vx,double vy){
  double dvx;
  dvx = -(g/l)*x+(((k*l*l)/(L*L*m))*(y-x));
  return dvx;
}
double fvy(double t,double x,double y, double vx,double vy){
  double dvy;
  dvy = -(g/l)*y+(((k*l*l)/(L*L*m))*(x-y));
  return dvy;
}  
  
  
double rk4(double (*dx)(double,double,double,double,double),
	   double (*dy)(double,double,double,double,double),
	   double (*dvx)(double,double,double,double,double),
	   double (*dvy)(double,double,double,double,double),
	   double ti,double xi,double yi,double vxi,double vyi,double tf,double& xf,double& yf,double& vxf,double& vyf){
  double h,k1x,k2x,k3x,k4x,k1y,k2y,k3y,k4y,k1vx,k2vx,k3vx,k4vx,k1vy,k2vy,k3vy,k4vy;
  h = tf - ti;
  k1x = h * dx(ti,xi,yi,vxi,vyi);
  k1y = h * dy(ti,xi,yi,vxi,vyi);
  k1vx = h * dvx(ti,xi,yi,vxi,vyi);
  k1vy = h * dvy(ti,xi,yi,vxi,vyi);
  k2x = h * dx(ti+h/2.0,xi+k1x/2.0,yi+k1y/2.0,vxi+k1vx/2.0,vyi+k1vy/2.0);
  k2y = h * dy(ti+h/2.0,xi+k1x/2.0,yi+k1y/2.0,vxi+k1vx/2.0,vyi+k1vy/2.0);
  k2vx = h * dvx(ti+h/2.0,xi+k1x/2.0,yi+k1y/2.0,vxi+k1vx/2.0,vyi+k1vy/2.0);
  k2vy = h * dvy(ti+h/2.0,xi+k1x/2.0,yi+k1y/2.0,vxi+k1vx/2.0,vyi+k1vy/2.0);
  k3x = h * dx(ti+h/2.0,xi+k2x/2.0,yi+k2y/2.0,vxi+k2vx/2.0,vyi+k2vy/2.0);
  k3y = h * dy(ti+h/2.0,xi+k2x/2.0,yi+k2y/2.0,vxi+k2vx/2.0,vyi+k2vy/2.0);
  k3vx = h * dvx(ti+h/2.0,xi+k2x/2.0,yi+k2y/2.0,vxi+k2vx/2.0,vyi+k2vy/2.0);
  k3vy = h * dvy(ti+h/2.0,xi+k2x/2.0,yi+k2y/2.0,vxi+k2vx/2.0,vyi+k2vy/2.0);
  k4x = h * dx(ti+h,xi+k3x,yi+k3y,vxi+k3x,vyi+k3y);
  k4y = h * dy(ti+h,xi+k3x,yi+k3y,vxi+k3x,vyi+k3y);
  k4vx = h * dvx(ti+h,xi+k3x,yi+k3y,vxi+k3x,vyi+k3y);
  k4vy = h * dvy(ti+h,xi+k3x,yi+k3y,vxi+k3x,vyi+k3y);
  xf = xi + (k1x + 2.0 * (k2x+k3x) + k4x)/6.0;
  yf = yi + (k1y + 2.0 * (k2y+k3y) + k4y)/6.0;
  vxf = vxi + (k1vx + 2.0 * (k2vx+k3vx) + k4vx)/6.0;
  vyf = vyi + (k1vy + 2.0 * (k2vy+k3vy) + k4vy)/6.0;
  
}


int main (){
  
  double xi,xf,yi,yf,vxi,vxf,vyi,vyf,ti,tf,dt,tmax;
      
  ofstream rk("pendulo.data");
  
  ti = 0.0, xi = 1.0, yi = 0 , vxi = 0 , vyi = 0, dt = 0.001, tmax = 100;
    
  for(tf=ti;tf<=tmax;tf=ti+dt){
    
    rk4(fx,fy,fvx,fvy,ti,xi,yi,vxi,vyi,tf,xf,yf,vxf,vyf);
         
    rk << tf <<"  " << xf <<"  " << yf  <<"  " << xf+yf <<"  " << xf-yf<< endl;

    
    xi = xf;
    yi = yf;
    vxi = vxf;
    vyi = vyf;
    ti = tf;
    
    
  }
  
}  
  
  
