//
//  ray_tracing.cpp
//  
//
//  Created by Andrew Mummery on 18/01/2019.
//

#include <iostream>
#include <cmath>
#include <vector>
#include <fstream>
#include <ctime>
#include <algorithm>
#include <iomanip>
#include <stdlib.h>
#include <sstream>
#include <stdexcept>

using namespace std;

class functions{
    
public:
    
    // These functions are short hand definitions in the metric .
    
    double Delta(double r, double a){
        double delta = r*r - 2*r + a*a;
        return delta;
    }
    
    double Sigma(double r, double a, double theta){
        double sigma = r*r + a*a*cos(theta)*cos(theta);
        return sigma;
    }
    
    double A(double r, double a, double theta){
        double ay = (r*r + a*a)*(r*r + a*a) - a*a*Delta(r,a)*sin(theta)*sin(theta);
        return ay;
    }
    
    // These functions are the Christoffel coefficients for the
    // Kerr metric in Boyer-Lindquist co-ordinates
    
    // We begin with the coefficients important for the evolution of the radial co-ordinate
    
    // These are the quadratic terms
    
    double crtt(double r, double a, double theta){
        double CRTT = Delta(r,a)*(r*r - a*a*cos(theta)*cos(theta))/(Sigma(r,a,theta)*Sigma(r,a,theta)*Sigma(r,a,theta));
        return CRTT;
    }
    
    double crrr(double r, double a, double theta){
        double CRRR = (r*a*a*sin(theta)*sin(theta)-(r*r - a*a*cos(theta)*cos(theta)))/(Sigma(r,a,theta)*Delta(r,a));
        return CRRR;
    }
    
    double croo(double r, double a, double theta){
        double CROO = - r*Delta(r,a)/Sigma(r,a,theta);
        return CROO;
    }
    
    double crpp(double r, double a, double theta){
        double CRPP =  Delta(r,a)*sin(theta)*sin(theta)*(-r*Sigma(r,a,theta)*Sigma(r,a,theta) + a*a*sin(theta)*sin(theta)*(r*r - a*a*cos(theta)*cos(theta)))/(Sigma(r,a,theta)*Sigma(r,a,theta)*Sigma(r,a,theta));
        return CRPP;
    }
    
    // and these are the cross terms
    
    double crtp(double r, double a, double theta){
        double CRTP = - Delta(r,a)*a*sin(theta)*sin(theta)*(r*r - a*a*cos(theta)*cos(theta))/(Sigma(r,a,theta)*Sigma(r,a,theta)*Sigma(r,a,theta));
        return CRTP;
    }
    
    double  crro(double r, double a, double theta){
        double CRRO = - a*a*sin(theta)*cos(theta)/(Sigma(r,a,theta));
        return CRRO;
    }
    
    // These are the coefficients important for the evolution of the theta co-ordinate
    
    // These are the quadratic terms
    
    double  cott(double r, double a, double theta){
        double COTT =  - 2*a*a*r*sin(theta)*cos(theta)/(Sigma(r,a,theta)*Sigma(r,a,theta)*Sigma(r,a,theta));
        return COTT;
    }
    
    double corr(double r, double a, double theta){
        double CORR = a*a*sin(theta)*cos(theta)/(Sigma(r,a,theta)*Delta(r,a));
        return CORR;
    }
    
    double cooo(double r, double a, double theta){
        double COOO = - a*a*sin(theta)*cos(theta)/Sigma(r,a,theta);
        return COOO;
    }
    
    double copp(double r, double a, double theta){
        double COPP = - sin(theta)*cos(theta)*(A(r,a,theta)*Sigma(r,a,theta) + 2*(r*r+a*a)*a*a*r*sin(theta)*sin(theta))/(Sigma(r,a,theta)*Sigma(r,a,theta)*Sigma(r,a,theta));
        return COPP;
    }
    
    // and these are the cross terms
    
    double cotp(double r, double a, double theta){
        double COTP = 2*a*r*(r*r + a*a)*sin(theta)*cos(theta)/(Sigma(r,a,theta)*Sigma(r,a,theta)*Sigma(r,a,theta));
        return COTP;
    }
    
    double coro(double r, double a, double theta){
        double CORO = r/Sigma(r,a,theta);
        return CORO;
    }
    
    // These are the metric coefficients for the Kerr metric
    // in Boyer-Lindquist co-ordinates
    
    double gtt(double r, double a, double theta){
        double GTT = -(1-2*r/Sigma(r,a,theta));
        return GTT;
    }
    
    double gtp(double r, double a, double theta){
        double GTP = - 2*a*r*sin(theta)*sin(theta)/Sigma(r,a,theta);
        return GTP;
    }
    
    double grr(double r, double a, double theta){
        double GRR = Sigma(r,a,theta)/Delta(r,a);
        return GRR;
    }
    
    double goo(double r, double a, double theta){
        double GOO = Sigma(r,a,theta);
        return GOO;
    }
    
    double gpp(double r, double a, double theta){
        double GPP = (r*r + a*a + 2*a*a*r*sin(theta)*sin(theta)/Sigma(r,a,theta))*sin(theta)*sin(theta);
        return GPP;
    }
    
    // This function tests for error propogating throughout the
    // integration, for a photon this should equal -1 for all times
    
    double squig(double r, double a , double theta, double td, double rd, double od, double pd){
        double SQUIG = (grr(r,a,theta)*rd*rd + gpp(r,a,theta)*pd*pd + goo(r,a,theta)*od*od + 2*gtp(r,a,theta)*td*pd)/(gtt(r,a,theta)*td*td);
        return SQUIG;
    }
    
    // these functions calculate the angular momentum and
    // energy of a given photon, these two quantities are
    // conserved along geodesics.
    
    double E(double r, double a, double theta, double td, double pd){
        double ee = - gtt(r,a,theta) * td - gtp(r,a,theta) * pd;
        return ee;
    }
    
    double L(double r, double a, double theta, double td, double pd){
        double ll = gpp(r,a,theta) * pd + gtp(r,a,theta) * td ;
        return ll;
    }
    
    double j(double r, double a, double theta, double td, double pd){
        double jj = L(r,a,theta,td,pd)/E(r,a,theta,td,pd);
        return jj;
    }
};



class Initial_conditions{
    
public:
    
    // These functions take a position in the image plane (alpha0, beta0)
    // at an angle theta_i from the black hole center,
    // and turn that into initial conditions on the Boyer-Lindquist
    // co-ordinates (r, theta, phi), and also generate the initial
    // photon 4-velocity (u^t, u^r, u^theta, u^phi)
    
    // Position initialisation:
    
    functions f;
    
    double r_init(double d, double a0, double b0){
        double R_INIT = sqrt(d*d + a0*a0 + b0*b0);
        return R_INIT;
    }
    
    double theta_init(double d, double a0, double b0, double theta_i){
        double THETA_INIT = acos( (d*cos(theta_i)+ b0*sin(theta_i))/r_init(d,a0,b0) );
        return THETA_INIT;
    }
    
    double phi_init(double d, double a0, double b0, double theta_i){
        double PHI_INIT = atan(a0/(d*sin(theta_i)-b0*cos(theta_i)));
        return PHI_INIT;
    }
    
    // 4-velcoity initialisation:
    
    double vr_init(double d, double a0, double b0){
        double VR_INIT =  - d/r_init(d,a0,b0);
        return VR_INIT;
    }
    
    double vo_init(double d, double a0, double b0, double theta_i){
        double VO_INIT =  - (-cos(theta_i) + d/(pow(r_init(d,a0,b0),2)) * (d*cos(theta_i)+b0*sin(theta_i)))/sqrt(pow(r_init(d,a0,b0),2) - pow((d*cos(theta_i)+b0*sin(theta_i)),2) );
        return VO_INIT;
    }
    
    double vp_init(double d, double a0, double b0, double theta_i){
        double VP_INIT = - (-a0*sin(theta_i)/(pow(d*sin(theta_i)-b0*cos(theta_i),2) + a0*a0));
        return VP_INIT;
    }
    
    double vt_init(double a, double d, double a0, double b0, double theta_i){
        double aaa = f.gtt(r_init(d,a0,b0),a,theta_i);
        double bbb = 2*f.gtp(r_init(d,a0,b0),a,theta_i)*vp_init(d,a0,b0,theta_i);
        double ccc  = f.grr(r_init(d,a0,b0),a,theta_i)*vr_init(d,a0,b0)*vr_init(d,a0,b0) + f.goo(r_init(d,a0,b0),a,theta_i)*vo_init(d,a0,b0,theta_i)*vo_init(d,a0,b0,theta_i) + f.gpp(r_init(d,a0,b0),a,theta_i)*vp_init(d,a0,b0,theta_i)*vp_init(d,a0,b0,theta_i);
        
        double VT_INIT = (- bbb - sqrt(bbb*bbb - 4*aaa*ccc))/(2*aaa); // this is the physical version as it is > 0
        return VT_INIT;
    }
};

class Equations_of_motion{
    
    
public:
    
    // These functions give the 4 equations of motion
    // for the photon on its path. There are two first
    // order equations, the ones which govern the evolution
    // of phi and t. r and theta evolve according to
    // the second-order geodesic equations.
    
    // For the second order equations we solve it as
    // a first order equation for the velocities,
    // and then use these velocities to calculate the
    // change of the co-ordinate positions.
    
    functions f;
    
    double tdot(double r, double a, double theta,double j){
        double TDOT = - (f.gpp(r,a,theta) + j*f.gtp(r,a,theta))/(f.gpp(r,a,theta)*f.gtt(r,a,theta) - f.gtp(r,a,theta)*f.gtp(r,a,theta));
        return TDOT;
    }
    
    double pdot(double r, double a, double theta, double j){
        double PDOT = (j*f.gtt(r,a,theta) + f.gtp(r,a,theta))/(f.gpp(r,a,theta)*f.gtt(r,a,theta) - f.gtp(r,a,theta)*f.gtp(r,a,theta));
        return PDOT;
    }
    
    double vrdot(double r, double a, double theta, double td, double rd, double od, double pd){
        double VRDOT = (- f.crtt(r,a,theta)*td*td - f.crrr(r,a,theta)*rd*rd - f.croo(r,a,theta)*od*od - f.crpp(r,a,theta)*pd*pd - 2*f.crtp(r,a,theta)*td*pd - 2*f.crro(r,a,theta)*rd*od);
        return VRDOT;
    }
    
    double vodot(double r, double a, double theta, double td, double rd, double od, double pd){
        double VODOT = - f.cott(r,a,theta)*td*td - f.corr(r,a,theta)*rd*rd - f.cooo(r,a,theta)*od*od - f.copp(r,a,theta)*pd*pd - 2*f.cotp(r,a,theta)*td*pd - 2*f.coro(r,a,theta)*rd*od;
        return VODOT;
    }
    
};


class Integrators{
    
    // This set of functions gives the implementation
    // of the Runge-Kutta 4 algorithm.
    
public:
    Equations_of_motion eom;
    functions f;
    
    void k1(double r, double a, double theta, double td, double rd, double od, double pd, double j, double h, double* kk1){
        kk1[0] = eom.vrdot(r,a,theta,td,rd,od,pd);
        kk1[1] = eom.vodot(r,a,theta,td,rd,od,pd);
        kk1[2] = eom.tdot(r,a,theta,j);
        kk1[3] = rd;
        kk1[4] = od;
        kk1[5] = eom.pdot(r,a,theta,j);
    }
    
    void k2(double r, double a, double theta, double td, double rd, double od, double pd, double j, double h, double* kk1, double* kk2){
        kk2[0] = eom.vrdot(r + h*kk1[3]/2 ,a ,theta + h*kk1[4]/2 , eom.tdot(r+ h*kk1[3]/2, a, theta + h*kk1[4]/2,j ),rd + h*kk1[0]/2,od + h*kk1[1]/2, eom.pdot(r + h*kk1[3]/2, a, theta + h*kk1[4]/2,j ));
        kk2[1] = eom.vodot(r + h*kk1[3]/2 ,a ,theta + h*kk1[4]/2 , eom.tdot(r+ h*kk1[3]/2, a, theta + h*kk1[4]/2,j ),rd + h*kk1[0]/2,od + h*kk1[1]/2, eom.pdot(r + h*kk1[3]/2, a, theta + h*kk1[4]/2,j ));
        kk2[2] = eom.tdot(r + h*kk1[3]/2, a, theta + h*kk1[4]/2,j );
        kk2[3] = rd + h*kk1[0]/2;
        kk2[4] = od + h*kk1[1]/2;
        kk2[5] = eom.pdot(r + h*kk1[3]/2, a, theta + h*kk1[4]/2,j );
    }
    
    void k3(double r, double a, double theta, double td, double rd, double od, double pd, double j, double h, double* kk2, double* kk3){
        kk3[0] = eom.vrdot(r + h*kk2[3]/2 ,a ,theta + h*kk2[4]/2 , eom.tdot(r+ h*kk2[3]/2, a, theta + h*kk2[4]/2,j ),rd + h*kk2[0]/2,od + h*kk2[1]/2, eom.pdot(r + h*kk2[3]/2, a, theta + h*kk2[4]/2,j ));
        kk3[1] = eom.vodot(r + h*kk2[3]/2 ,a ,theta + h*kk2[4]/2 , eom.tdot(r+ h*kk2[3]/2, a, theta + h*kk2[4]/2,j ),rd + h*kk2[0]/2,od + h*kk2[1]/2, eom.pdot(r + h*kk2[3]/2, a, theta + h*kk2[4]/2,j ));
        kk3[2] = eom.tdot(r + h*kk2[3]/2, a, theta + h*kk2[4]/2,j );
        kk3[3] = rd + h*kk2[0]/2;
        kk3[4] = od + h*kk2[1]/2;
        kk3[5] = eom.pdot(r + h*kk2[3]/2, a, theta + h*kk2[4]/2,j );
    }
    
    void k4(double r, double a, double theta, double td, double rd, double od, double pd, double j, double h, double* kk3, double* kk4){
        kk4[0] = eom.vrdot(r + h*kk3[3] ,a ,theta + h*kk3[4] , eom.tdot(r+ h*kk3[3], a, theta + h*kk3[4],j ),rd + h*kk3[0],od + h*kk3[1], eom.pdot(r + h*kk3[3], a, theta + h*kk3[4],j ));
        kk4[1] = eom.vodot(r + h*kk3[3] ,a ,theta + h*kk3[4] , eom.tdot(r+ h*kk3[3], a, theta + h*kk3[4],j ),rd + h*kk3[0],od + h*kk3[1], eom.pdot(r + h*kk3[3], a, theta + h*kk3[4],j ));
        kk4[2] = eom.tdot(r + h*kk3[3], a, theta + h*kk3[4],j );
        kk4[3] = rd + h*kk3[0];
        kk4[4] = od + h*kk3[1];
        kk4[5] = eom.pdot(r + h*kk3[3], a, theta + h*kk3[4],j );
    }
    
    void evolve_one_time_step(double* Y, double a, double j, double h, double* kk1, double* kk2, double* kk3, double* kk4){
        double r = Y[3];
        double theta = Y[4];
        double rd = Y[0];
        double od = Y[1];
        
        double td = eom.tdot(r, a, theta, j);
        double pd = eom.pdot(r, a, theta, j);
        
        k1(r, a, theta, td, rd, od, pd, j, h, kk1);
        k2(r, a, theta, td, rd, od, pd, j, h, kk1, kk2);
        k3(r, a, theta, td, rd, od, pd, j, h, kk2, kk3);
        k4(r, a, theta, td, rd, od, pd, j, h, kk3, kk4);
        
        Y[0] = Y[0] + 1.0/6.0*h*(kk1[0] + 2*kk2[0] + 2*kk3[0] + kk4[0]);
        Y[1] = Y[1] + 1.0/6.0*h*(kk1[1] + 2*kk2[1] + 2*kk3[1] + kk4[1]);
        Y[2] = Y[2] + 1.0/6.0*h*(kk1[2] + 2*kk2[2] + 2*kk3[2] + kk4[2]);
        Y[3] = Y[3] + 1.0/6.0*h*(kk1[3] + 2*kk2[3] + 2*kk3[3] + kk4[3]);
        Y[4] = Y[4] + 1.0/6.0*h*(kk1[4] + 2*kk2[4] + 2*kk3[4] + kk4[4]);
        Y[5] = Y[5] + 1.0/6.0*h*(kk1[5] + 2*kk2[5] + 2*kk3[5] + kk4[5]);
    }
    
    double get_next_time_step(double a0, double rad, double vrad, double thet, double vthet, double ph, double vph, double F){
        // gets the next time-step in the algorithm.
        double elm[3];
        double min;
        double dT;
        if (a0 != 0){// dT is set by F multiplied by the shortest variability time-scale in the problem.
            elm[0] = abs(rad/vrad);
            elm[1] = abs(thet/vthet);
            elm[2] = abs(ph/vph);
            min = elm[0];
            for (int m = 0; m < 2; m ++){
                if (elm[m] < min){
                    min = elm[m];
                }
            }
            dT = F*min;
        }
        else{//If a0 = 0, then \phi = 0, and we will pass over the pole \theta = 0, at some point in the trajectory.
            // We must only consider the radial time-scale in that case.
            dT = F*abs(rad/vrad);
        }
    
        return dT;
    }
};

class One_Photon{
    
  // This will solve the geodesic equations for a single photon
    
public:
    functions f;
    Equations_of_motion eom;
    Integrators in;
    Initial_conditions start;
    
    void ray_trace_output(vector<double> arr, int len_arr){
        // Takes a vector input arr, with length len_arr
        // and outputs it as one big string.
        //
        // This string can be converted into a numpy array
        // quite easily within python.
        std::stringstream ss;
        for(int i = 0; i < len_arr; ++i)
        {
          if(i != 0)
            ss << ",";
          ss << arr[i];
        }
        std::string s = ss.str();
        std::cout << s << std::endl;
    }
    
    double isco_a_pos(double a){
        // Returns the ISCO radius for prograde spins, a > 0.
        double Z_1 = 1 + pow(1-a*a,1.0/3.0) * (pow(1+a,1.0/3.0) + pow(1-a,1.0/3.0));
        double Z_2 = sqrt(3*a*a + Z_1*Z_1);
        double rISCO =  (3 + Z_2 -  sqrt((3-Z_1)*(3 + Z_1 + 2 * Z_2)));
        return rISCO;
    }

    double isco_a_neg(double a){
        // Returns the ISCO radius for retrograde spins, a < 0.
        double Z_1 = 1 + pow(1-a*a,1.0/3.0) * (pow(1+a,1.0/3.0) + pow(1-a,1.0/3.0));
        double Z_2 = sqrt(3*a*a + Z_1*Z_1);
        double rISCO =  (3 + Z_2 + sqrt((3-Z_1)*(3 + Z_1 + 2 * Z_2)));
        return rISCO;
    }

    double get_isco(double a){
    // Gets the ISCO radius for a given BH spin parameter
        double ri = 1.0;
        if (a == 0){
            ri = 6;// Schwarzschild result.
        }
        else if (a > 0){
            ri = isco_a_pos(a);// Analytic formula for the ISCO depends on the sign of the spin parameter.
        }
        else {
            ri = isco_a_neg(a);// Analytic formula for the ISCO depends on the sign of the spin parameter.
        }
        return ri;
    }
    
    double get_event_horizon(double a){
        // Gets the event horizon for a given BH spin parameter.
        double rH =  1 + sqrt(1-a*a);
        return rH;
    }

    void one_orbit(double a, double a0, double b0, double theta_i, double d, double F, const int n, double& rf, double& phif, double& thetaf, std::string technique){
        // First we initialise the photons 4-velocities.
        
        double vt0 = start.vt_init(a,d,a0,b0,theta_i);
        double vr0 = start.vr_init(d,a0,b0);
        double vo0 = start.vo_init(d,a0,b0,theta_i);
        double vp0 = start.vp_init(d,a0,b0,theta_i);
        
        // then we initialise the photons position
        
        double t0 = 1;
        double r0 = start.r_init(d,a0,b0);
        double theta0 = start.theta_init(d,a0,b0,theta_i);
        double phi0 = start.phi_init(d,a0,b0,theta_i);
        
        // calculate the photons normalised angular momentum
        // the important quantity for red-shift calculations
        
        double j = f.j(r0,a,theta0,vt0,vp0);
        
        // calculate initial r and theta acceleration
        
        double vrdot0 = eom.vrdot(r0,a,theta0,vt0,vr0,vo0,vp0);
        double vodot0 = eom.vodot(r0,a,theta0,vt0,vr0,vo0,vp0);
        
        // important vectors for algorithm implimentation
        
        vector<double> t(1);
        double dT = 0;
        
        // Varying time-stepping, time-step in integration is
        // a fixed fraction of the time-scale of the
        // rapidest changing variable
        
        // d tau == F * min( r/r_dot, theta/theta_dot, phi/phi_dot )
        
        dT = in.get_next_time_step(a0, r0, vr0, theta0, vo0, phi0, vp0, F);
        t[0] = dT;
        
        // set-up for the algorithm implimentation
        
        double rk4y[6];
        double kk1[6];
        double kk2[6];
        double kk3[6];
        double kk4[6];
        
        
        rk4y[0] = vr0;
        rk4y[1] = vo0;
        rk4y[2] = t0;
        rk4y[3] = r0;
        rk4y[4] = theta0;
        rk4y[5] = phi0;
        
        
        kk1[0] = kk1[1] = kk1[2] = kk1[3] = kk1[4] = kk1[5] = 0;
        kk2[0] = kk2[1] = kk2[2] = kk2[3] = kk2[4] = kk2[5] = 0;
        kk3[0] = kk3[1] = kk3[2] = kk3[3] = kk3[4] = kk3[5] = 0;
        kk4[0] = kk4[1] = kk4[2] = kk4[3] = kk4[4] = kk4[5] = 0;
        
        // Quantities calculated at every time-step
        
        vector<double> r(1);
        vector<double> theta(1);
        vector<double> phi(1);


        r[0] = r0;
        theta[0] = theta0;
        phi[0] = phi0;

        int m = 0;// For watching the time-steps.
        
        // First technique, all photon paths terminated when they move through the disc plane (z = 0).
        // Even if they pass between the ISCO and event horizon.
        if (technique == "Simple"){
            double rH = get_event_horizon(a);
            while (rk4y[3]*cos(rk4y[4]) > 0){ // stops when z < 0
                if (m < n){ // unless algorithm has got stuck (generally happens when photon hits black-hole)
                    if (rk4y[3] > rH){ // if r < rH then  in the black-hole
                        in.evolve_one_time_step(rk4y,a,j,dT,kk1,kk2,kk3,kk4); // implements the RK4 algorithm
                        r.push_back(rk4y[3]);
                        theta.push_back(rk4y[4]); // updates the variables
                        phi.push_back(rk4y[5]);
                        dT = in.get_next_time_step(a0, rk4y[3], rk4y[0], rk4y[4], rk4y[1], rk4y[5], eom.pdot(rk4y[3],a,rk4y[4],j), F);
                        t.push_back(t[m]+dT);
                    }
                    else{//I am in the black hole.
                        break;
                    }
                    m = m + 1; // in case algorithm gets stuck
                }
                else{// Too many steps.
                    break;
                }
            }
        }
        
        // Second technique, allows light rays to pass between the disc and event horizon: rH < r < rI.
        if (technique == "Disc"){
            double rI = get_isco(a);
            double rH = get_event_horizon(a);
            while (rk4y[3]*cos(rk4y[4]) > -20){
                if (m < n){
                    if  (rk4y[0]*cos(rk4y[4]) - rk4y[3]*sin(rk4y[4])*rk4y[1] < 0){
                        if (rk4y[3] > rH){
                            if (abs(rk4y[3]*cos(rk4y[4])) > 0.05){
                                    in.evolve_one_time_step(rk4y,a,j,dT,kk1,kk2,kk3,kk4); // implements the RK4 algorithm
                                    r.push_back(rk4y[3]);
                                    theta.push_back(rk4y[4]); // updates the variables
                                    phi.push_back(rk4y[5]);
                                    dT = in.get_next_time_step(a0, rk4y[3], rk4y[0], rk4y[4], rk4y[1], rk4y[5], eom.pdot(rk4y[3],a,rk4y[4],j), F);
                                    t.push_back(t[m]+dT);
                                    m = m + 1;
                            }
                            else if (rk4y[3] < rI){
                                    in.evolve_one_time_step(rk4y,a,j,dT,kk1,kk2,kk3,kk4); // implements the RK4 algorithm
                                    r.push_back(rk4y[3]);
                                    theta.push_back(rk4y[4]); // updates the variables
                                    phi.push_back(rk4y[5]);
                                    dT = in.get_next_time_step(a0, rk4y[3], rk4y[0], rk4y[4], rk4y[1], rk4y[5], eom.pdot(rk4y[3],a,rk4y[4],j), F);
                                    t.push_back(t[m]+dT);
                                    m = m + 1;
                                    }
                            else{
                                break;
                            }
                        }
                        else{
                            break;
                        }
                    }
                    else if (rk4y[3]*cos(rk4y[4]) < 20){
                        if (rk4y[3] > rH){
                            if (abs(rk4y[3]*cos(rk4y[4])) > 0.1){
                                    in.evolve_one_time_step(rk4y,a,j,dT,kk1,kk2,kk3,kk4); // implements the RK4 algorithm
                                    r.push_back(rk4y[3]);
                                    theta.push_back(rk4y[4]); // updates the variables
                                    phi.push_back(rk4y[5]);
                                    dT = in.get_next_time_step(a0, rk4y[3], rk4y[0], rk4y[4], rk4y[1], rk4y[5], eom.pdot(rk4y[3],a,rk4y[4],j), F);
                                    t.push_back(t[m]+dT);
                                    m = m + 1;
                            }
                                else if (rk4y[3] < rI){
                                    in.evolve_one_time_step(rk4y,a,j,dT,kk1,kk2,kk3,kk4); // implements the RK4 algorithm
                                    r.push_back(rk4y[3]);
                                    theta.push_back(rk4y[4]); // updates the variables
                                    phi.push_back(rk4y[5]);
                                    dT = in.get_next_time_step(a0, rk4y[3], rk4y[0], rk4y[4], rk4y[1], rk4y[5], eom.pdot(rk4y[3],a,rk4y[4],j), F);
                                    t.push_back(t[m]+dT);
                                    m = m + 1;
                                }
                                else{
                                    break;
                                }
                        }
                    else{
                        break;
                    }
                    }
                else{
                    break;
                }
            

            }
            }
        }
        
        // Third technique, assumes no disc in the z = 0 plane.
        if (technique == "NoDisc"){
            double rH = get_event_horizon(a);
            while ( (rk4y[3]*cos(rk4y[4]) > -20) && (rk4y[3]*sin(rk4y[4])*cos(rk4y[5]) > -30)){// stops when photon reaches z = -20 or x = -30, or x > 30 and v_x > 0.
                if ( ((rk4y[0]*sin(rk4y[4])*cos(rk4y[5]) + rk4y[3] * rk4y[1] * cos(rk4y[4])*cos(rk4y[5]) - rk4y[3] * eom.pdot(rk4y[3],a,rk4y[4],j) * sin(rk4y[4])*sin(rk4y[5]) > 0)  && rk4y[3]*sin(rk4y[4])*cos(rk4y[5]) > 30) ){
                    break;// Breaks if x > 30 and v_x > 0.
                }
                if (m < n){
                    if  (rk4y[0]*cos(rk4y[4]) - rk4y[3]*sin(rk4y[4])*rk4y[1] < 0){// if z-velocity is negative.
                        if (rk4y[3] > rH){
                            in.evolve_one_time_step(rk4y,a,j,dT,kk1,kk2,kk3,kk4); // implements the RK4 algorithm
                            r.push_back(rk4y[3]);
                            theta.push_back(rk4y[4]); // updates the variables
                            phi.push_back(rk4y[5]);
                            dT = in.get_next_time_step(a0, rk4y[3], rk4y[0], rk4y[4], rk4y[1], rk4y[5], eom.pdot(rk4y[3],a,rk4y[4],j), F);
                            t.push_back(t[m]+dT);
                            m = m + 1;
                        }
                        else{
                            break;
                        }
                    }
                    else if (rk4y[3]*cos(rk4y[4]) < 20){// else I am going up, and z < + 20.
                        if (rk4y[3] > rH){
                            in.evolve_one_time_step(rk4y,a,j,dT,kk1,kk2,kk3,kk4); // implements the RK4 algorithm
                            r.push_back(rk4y[3]);
                            theta.push_back(rk4y[4]); // updates the variables
                            phi.push_back(rk4y[5]);
                            dT = in.get_next_time_step(a0, rk4y[3], rk4y[0], rk4y[4], rk4y[1], rk4y[5], eom.pdot(rk4y[3],a,rk4y[4],j), F);
                            t.push_back(t[m]+dT);
                            m = m + 1;
                        }
                    else{
                        break;
                    }
                    }
                else{
                    break;
                }
            

            }
            }
        }

        
        ray_trace_output(r, m);
        ray_trace_output(theta, m);
        ray_trace_output(phi, m);
        ray_trace_output(t, m);
        
        std::cout << j << std::endl;
        
        rf = r[m];
        thetaf = theta[m]; // Edits these variables in the main script
        phif = phi[m];
        
        
    }
};

int main(int argc, const char * argv[]) {
    // We import the class of one photon orbit integrators
    
    One_Photon onp;
    
    if (argc != 6){// Should get 1: C++ file name (always happens), 2-5: Parameters alpha0, beta0, theta0, a & 6: technique to be used.
        throw runtime_error("Not enough input parameters");
        return -1;
    }
    
    double a0 = atof(argv[1]);
    double b0 = atof(argv[2]);
    double theta_i = atof(argv[3])*M_PI/180;
    double a = atof(argv[4]);
    std::string technique =  argv[5];
    
    if (abs(a) > 1){// Unphysical black hole spin parameter.
        throw runtime_error("Unphysical spin value");
        return -1;
    }
    
    if (technique != "Simple" && technique != "Disc" && technique != "NoDisc"){
        throw runtime_error("Photon logic 'technique' must be one of: Simple, Disc, NoDisc");
        return -1;
    }
    
    // These parameters are unchanging for the Photon orbit
    // d = observer distance (approximation as in reality -> infinity)
    // F involved in time step in algorithm
    // n is a maximum number of steps per integration per photon, exists to avoid algorithm getting stuck. (Happens when the photon reaches the event horizon).
    
    double d = 10000;
    double F = 1.0/128.0;
    const int n = 50000;
    
    // these will be updated after photon integration with the
    // points in the disc where the photon originated from
    
    double rf = 0;
    double phif = 0;
    double thetaf = 0;
        
    onp.one_orbit(a,a0,b0,theta_i,d,F,n,rf,phif,thetaf,technique); // solve the orbit equations for a single photon
    
    std::cout << rf << std::endl;
    std::cout << thetaf << std::endl;
    std::cout << phif << std::endl;
    return 0;
};
