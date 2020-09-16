#include <cstdio>
#include <cmath>
#include "Control.d/ControlInterface.h"

class MyControl : public ControlInterface {

  protected:
    double dt;	       // time step size
    
  public:
  
    // set time step
    void setDt(double h) { dt = h; }

    // Control force initialization routine
    void init(double *displacement, double *velocity, double *acc,
                      SingleDomainDynamic * probDesc=0);

    // Control force routine
    void ctrl(double *dis, double *vel, double *acc, double *f,
                      double time=0, 
		      SysState<Vector> *state=0, Vector *ext_f=0);

    // User defined force routine
    void usd_forc(double time, double *userDefineForc);

    // User defined displacement routine
    void usd_disp(double time, double *userDefineDisp,
                          double *userDefineVel, double *userDefineAcc);
};

ControlInterface *controlObj = new MyControl();

void
MyControl::ctrl(double *displacement, double *velocity, double *acceleration, 
              double *force, double time, SysState<Vector> *state, 
	      Vector *ext_f)
{
 // blank intentionally
}

void
MyControl::init(double *displacement, double *velocity, double *acceleration,
                SingleDomainDynamic * probDesc)
{
 // blank intentionally
}

void
MyControl::usd_disp(double time, double *userDefineDisp, double *userDefineVel, double *userDefineAcc)
{
  const double f = 14.4; //Hz
  const double theta0 = 22.5; //Degrees
  const double theta0_rad = (theta0 * M_PI) / 180.0;
  const double w = 2.0 * M_PI * f;
  const double phi = -0.759955086;
  const double offset = 15.5*M_PI/180.;
  
  //x-aligned Axis of Rotation
  userDefineDisp[0] = (theta0_rad    * sin(w*time+phi) + offset);
  userDefineDisp[1] = -(theta0_rad  * sin(w*time+phi) + offset);
  userDefineVel[0]  = w*theta0_rad  * cos(w*time+phi);
  userDefineVel[1]  = -(w*theta0_rad * cos(w*time+phi));
}

void
MyControl::usd_forc(double time, double *usdForce)
{
 // blank intentionally
}

