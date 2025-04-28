#pragma once

/*
TO DO:
- create class for fluid profile to  write to disk, streamplot, print (maybe just have this in vector class?)
*/

double dxi_dtau(double xi, double v, double csq);

double dv_dtau(double xi, double v, double csq);

double dw_dtau(double xi, double v, double w, double csq);

RK4::State profile(RK4::State& init, double t0, double tf, int n);