//in this file we write all th classes and the methods used
#include "Planetesimal.h"
#include <cmath>
#include <iostream>
#include <math.h>
using namespace std;
double const NewtG = 6.67259e-8;
double const au = 1.496*pow(10,13);


Planetesimal::Planetesimal() {

	m_mass = 0;
	m_vx = 0;
	m_vy = 0;
	m_v = 0;
	m_x = 0;
	m_y = 0;
	m_r = 0;

}

Planetesimal::Planetesimal(double mass,double vx, double vy,double x,double y) {

	m_mass = mass;
	m_vx = vx;
	m_vy = vy;
	m_v = sqrt(vx*vx+vy*vy);
	m_x = x;
	m_y = y;
	m_r = sqrt(x*x+y*y);

}


void Planetesimal::setmass(double mass){
	m_mass=mass;
}

void Planetesimal::setvx(double vx){
	m_vx=vx;
}
void Planetesimal::setvy(double vy){
	m_vy=vy;
	m_v = sqrt(m_vx*m_vx+m_vy*m_vy);
	
}
void Planetesimal::setx(double x){
	m_x=x;
}
void Planetesimal::sety(double y){
	m_y=y;
	m_r = sqrt(m_x*m_x+m_y*m_y);

}

void Planetesimal::setr(){
	m_r=sqrt(pow(m_x,2)+pow(m_y,2));
}
double Planetesimal::getx(){
	return m_x;
}
double Planetesimal::gety(){
	return m_y;
}
double Planetesimal::getr(){
	return m_r;
}
void Planetesimal::setv(){
	m_v=sqrt(pow(m_vx,2)+pow(m_vy,2));
}
double Planetesimal::getv(){
	return m_v;
}
double Planetesimal::getvx(){
	return m_vx;
}
double Planetesimal::getvy(){
	return m_vy;
}
double Planetesimal::getmass(){
	return m_mass;
}

double Energy(Planetesimal Halley,Sun sole,Planet jupiter){
	return 0.5*(pow(Halley.getvx(),2)+pow(Halley.getvy(),2))-sole.getmass()*NewtG/Halley.getr()-jupiter.getmass()*NewtG/Distance(Halley,jupiter);
}

Runge_Kutta::Runge_Kutta(){
	m_t=0;
	m_dx=0;
	m_dy=0;
	m_dvx=0;
	m_dvy=0;
}
Runge_Kutta::Runge_Kutta(double t, double x,double y,double vx, double vy){
	m_t=t;
	m_dx=0;
	m_dy=0;
	m_dvx=0;
	m_dvy=0;
}
double Runge_Kutta::Planetesimal_Gravitationalplusdrag_Equation_xaxe(Planetesimal Halley, Sun sole, Planet jupiter, double x,double y)
{
	return -NewtG*(sole.getmass()*x/pow(x*x+y*y,1.5)+jupiter.getmass()*(x-jupiter.getx())/pow(Distance(Halley,jupiter),3));
}
void Runge_Kutta::sett(double t)
{
	m_t=t;
}
double Runge_Kutta::Planetesimal_Gravitationalplusdrag_Equation_yaxe(Planetesimal Halley, Sun sole,  Planet jupiter,double x, double y)
{
	return -NewtG*(sole.getmass()*y/pow(x*x+y*y,1.5)+jupiter.getmass()*(y-jupiter.gety())/pow(Distance(Halley,jupiter),3));
}
void Runge_Kutta::Solve(Planetesimal Halley,Sun sole,  Planet jupiter,double precision,int& counter){
	double x0,x1,x2,x3,x4=0;
	double vx0,vx1,vx2,vx3,vx4=0;
	double y0,y1,y2,y3,y4=0;
	double vy0,vy1,vy2,vy3,vy4=0;
	double dxdouble,dydouble, dvxdouble, dvydouble;
	double dx2,dy2,dvx2,dvy2;
	double dx1,dy1,dvx1,dvy1;
	double dxnormal,dynormal,dvxnormal,dvynormal;
	counter=0;
//	do
//	{
	x0=Halley.getx();
	vx0=Halley.getvx();
	vy0=Halley.getvy();
	y0=Halley.gety();

	vx1=m_t*Planetesimal_Gravitationalplusdrag_Equation_xaxe(Halley,sole, jupiter,x0,y0);
        vy1=m_t*Planetesimal_Gravitationalplusdrag_Equation_yaxe(Halley,sole, jupiter,x0,y0);
	x1=m_t*vx0;
	y1=m_t*vy0;

	vx2=m_t*Planetesimal_Gravitationalplusdrag_Equation_xaxe(Halley,sole, jupiter,x0+x1/2.,y0+y1/2.);
        vy2=m_t*Planetesimal_Gravitationalplusdrag_Equation_yaxe(Halley,sole, jupiter,x0+x1/2.,y0+y1/2.);
     	x2=m_t*(vx0+vx1/2.);
        y2=m_t*(vy0+vy1/2.);

        vx3=m_t*Planetesimal_Gravitationalplusdrag_Equation_xaxe(Halley,sole,jupiter,x0+x2/2.,y0+y2/2.);
        vy3=m_t*Planetesimal_Gravitationalplusdrag_Equation_yaxe(Halley,sole, jupiter,x0+x2/2.,y0+y2/2.);
	x3=m_t*(vx0+vx2/2.);
        y3=m_t*(vy0+vy2/2.);

  	vx4=m_t*Planetesimal_Gravitationalplusdrag_Equation_xaxe(Halley,sole, jupiter,x0+x3,y0+y3);
        vy4=m_t*Planetesimal_Gravitationalplusdrag_Equation_yaxe(Halley,sole, jupiter,x0+x3,y0+y3);
 	x4=m_t*(vx0+vx3);
	y4=m_t*(vy0+vy3);

	dxdouble=(x1+2.*x2+2.*x3+x4)/6;
       	dydouble=(y1+2.*y2+2.*y3+y4)/6;
  	dvxdouble=(vx1+2.*vx2+2.*vx3+vx4)/6;
       	dvydouble=(vy1+2.*vy2+2.*vy3+vy4)/6;


/*
	//ora ho rifaccio con la metà del passo

	vx1=m_t/2.*Planetesimal_Gravitationalplusdrag_Equation_xaxe(Halley,sole, jupiter,x0,y0);
        vy1=m_t/2.*Planetesimal_Gravitationalplusdrag_Equation_yaxe(Halley,sole, jupiter,x0,y0);
	x1=m_t/2.*vx0;
	y1=m_t/2.*vy0;

	vx2=m_t/2.*Planetesimal_Gravitationalplusdrag_Equation_xaxe(Halley,sole, jupiter,x0+x1/2.,y0+y1/2.);
        vy2=m_t/2.*Planetesimal_Gravitationalplusdrag_Equation_yaxe(Halley,sole, jupiter,x0+x1/2.,y0+y1/2.);
     	x2=m_t/2.*(vx0+vx1/2.);
        y2=m_t/2.*(vy0+vy1/2.);

        vx3=m_t/2.*Planetesimal_Gravitationalplusdrag_Equation_xaxe(Halley,sole, jupiter,x0+x2/2.,y0+y2/2.);
        vy3=m_t/2.*Planetesimal_Gravitationalplusdrag_Equation_yaxe(Halley,sole, jupiter,x0+x2/2.,y0+y2/2.);
	x3=m_t/2.*(vx0+vx2/2.);
        y3=m_t/2.*(vy0+vy2/2.);
  	vx4=m_t/2.*Planetesimal_Gravitationalplusdrag_Equation_xaxe(Halley,sole, jupiter,x0+x3,y0+y3);
        vy4=m_t/2.*Planetesimal_Gravitationalplusdrag_Equation_yaxe(Halley,sole, jupiter,x0+x3,y0+y3);
 	x4=m_t/2.*(vx0+vx3);
	y4=m_t/2.*(vy0+vy3);
	dx1=(x1+2.*x2+2.*x3+x4)/6;
       	dy1=(y1+2.*y2+2.*y3+y4)/6;
  	dvx1=(vx1+2.*vx2+2.*vx3+vx4)/6;
       	dvy1=(vy1+2.*vy2+2.*vy3+vy4)/6;
	x0=x0+dx1;
	y0=y0+dy1;
	vx0=vx0+dvx1;
	vy0=vy0+dvy1;

	//ecco fatto primo
	//mo faccio il secondo
	vx1=m_t/2.*Planetesimal_Gravitationalplusdrag_Equation_xaxe(Halley,sole, jupiter,x0,y0);
        vy1=m_t/2.*Planetesimal_Gravitationalplusdrag_Equation_yaxe(Halley,sole, jupiter,x0,y0);
	x1=m_t/2.*vx0;
	y1=m_t/2.*vy0;

	vx2=m_t/2.*Planetesimal_Gravitationalplusdrag_Equation_xaxe(Halley,sole, jupiter,x0+x1/2.,y0+y1/2.);
        vy2=m_t/2.*Planetesimal_Gravitationalplusdrag_Equation_yaxe(Halley,sole, jupiter,x0+x1/2.,y0+y1/2.);
     	x2=m_t/2.*(vx0+vx1/2.);
        y2=m_t/2.*(vy0+vy1/2.);

        vx3=m_t/2.*Planetesimal_Gravitationalplusdrag_Equation_xaxe(Halley,sole, jupiter,x0+x2/2.,y0+y2/2.);
        vy3=m_t/2.*Planetesimal_Gravitationalplusdrag_Equation_yaxe(Halley,sole, jupiter,x0+x2/2.,y0+y2/2.);
	x3=m_t/2.*(vx0+vx2/2.);
        y3=m_t/2.*(vy0+vy2/2.);
  	vx4=m_t/2.*Planetesimal_Gravitationalplusdrag_Equation_xaxe(Halley,sole, jupiter,x0+x3,y0+y3);
        vy4=m_t/2.*Planetesimal_Gravitationalplusdrag_Equation_yaxe(Halley,sole, jupiter,x0+x3,y0+y3);
 	x4=m_t/2.*(vx0+vx3);
	y4=m_t/2.*(vy0+vy3);
	dx2=(x1+2.*x2+2.*x3+x4)/6;
       	dy2=(y1+2.*y2+2.*y3+y4)/6;
  	dvx2=(vx1+2.*vx2+2.*vx3+vx4)/6;
       	dvy2=(vy1+2.*vy2+2.*vy3+vy4)/6;

	counter++;
	dxnormal=dx1+dx2;
	dynormal=dy1+dy2;
	
	if(abs(sqrt(dxdouble*dxdouble+dydouble*dydouble)/sqrt(dxnormal*dxnormal+dynormal*dynormal))/15.>precision)
	{
		m_t=m_t/2.;
	}

	} while(abs(sqrt(dxdouble*dxdouble+dydouble*dydouble)-sqrt(dxnormal*dxnormal+dynormal*dynormal))/15.>precision);

*/
	m_dx=dxdouble;
       	m_dy=dydouble;
  	m_dvx=dvxdouble;
       	m_dvy=dvydouble;
		   

}

double Runge_Kutta::getdx(){
	return m_dx;
}
double Runge_Kutta::getdy(){
	return m_dy;
}
double Runge_Kutta::getdvx(){
	return m_dvx;
}
double Runge_Kutta::getdvy(){
	return m_dvy;
}
double Runge_Kutta::getds(){
	return sqrt(pow(m_dx,2)+pow(m_dy,2));
}
double Runge_Kutta::getdv(){
	return sqrt(pow(m_dvx,2)+pow(m_dvy,2));
}
double Runge_Kutta::gett(){
	return m_t;
}

double Time_Step(Planetesimal Halley){

	return 100;
}

double Planetesimal::angolar_momentum(){
	return m_x*m_vy-m_y*m_vx;
}


Sun::Sun() {
	m_mass = 0;
	m_x = 0;
	m_y = 0;
}

Sun::Sun(double mass,double x,double y) {
	m_mass = mass;
	m_x = x;
	m_y = y;
}
double Sun::getx(){
	return m_x;
}
double Sun::gety(){
	return m_y;
}
double Sun::getmass(){
	return m_mass;
}

void Sun::setx(double x){
	m_x = x;
}
void Sun::sety(double y){
	m_y = y;
}
void Sun::setmass(double mass){
	m_mass=mass;
}







Planet::Planet() {
	m_mass = 0;
	m_x = 0;
	m_y = 0;
	m_omega = 0;
	m_theta = 0;
	m_a = 0;
}

Planet::Planet(double mass,double omega,double theta,double a) {
	m_mass = mass;
	m_a = a;
	m_x = a*cos(theta);
	m_y = a*sin(theta);
	m_omega = omega;
	m_theta = theta;
}
double Planet::getx(){
	return m_x;
}
double Planet::gety(){
	return m_y;
}
double Planet::getmass(){
	return m_mass;
}
double Planet::getr(){
	return sqrt(m_x*m_x+m_y*m_y);
}
double Planet::getomega(){
	return m_omega;
}
double Planet::gettheta(){
	return m_theta;
}
void Planet::setx(){
	m_x = m_a*cos(m_theta);
}
void Planet::sety(){
		m_y = m_a*sin(m_theta);
}
void Planet::setmass(double mass){
	mass = m_mass;
}
void Planet::setomega(double omega){
	m_omega = omega;
}
void Planet::settheta(double theta){
	m_theta = theta;
}
void Planet::uploadtheta(double t){
	m_theta = m_theta+m_omega*t;
	m_x = m_a*cos(m_theta);
	m_y = m_a*sin(m_theta);
}
double Distance(Planetesimal Halley,Planet jupiter){
	double xdist,ydist;
	xdist = Halley.getx()-jupiter.getx();
	ydist = Halley.gety()-jupiter.gety();
	return sqrt(xdist*xdist+ydist*ydist);
	
}

double fRand(double fMin, double fMax)
{
    double f = (double)rand() / RAND_MAX;
    return fMin + f * (fMax - fMin);
}