using namespace std							//usual declaration in c++
#include "Planetesimal.h"						//file with the definition of the classes and functions
#include <iostream>
#include <fstream>
#include <iomanip>

#include <cmath>
#include <stdlib.h>
#include <string>

int main (int argn, char *argv[])
{
	double newx,newy,newvx,newvy;
	double timestep;
	double randau,randtheta;
	int plannumber = 1;
	ofstream out;
	out.open("example.txt");
	int counter;
	string atmostring,outfilestring,outp;
	double energy,initialenergy,initialv;
	double const au = 1.496*pow(10,13);
	double const sunmass = 1.989*pow(10,33);
	double const sunr = 6.957*pow(10,10);
	double const secinyear = 3.154e+7;
	double const earthorbitalvelocity = 29.78*pow(10,5);
	double const secinday =86400;
	double time = 0;
	stringstream ss;

	double const jupomega = 2*M_PI/(secinyear*11.86);
	Sun sole(sunmass,0,0);
	Planet jupiter(1.97*pow(10,30),jupomega,0,5.2*au);
	Runge_Kutta SolveEqMotion;
	Planetesimal Halley2(1,0.01,-3000000,1.002*au,0);
	Planetesimal Halley[plannumber]= {Halley2};
	initialenergy = 0;

	for(int i = 0; i<plannumber; i++){
		do{
			randau = fRand(1,10);
		}while(randau > 6.);
		randtheta = fRand(0,2*M_PI);
		initialv = earthorbitalvelocity/sqrt(randau);
		Halley[i].setmass(1);
		Halley[i].setvx(initialv*sin(randtheta));
		Halley[i].setvy(initialv*cos(randtheta));
		Halley[i].setx(cos(randtheta)*randau*au);
		Halley[i].sety(sin(randtheta)*randau*au);

		Halley[i].setvx(0.1);
		Halley[i].setvy(-3000000);
		Halley[i].setx(1*au);
		Halley[i].sety(0.1);
		
		initialenergy = initialenergy+Energy(Halley[i],sole,jupiter);
	}



	do
	{
	
	
	timestep = Time_Step(Halley[0]);
	time = time +timestep;
	jupiter.uploadtheta(timestep);
	SolveEqMotion.sett(timestep);
	for(int i = 0; i<plannumber; i++){
		SolveEqMotion.Solve(Halley[i],sole,jupiter,0.01,counter);
		newx=Halley[i].getx()+SolveEqMotion.getdx();
		newvx=Halley[i].getvx()+SolveEqMotion.getdvx();
		newy=Halley[i].gety()+SolveEqMotion.getdy();
		newvy=Halley[i].getvy()+SolveEqMotion.getdvy();
		Halley[i].setx(newx);
		Halley[i].sety(newy);
		Halley[i].setvx(newvx);
		Halley[i].setvy(newvy);

		Halley[i].setr();
		Halley[i].setv();

	if(Halley[i].getr()<sunr){
		break;
	}
	if(isnan(Halley[i].getr())){
		break;
	}
	//	cout<<Halley[i].getx()<<"   "<<Halley[i].gety()<<endl;

	
	//cout<<Halley[i].getr()/au<<endl;
	}
	//cout<<time/secinday<<endl;

	if(int(time/secinday) % 10 == 0){
		stringstream ss;
		ss << int(0.1*time/secinday);
		atmostring = ss.str();
		outfilestring = "snapshot" + atmostring +".data";
		outp=outfilestring;

		out.open(outp.c_str(),ios::out);
	}
	
	
	
	energy = 0;
	for(int i = 0; i<plannumber; i++){
		energy = energy+Energy(Halley[i],sole,jupiter);
		if(int(time/secinday) % 10 == 0){
			out<<Halley[i].getx()/au<<"  "<<Halley[i].gety()/au<<"  "<<jupiter.getx()/au<<"  "<<jupiter.gety()/au<<endl;
	}
	}
	out.close();
cout<<Halley[0].getr()/au<<endl;
cout<<(initialenergy-energy)/initialenergy<<"  "<<time/secinday<<endl;


}while(time/secinday <1e4);

return 0;
}
