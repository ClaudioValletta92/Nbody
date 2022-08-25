
#ifndef __planetesimal_h__
#define __planetesimal_h__
//Definizions of all the classes and the functions
;class Sun{
	public:
	Sun();
	Sun(double mass,double x,double y);
	double getx();
	double gety();
	double getmass();
	void setx(double x);
	void sety(double y);
	void setmass(double mass);

	private:
	double m_mass,m_x,m_y;

};

;class Planet{
	public:
	Planet();
	Planet(double mass,double omega,double theta,double a);
	double getx();
	double gety();
	double getmass();
	double getr();
	double getomega();
	double gettheta();
	void setx();
	void sety();
	void setmass(double mass);
	void setomega(double omega);
	void settheta(double theta);
	void uploadtheta(double t);
    

	private:
	double m_mass,m_x,m_y,m_omega,m_theta,m_a;

};

class Planetesimal {
public:
	//costruttore di default
	Planetesimal();
/* mean molecular weight xmu,planetesimal density pden,material strength strength, energy of vaporization E0, heat of fusion particle Ef, size radius of planetesimal  */

	Planetesimal(double mass,double vx,double vy,double x,double y); //costruttore completo
	//singoli costruttori e metodi per ottenere il valore delle variabili
	void setmass(double size);
	void setvx(double vx);
	void setvy(double vy);
	void setv();
	void setx(double x);
	void sety(double y);
	void setr();
	void print();

	double getmass();
	double getvx();
	double getvy();
	double getv();
	double getx();
	double gety();
	double getr();

	double angolar_momentum();


private:
	double m_mass,m_vx, m_vy, m_v, m_x, m_y, m_r;
};



class Runge_Kutta{
	public:
	Runge_Kutta();
	Runge_Kutta(double t,double x,double y,double vx,double vy);
	double Planetesimal_Gravitationalplusdrag_Equation_xaxe(Planetesimal Halley,Sun sole, Planet jupiter,double x, double y);
	double Planetesimal_Gravitationalplusdrag_Equation_yaxe(Planetesimal Halley,Sun sole, Planet jupiter,double x, double y);
	void Solve(Planetesimal Halley,Sun sole,   Planet jupiter,double precision,int& counter);
	double getdx();
	double getdy();
	double getdvx();
	double getdvy();
	double getds();
	double getdv();
	double gett();
	void sett(double t);
	private:
	double m_t, m_dx, m_dy, m_dvx, m_dvy;

};



double fRand(double fMin, double fMax);
double Distance(Planetesimal Halley,Planet jupiter);
double Time_Step(Planetesimal Halley);
double Energy(Planetesimal Halley, Sun sole,Planet jupiter);

#endif
