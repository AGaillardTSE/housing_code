/* **************************************************** */
//      REPLICATION SOMMER ET SULLIVAN (AER, 2017)      //
/*                  February 2018                       */
//               Alexandre Gaillard                     //
/* **************************************************** */


// POLICY :: remove mortgage reduction in tax //
const double policy_mortgage = 1.0;


// OPTION //
#define OMP 1
int const nbthread = 7;


/****************/
//    INCLUDE   //
/****************/

#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>
#include <math.h>
#include <limits>
#include <assert.h>
#include <cstring>
#include <stdio.h>
#include <fstream>
#include <iostream>
#include <random>
#if OMP == 1 
#include <omp.h>
#endif



// TYPE 1 SMOOTHER //
double smoothpar = 0.5;


#define deriv(val1,val2,val3,x1,x2,x3) ((1.0 - (x3 - x2)/(x3 - x1))*((val3 - val2)/(x3-x2)) + ((x3 - x2)/(x3 - x1))*((val2 - val1)/(x2-x1)))


/******************************/
//  GRID DEFINITION AND INDEX //
/******************************/


// LENGTH OF THE GRIDS //
#define length_k 500 // must be an even number !!
#define length_y 7
#define length_h 7
#define length_s 17
#define length_c 500


// DIMENSION FOR POINTERS //
#define dimKHY (length_k*length_h*length_y)
#define dimHYS (length_h*length_y*length_s)
#define dimKHYS (length_k*length_h*length_y*length_s)
#define dimCHYS (length_c*length_h*length_y*length_s)
#define dimKH (length_k*length_h)
#define dimCH (length_c*length_h)



// INDEX FOR POINTERS //
#define inxKHY(ki,hi,yi) ((length_k)*(length_h)*(yi)+(length_k)*(hi)+(ki))
#define inxHYS(hi,yi,si) ((length_h)*(length_y)*(si)+(length_h)*(yi)+(hi))
#define inxKHYS(ki,hi,yi,si) ((length_k)*(length_h)*(length_y)*(si)+(length_k)*(length_h)*(yi)+(length_k)*(hi)+(ki))
#define inxCHYS(ci,hi,yi,si) ((length_c)*(length_h)*(length_y)*(si)+(length_c)*(length_h)*(yi)+(length_c)*(hi)+(ci))
#define inxKH(ki,hi) ((length_k)*(hi)+(ki))
#define inxCH(ci,hi) ((length_c)*(hi)+(ci))



// DEPOSIT, MORTGAGE //
double kmax = 20;
double cmax = 18;

double K[dimKH], Kprev[dimKH], COH[dimCH], cmin[length_h];


// HOUSING //
double H[length_h] = {0.0, 1.0, 1.25, 1.5, 1.75, 2.0, 5.0};
double S[length_s] = {0.5, 1.0, 1.25, 1.5, 1.75, 2.0, 2.25, 2.5, 3.0, 3.25, 3.5, 3.75, 4.0, 4.25, 4.5, 4.75, 5.0};



// STOCHASTIC AGING ECONOMY //
// Following Heathcote (2005), Castaneda, Diaz-Gimenez and Rios-Rull (2003) //

double Y[length_y] = {0.34250,0.48950,0.69970,1.0,1.42920,2.04270,2.91950};
double PI_tot[length_y][length_y] = {
{0.37390, 0.61670, 0.00940, 0.0, 0.0, 0.0, 0.0},
{0.08690, 0.49170, 0.41310, 0.00830, 0.0, 0.0, 0.0},
{0.00190, 0.12650, 0.56130, 0.30470, 0.00560, 0.0, 0.0},
{0.0, 0.00340, 0.16830, 0.57980, 0.24510, 0.00340, 0.0},
{0.0, 0.0, 0.00560, 0.20920, 0.56120, 0.22210, 0.00190},
{0.0, 0.0, 0.0, 0.00830, 0.2350, 0.49170, 0.2650},
{0.41530, 0.0, 0.0, 0.0, 0.00940, 0.20140, 0.37390}
};






/**********/
// PRICES //
/**********/

double rent = 0.2477384507660;       // rental price implied market clearing
double price = 3.0521774291990;        // housing price implied by market clearing

/**************/
// PARAMETERS //
/**************/

// Households //
const double nper = 0.01;       // 2% of population growth
double betapar = 0.985;         //  discount rate
double sigmaCRRA = 2.5;         // CRRA PARAMETER
double sigmaCOBB = 0.685;       // COBB-doublas parameter
double theta = 0.2;             // Down payment requirement
#define U(x,y) (pow(pow((x),sigmaCOBB)*pow((y),(1-sigmaCOBB)),(1.0-sigmaCRRA))/(1.0-sigmaCRRA))
#define invMU(mu,ss) (pow((mu/(sigmaCOBB*pow(ss,(1-sigmaCOBB)*(1-sigmaCRRA)))),(1/(sigmaCOBB*(1-sigmaCRRA)-1))))

// fixed prices //
double wage = 1.0;              // wage rate implied by market clearing
double irateS = 0.04;           // interest rate implied by market clearing
double irateB = irateS + 0.015;

// Government //
double tau_s = 0.07;            // selling cost
double tau_b = 0.025;           // buying cost
double tau_p = 0.075;           // payroll tax
double tau_h = 0.01;            // property tax
double tau_r = 0.01;            // tax on income rent
double tau_m = 1.00;            // mortgage deductibility rate
double tau_LL = 0.023;          // deductibility rate for mortgage
double epar = 0.09498282502;
double xipar = 0.14832934318;
double threshold[7] = {0.0,0.21728947642,0.88347038617,2.1403663995,4.4641927761,9.70516290205,999999};
double prog_rate[6] = {0.10,0.15,0.25,0.28,0.33,0.35};

// Housing sector //
double deltaH = 0.015;          // maintenance cost rate
double phipar = 0.056;          // fixed cost for landlord






// Linspace functions //
#define linspace2(xmin,xmax,n,i) ((i)>=(100))?((i-100)*(((xmax)-(Gridmiddle))/(n/2))+(Gridmiddle)):(((i)*(((Gridmiddle)-(xmin))/((n/2)))+(xmin)))
#define invlinspace2(xmin,xmax,n,x) ((x)>=(Gridmiddle))?(((((x)-(Gridmiddle))/((xmax)-(Gridmiddle)))*((n/2)))+100):(((((x)-(xmin))/((Gridmiddle)-(xmin)))*((n/2))))















/*************************/
//  NUMERICAL PARAMETERS //
/*************************/

// distibution criterion
#define epsdist 0.0000000001

//value function cnvg criterion
#define epsilonValue 0.000001
const int maxiterVF = 1000;

//aggregate capital cnvg criterion
#define epsilonprice 0.0001
#define itermax 20

//aggregate capital relaxation
#define relaxK 0.2

//aggregate capital relaxation
#define relaxL 0.2

//aggregate tax relaxation
#define relaxTax 0.2

// TOLERANCE MNBRACK
const double GOLD=1.618034,GLIMIT=100.0,TINY=1.0e-20;

// TOLERANCE MNGOLDEN
const double TOL= 1.0e-5;
const double R=0.61803399,C=1.0-R;


// FILE OUT //
const char distfile[]="dist.out";
const char valfile[]="val.out";
const char distassetfile[]="distwealth.out";


