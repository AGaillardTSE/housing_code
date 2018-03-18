/* **************************************************** */
//      REPLICATION SOMMER ET SULLIVAN (AER, 2017)      //
/*                  February 2018                       */
//               Alexandre Gaillard                     //
/* **************************************************** */


// METHOD used:: EGM method as in Rust et al. (2017) //
// we define the Endogenous Grid over a value of Cash-on-hand, this reduce substantially the size of loops //



/* ********* */
/* CPP FILES */
/* ********* */

// HEADER
#include "header.h"

// USEFUL
#include "/Users/alexandregaillard/Documents/Compiler/lib/libperso/useful.cpp"
#include "/Users/alexandregaillard/Documents/Compiler/lib/libperso/tauchen.cpp"

// CPP FILES //
#include "cohfun.cpp"
#include "upper_envelope_condition.cpp"
#include "policy.cpp"
#include "SIMULATION.cpp"






int main(int argc, char* argv[])
{


/********************/
/** INITIALIZATION **/
/********************/

int k, y, h, s, c, maxS;
double coh_temp;

// VALUE AND POLICY FUNCTIONS //
double *VF, *Know, *eVF, *deVF, *MUtilde, *VFtilde;
int *Hnow, *Snow;
deVF = (double *) calloc((dimKHY), sizeof(double));
eVF = (double *) calloc((dimKHY), sizeof(double));
VF = (double *) calloc((dimKHY), sizeof(double));
Know = (double *) calloc((dimKHY), sizeof(double));
Hnow = (int *) calloc((dimKHY), sizeof(int));
Snow = (int *) calloc((dimKHY), sizeof(int));

double *COHendo, *COHendo_min;
COHendo = (double *) calloc((dimKHYS), sizeof(double));
COHendo_min = (double *) calloc((dimHYS), sizeof(double));



// DISTRIBUTION //
double *dist;
dist = (double *) calloc((dimKHY), sizeof(double));
dist[0] = 1.0; // guess



double critereprice, critrent, aggRout, aggHout, currentH, pricemin, pricemax, rentmin, rentmax;
double epsprice = 0.0001;
double relaxP = 0.01;
double relaxR = 0.01;
double epsrent = 0.0001;
double Hsupply = 1.2;


pricemax = 3.06;
pricemin = 3.04;
price = (pricemax + pricemin)/2;


// GUESS //

/** define COH grid and MUtilde and VFtilde **/
for(h = 0; h < length_h; h++){ // chosen h
    for(c = 0; c < length_c; c++){ // cash on hand
        COH[inxCH(c,h)] = linspace(((-(1-theta)*price*H[h]) + 0.00000001),cmax,length_c,c); // the dependence on h to constraint the set of cash on hand to be feasible.
     }
    
     cmin[h] = (-(1-theta)*price*H[h]) + 0.00000001;
}

/** define VF, K, cmin **/
for(h = 0; h < length_h; h++){ // this is the chosen h
    for(k = 0; k < length_k; k++){  // this is the current y

        if(h > 0){K[inxKH(k,h)] = linspace((-(1-theta)*price*H[h]),kmax,length_k,k);}
        if(h == 0){K[inxKH(k,h)] = linspace(0,kmax,length_k,k);}
        
        for(y = 0; y < length_y; y++){ // current shock income
        
            coh_temp = COHfun(k,y,0,0,0); // (first guess)
            
            VF[inxKHY(k,h,y)] = U(coh_temp,S[0])/(1-betapar);  // (has to be increasing in k)
            
        }
    }
}
    

/** 1. copy VF into eVF **/
bascule(VF,eVF,dimKHY);


/** 2. COMPUTE first derivatives + COH_endo **/  // (OK -- first time == no kinks)
for(k=0;k<length_k;k++){ // next k
    for(h=0;h<length_h;h++){ // next h
        for(y=0;y<length_y;y++){ // current y
        
            // 2.1 COMPUTE DERIVATIVES: with respect to k //
            if(k == 1){
                deVF[inxKHY(0,h,y)] = (eVF[inxKHY(1,h,y)]-eVF[inxKHY(0,h,y)])/(K[inxKH(1,h)] - K[inxKH(0,h)]);
            }
            
            if(k >= 2){
                deVF[inxKHY((k-1),h,y)] = deriv(eVF[inxKHY((k-2),h,y)],eVF[inxKHY((k-1),h,y)],eVF[inxKHY(k,h,y)],K[inxKH((k-2),h)],K[inxKH((k-1),h)],K[inxKH(k,h)]);
            }
            
            if(k == (length_k-1)){
                deVF[inxKHY((length_k-1),h,y)] = (eVF[inxKHY((length_k-1),h,y)]-eVF[inxKHY((length_k-2),h,y)])/(K[inxKH((length_k-1),h)] - K[inxKH((length_k-2),h)]);
            }

            // COMPUTE THE COH_ENDO //
            if(h == 0){maxS = length_s;} // case where h=0, car choose any s
            if(h > 0 && h < 6){maxS = h + 1;} // case  where Sgrid is identical to hgrid
            if(h >= 6){maxS = length_s;} // have so much housing capital, so can be landlord.

            for(s = 0; s<maxS; s++){ // s current
            
                if(k == 1){
                    COHendo_min[inxHYS(h,y,s)] = K[inxKH(0,h)] + invMU(deVF[inxKHY(0,h,y)],S[s]);
                    COHendo[inxKHYS(0,h,y,s)] = K[inxKH(0,h)] + invMU(deVF[inxKHY(0,h,y)],S[s]);
                }
                if(k >= 2){
                    COHendo[inxKHYS((k-1),h,y,s)] = K[inxKH((k-1),h)] + invMU(deVF[inxKHY((k-1),h,y)],S[s]);
                    
                    if(COHendo[inxKHYS(k-1,h,y,s)] < COHendo[inxKHYS((k-2),h,y,s)]){printf("BUG COHk < COHk-1: %f %f %d %d %d %d %f %f %f %f ", COHendo[inxKHYS(k-1,h,y,s)], COHendo[inxKHYS((k-2),h,y,s)], k-1, k-2, y, h, K[inxKH(k-1,h)], invMU(deVF[inxKHY((k-1),h,y)],S[s]),K[inxKH(k-2,h)], invMU(deVF[inxKHY(k-2,h,y)],S[s]));getchar();
                    }
                }  // here "k" is the end of period capital value.
                if(k == (length_k-1)){COHendo[inxKHYS(k,h,y,s)] = K[inxKH(k,h)] + invMU(deVF[inxKHY(k,h,y)],S[s]);}
                
                
            } // end s
            
 
        } // end y
    } // end h
} // end k








/** START EQUILIBRIUM FIXED POINT **/
printf("STARTING EQUILIBRIUM FIXED POINT...\n");

timeval t1, t2, t3, t4;
double elapsedTime;


// start timer
gettimeofday(&t1, NULL);




pricemax = 3.06;
pricemin = 3.04;
price = (pricemax + pricemin)/2;
//critereprice = 1.0;
aggHout= 1.0;
// STARTING EQUILIBRIUM ON HOUSING PRICES //
// rental and housing market are too much dependent --> nested loop //
while(fabs(aggHout) > epsprice){
    
    // check time //
    gettimeofday(&t3, NULL);
    
    
    // STARTING EQUILIBRIUM ON RENTAL MARKET //
    rentmax = 0.2425;
    rentmin = 0.2375;
    rent = (rentmax + rentmin)/2;
    
    aggRout = 1.0;
    while(fabs(aggRout) > epsrent){
        
        rent = (rentmax + rentmin)/2;
        
        // POLICY ITERATION (using EGM) //
        POLICY(VF,eVF,deVF,COHendo,COHendo_min,Know,Hnow,Snow);

        // SIMULATION (using histogram method) //
        SIMULATION(Know,Hnow,Snow,dist,&aggRout,&aggHout,&currentH);


//        // CHECK EQUILIBRIUM CONDITIONS FOR RENTAL MARKET //
//        if(aggRout > 0.0001 || aggRout < 0.0001){ // then update rent
//            newrent = rent + rent*(aggRout);
//        }
//
//        // New prices //
//        newrent = relaxR*newrent + (1-relaxR)*rent;
//        newprice = relaxP*newprice + (1-relaxP)*price;
        
        
//        // criterion //
//        critereprice = (abs(newrent - rent)/rent);
//        critereprice = max(critereprice,(abs(newprice - price)/price));
        
        printf("EQUILIBRIUM RENTAL MARKET :: aggRout=%f, (aggHout-Hsupply)=%f,  rent=%20.15f price=%20.15f \n",  aggRout, aggHout - Hsupply,rent , price);
        
        
        // update price //
        if(aggRout > 0.0){ // then rent is too low
            rentmin = rent;
        } else {
            rentmax = rent;
        }
        
        
    } // end looking for rental market clearing condition
    
    
    
    gettimeofday(&t4, NULL);

    elapsedTime = (t4.tv_sec - t3.tv_sec) * 1000.0;
    elapsedTime += (t4.tv_usec - t3.tv_usec) / 1000.0;

    
    // Housing supply //
    //Hsupply = currentH*(nper);

//    // CHECK EQUILIBRIUM CONDITIONS FOR HOUSING MARKET //
//    if(aggHout > Hsupply || aggHout < Hsupply){ // then update rent
//        newprice = price + price*(aggHout - Hsupply);
//    }
//
//    newprice = relaxP*newprice + (1-relaxP)*price;
//
//    // criterion //
//    critereprice = (abs(newrent - rent)/rent);
//    critereprice = max(critereprice,(abs(newprice - price)/price));

    printf("EQUILIBRIUM HOUSING MARKET ::  aggRout=%f, (aggHout-Hsupply)=%f, rent=%f, price=%f, time=%fs\n",  aggRout, aggHout - Hsupply,rent, price, elapsedTime/1000);

    
    if(aggHout > Hsupply){ // demand is too high -> increase price
        pricemin = price;
    } else {
        pricemax = price;
    }
    
    price = (pricemax + pricemin)/2;


    /** REDEFINE THE GRID **/
    
    // COMPUTE NEW COH GRID GIVEN PRICES & NEW GUESS //
    /** define COH grid and MUtilde and VFtilde **/
    for(h = 0; h < length_h; h++){ // chosen h
        for(c = 0; c < length_c; c++){ // cash on hand
            COH[inxCH(c,h)] = linspace(((-(1-theta)*price*H[h]) + 0.00000001),cmax,length_c,c); // the dependence on h to constraint the set of cash on hand to be feasible.
         }
        
         cmin[h] = (-(1-theta)*price*H[h]) + 0.00000001;
    }


    double slope, Kcase, remove, weight;
    int threshold_k, kk;

    /** INTERPOLATE NEW VF ON THE NEW K GRID (Knew is increasing with k == monotonicity) **/
    bascule(K,Kprev,dimKH);

    for(h = 0; h < length_h; h++){
        for(y = 0; y < length_y; y++){
            
            threshold_k = 0;
            
            for(k = 0; k < length_k; k++){
            
                // compute the new grid point.
                if(h > 0){K[inxKH(k,h)] = linspace((-(1-theta)*price*H[h]),kmax,length_k,k);}
                if(h == 0){K[inxKH(k,h)] = linspace(0,kmax,length_k,k);}
        

                kk = max(threshold_k,0.0);
                Kcase = 0.0;
                remove = 0.0;
            
                if(K[inxKH(k,h)] <= Kprev[inxKH(0,h)]){Kcase = 1;remove=1;} // case where the new K grid is before the min previous grid
                
                // search the previous closest point //
                while((remove == 0) && (K[inxKH(k,h)] >= Kprev[inxKH(kk,h)])){
                    if(kk == (length_k-1)){Kcase = 2; remove = 1;}  // case where the new K grid is after the max previous grid
                    kk++;
                }
                
                if(Kcase == 0){ // normal case
                    weight = (K[inxKH(k,h)] - Kprev[inxKH((kk-1),h)])/(Kprev[inxKH(kk,h)] - Kprev[inxKH((kk-1),h)]);
                    eVF[inxKHY(k,h,y)] = inter1d(weight,VF[inxKHY((kk-1),h,y)],VF[inxKHY(kk,h,y)]);
                }
                
                if(Kcase == 2){
                    slope = (VF[inxKHY((length_k-1),h,y)] - VF[inxKHY((length_k-2),h,y)])/(Kprev[inxKH((length_k-1),h)] - Kprev[inxKH((length_k-2),h)]);
                    eVF[inxKHY(k,h,y)] = (K[inxKH(k,h)] - Kprev[inxKH((length_k-1),h)])*slope + VF[inxKHY((length_k-1),h,y)];
                }
                
                if(Kcase == 1){
                    slope = (VF[inxKHY(1,h,y)] - VF[inxKHY(0,h,y)])/(Kprev[inxKH(1,h)] - Kprev[inxKH(0,h)]);
                    eVF[inxKHY(k,h,y)] = (K[inxKH(k,h)] - Kprev[inxKH(0,h)])*slope + VF[inxKHY(0,h,y)];
                }
                
                threshold_k = kk - 1;
            }
        }
    }


    /** 2. COMPUTE first derivatives + COH_endo **/  // (OK -- first time == no kinks)
    for(k=0;k<length_k;k++){ // next k
        for(h=0;h<length_h;h++){ // next h
            for(y=0;y<length_y;y++){ // current y
            
                // 2.1 COMPUTE DERIVATIVES: with respect to k //
                if(k == 1){
                    deVF[inxKHY(0,h,y)] = (eVF[inxKHY(1,h,y)]-eVF[inxKHY(0,h,y)])/(K[inxKH(1,h)] - K[inxKH(0,h)]);
                }
                
                if(k >= 2){
                    deVF[inxKHY((k-1),h,y)] = deriv(eVF[inxKHY((k-2),h,y)],eVF[inxKHY((k-1),h,y)],eVF[inxKHY(k,h,y)],K[inxKH((k-2),h)],K[inxKH((k-1),h)],K[inxKH(k,h)]);
                }
                
                if(k == (length_k-1)){
                    deVF[inxKHY((length_k-1),h,y)] = (eVF[inxKHY((length_k-1),h,y)]-eVF[inxKHY((length_k-2),h,y)])/(K[inxKH((length_k-1),h)] - K[inxKH((length_k-2),h)]);
                }
                
                
                // COMPUTE THE COH_ENDO //
                if(h == 0){maxS = length_s;} // case where h=0, car choose any s
                if(h > 0 && h < 6){maxS = h + 1;} // case  where Sgrid is identical to hgrid
                if(h >= 6){maxS = length_s;} // have so much housing capital, so can be landlord.

                for(s = 0; s<maxS; s++){ // s current
                
                    if(k == 1){
                        COHendo_min[inxHYS(h,y,s)] = K[inxKH(0,h)] + invMU(deVF[inxKHY(0,h,y)],S[s]);
                        COHendo[inxKHYS(0,h,y,s)] = K[inxKH(0,h)] + invMU(deVF[inxKHY(0,h,y)],S[s]);
                    }
                    if(k >= 2){COHendo[inxKHYS((k-1),h,y,s)] = K[inxKH((k-1),h)] + invMU(deVF[inxKHY((k-1),h,y)],S[s]);}  // here "k" is the end of period capital value.
                    if(k == (length_k-1)){COHendo[inxKHYS(k,h,y,s)] = K[inxKH((length_k-1),h)] + invMU(deVF[inxKHY((length_k-1),h,y)],S[s]);}
                    
                } // end s
            } // end y
        } // end h
    } // end k

    UEC(deVF, eVF, COHendo);
    
} // end while equilibrium.




gettimeofday(&t2, NULL);
elapsedTime = (t2.tv_sec - t1.tv_sec) * 1000.0;
elapsedTime += (t2.tv_usec - t1.tv_usec) / 1000.0;

printf("===================Program ended in %f seconds\n",elapsedTime/1000);
    
return 0;
}




