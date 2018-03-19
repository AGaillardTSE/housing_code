/* **************************************************** */
//      REPLICATION SOMMER ET SULLIVAN (AER, 2017)      //
/*                  February 2018                       */
//               Alexandre Gaillard                     //
/* **************************************************** */


// Description of the algorithm:: use of a DC-EGM method //
// difficulty: because the borrowing constraint depends on h and h', then should track h and h' always (relax this can increase the speed)
// for simplicity: assume a pure static environment, such that borrowing constraint is simpler (prices are constant).

// INPUT: VF (next period), eVF (next period), deVF (next period) [(k,h,y)] and COHendo [(k,h',y,s)], COHendo_min [(h',y,s)] //
// 1. compute COHendo and COHendo_min (for given decision of h',s and k') DONE
// 1bis. compute Ktilde and VFtilde on a grid of COH [(c,h',y,s)] DONE
// 2. Given COHendo and COHendo_min, find all the solutions of COH[c,h'] that correspond to deVF (NOT DONE)
// 3. Test all the solution using the Upper Envelope condition -> keep the best solution (NOT DONE)
// 4. Choice the best value function [save s, hp], and compute VF (this is the original value function) (DONE)
// 5. Smooth the value function VF:: when compute VF = max{next period h / s}, then compute a smoothed VF = Exp{...) //  [this removes the secondary kinks]
// 6. Compute eVF, deVF and go back step 1.


void POLICY(double *VF, double *eVF, double *deVF, double *COHendo, double *COHendo_min, double *Knext, int *Hnow, int *Snow){


/* INITIALIZATION OF VARIABLES / POINTERS */

int y, k, h, s, hp, cc, c, ynext, Kcase, kk, threshold_k, iter, maxS, kkk;

double vmax, coh_temp, K_temp, VF_temp, tempvf, weight, slope1, slope2, slope, cons, klim, critereVF, extrapole, remove;

double *Ktilde, *VFtilde, *VFnew, *PHI, *VFmax;
VFnew = (double *) calloc((dimKHY), sizeof(double));    // value function
VFmax = (double *) calloc((dimKHY), sizeof(double));
PHI = (double *) calloc((dimKHY), sizeof(double));      // smoothed version of value function

// TEMP VALUE, CORRESPONDING TO THE COHENDO VALUE //
VFtilde = (double *) calloc((dimCHYS), sizeof(double)); // value function on the grid of COH
Ktilde = (double *) calloc((dimCHYS), sizeof(double));  // saving function on the grid of COH (or consumption)


critereVF=1.0;
iter=0;


while ((critereVF > epsilonValue)  && (iter < maxiterVF))
{
    
critereVF = 0.0;


/** ||-- EGM STEP --|| Find the COH corresponding to the choice of K', for given choices (s/h') and state (y) **/ // (ok)
// AT this stage, COHendo should be an increasing grid in k (monotonicity of cash on hand) [some value of "k'" are inconsistent with a COH]


#if OMP == 1
omp_set_dynamic(0);     // Explicitly disable dynamic teams
omp_set_num_threads(nbthread); // UiE 4 threads for all consecutive parallel regions

#endif


#if OMP == 1

#pragma omp parallel
{
#pragma omp for private(y,h,maxS,s,threshold_k,c, tempvf, kk, Kcase, remove, cons, weight, slope1, slope2)

#endif

for(y=0;y<length_y;y++){ // current y
    for(h=0; h<length_h; h++){ // choice of h
    
        // can not consume more service than house available //
        // restrict s <= h' if h > 0
        if(h == 0){maxS = length_s;} // case where h=0, car choose any s
        if(h > 0 && h < 6){maxS = h + 1;} // case  where Sgrid is identical to hgrid
        if(h >= 6){maxS = length_s;} // have so much housing capital, so can be landlord.

        for(s=0; s<maxS; s++){ // current decision of s
        
            threshold_k = 0;
            
            
            for(c = 0; c < length_c; c++){ // coh grid (can take into account everything)
            
                if(COH[inxCH(c,h)] <= COHendo_min[inxHYS(h,y,s)]){   // case where you eat the borrowing constraint
                    Ktilde[inxCHYS(c,h,y,s)] = K[inxKH(0,h)];
                    tempvf = eVF[inxKHY(0,h,y)];    // should be always a smoothed function.
                }
                
                if(COH[inxCH(c,h)] > COHendo_min[inxHYS(h,y,s)]){    // case where you don't eat the borrowing constraint
                    
                    kk = max(threshold_k,0); // because richer means save more.
                    Kcase = 0.0;
                    remove = 0.0;
                
                    // first set of search
                    while((remove == 0) && (COH[inxCH(c,h)] > COHendo[inxKHYS(kk,h,y,s)])){
                        if(kk == (length_k-1)){Kcase = 2; remove = 1;}
                        kk++;
                    }
                    
                    if(kk == 0){printf("IMPOSSIBLE");getchar();}
                    
                    // case where saving is above the max value --> extrapolate the point
                    if(Kcase == 2){
                        slope1 = (eVF[inxKHY((length_k-1),h,y)] - eVF[inxKHY((length_k-3),h,y)])/(COHendo[inxKHYS((length_k-1),h,y,s)] - COHendo[inxKHYS((length_k-2),h,y,s)]);
                        slope2 = (K[inxKH((length_k-1),h)] - K[inxKH((length_k-2),h)])/(COHendo[inxKHYS((length_k-1),h,y,s)] - COHendo[inxKHYS((length_k-2),h,y,s)]);
 
                        Ktilde[inxCHYS(c,h,y,s)] = (COH[inxCH(c,h)] - COHendo[inxKHYS((length_k-1),h,y,s)])*slope2/3 + K[inxKH((length_k-1),h)];
                        tempvf = (COH[inxCH(c,h)] - COHendo[inxKHYS((length_k-1),h,y,s)])*slope1 + eVF[inxKHY((length_k-1),h,y)];
                    }

                    // case where saving is below the max value --> interpolate the point
                    if(Kcase == 0){
                        weight = (COH[inxCH(c,h)] - COHendo[inxKHYS((kk-1),h,y,s)])/(COHendo[inxKHYS(kk,h,y,s)] - COHendo[inxKHYS((kk-1),h,y,s)]);
                        
                        Ktilde[inxCHYS(c,h,y,s)] = inter1d(weight,K[inxKH((kk-1),h)],K[inxKH(kk,h)]);
                        tempvf = inter1d(weight,eVF[inxKHY((kk-1),h,y)],eVF[inxKHY(kk,h,y)]);
                    }
                    
                    threshold_k = kk-1;
                    
                }
                
                cons = COH[inxCH(c,h)] - Ktilde[inxCHYS(c,h,y,s)];
                if(cons <= 0.0){printf("consumption negative: %f %f %f cohendo=%f slope=%f gap=%f %f %f %f %f %f %f %d %d %d", cons, Ktilde[inxCHYS(c,h,y,s)], COH[inxCH(c,h)], COHendo[inxKHYS((length_k-1),h,y,s)], slope2, (COH[inxCH(c,h)] - COHendo[inxKHYS((length_k-1),h,y,s)]), K[inxKH((length_k-1),h)] - K[inxKH((length_k-2),h)], (COHendo[inxKHYS((length_k-1),h,y,s)] - COHendo[inxKHYS((length_k-2),h,y,s)]), COHendo[inxKHYS((k-1),h,y,s)], COHendo[inxKHYS((k-2),h,y,s)], COHendo[inxKHYS((k-3),h,y,s)],  COHendo_min[inxHYS(h,y,s)], c, kk, Kcase); getchar();}
                VFtilde[inxCHYS(c,h,y,s)] = U((cons),S[s]) + tempvf;

                
            } // end c
        } // end s
    } // end h
} // end y


#if OMP == 1
#pragma omp for private(y, k, h, vmax, hp, s, maxS, coh_temp, K_temp, VF_temp, slope1, slope2, weight, cc)
#endif

/** Find the best choice of S and H **/
for(y=0;y<length_y;y++){ // current y
    for(k=0;k<length_k;k++){ // previous k
        for(h=0;h<length_h;h++){ // previous h
        
            // DISCRETE CHOICE //
            vmax = -10000;  // set the best value to -10000 first.

            //PHI[inxKHY(k,h,y)] = 0.0; // to smooth the function
            
            for(hp=0;hp<length_h;hp++){ // current h
            
                // restrict s <= h' if h > 0
                if(hp == 0){maxS = length_s;} // case where h=0, car choose any s
                if(hp > 0 && hp < 6){maxS = hp + 1;} // case  where Sgrid is identical to hgrid
                if(hp >= 6){maxS = length_s;} // have so much housing capital, so can be landlord.

                for(s=0;s<maxS;s++){ // current s
                
                
                    coh_temp = COHfun(k,y,h,hp,s); // give the coh value to compute the optimal saving decision.
                    
                    
                    if(coh_temp < cmin[hp]){ // in that case, can not consume positive amount, so not possible
                        K_temp = 0;
                        VF_temp = -10000;
                    }
                    
                    
                    if(coh_temp >= cmin[hp]){
                    
                        // case where you extrapolate the cash on hand because outside the grid //
                        if(coh_temp > COH[inxCH((length_c-1),hp)]){
                        
                            slope1 = (VFtilde[inxCHYS((length_c-1),hp,y,s)] - VFtilde[inxCHYS((length_c-2),hp,y,s)])/(COH[inxCH((length_c-1),hp)] - COH[inxCH((length_c-2),hp)]);
                            slope2 = (Ktilde[inxCHYS((length_c-1),hp,y,s)] - Ktilde[inxCHYS((length_c-2),hp,y,s)])/(COH[inxCH((length_c-1),hp)] - COH[inxCH((length_c-2),hp)]);
 
                            K_temp = (coh_temp - COH[inxCH((length_c-1),hp)])*slope2 + Ktilde[inxCHYS((length_c-1),hp,y,s)];
                            VF_temp = (coh_temp - COH[inxCH((length_c-1),hp)])*slope1 + VFtilde[inxCHYS((length_c-1),hp,y,s)];
                            
                        }
                        
                        
                        if(coh_temp <= COH[inxCH((length_c-1),hp)]){
                        
                            cc = (int)(floor(invlinspace(cmin[hp],cmax,length_c,coh_temp)));
                            
                            // verif, but should not happen, due to step 0. and because coh_temp < COH[inxCH((length_c-1),hp)]
                            if(cc >= (length_c-1)){printf("error cc > cmax");getchar();}
                            if(cc < 0){printf("error cc < 0");getchar();}
                            
                            weight = (coh_temp - COH[inxCH(cc,hp)])/(COH[inxCH((cc+1),hp)] - COH[inxCH(cc,hp)]);
                            
                            VF_temp = inter1d(weight,VFtilde[inxCHYS(cc,hp,y,s)],VFtilde[inxCHYS((cc+1),hp,y,s)]);
                            K_temp = inter1d(weight,Ktilde[inxCHYS(cc,hp,y,s)],Ktilde[inxCHYS((cc+1),hp,y,s)]);
                        }

                    }
                    
                    
                    if(VF_temp > vmax){
                        VFmax[inxKHY(k,h,y)] = VF_temp;
                        Snow[inxKHY(k,h,y)] = s;
                        Hnow[inxKHY(k,h,y)] = hp;
                        Knext[inxKHY(k,h,y)] = K_temp;
                        vmax = VF_temp;
                    }

                    //cons_temp = max(coh_temp - K_temp,0.0000001);
                    //cons[inxKHYHS(k,h,y,hp,s)] = cons_temp;
                    //MUnext[inxKHYHS(k,h,y,hp,s)] = MUc(cons[inxKHYHS(k,h,y,hp,s)],s);
                    //LOGIT[inxKHYHS(k,h,y,hp,s)] = exp(VF_temp/smoothpar);
//                    if(VF_temp > -500){
//                        PHI[inxKHY(k,h,y)] += exp(VF_temp/smoothpar); // compute PHI for smoothing.
//                    }
                    //printf("SOLUTION: %d %d %d %d %d %f %f %f %f %f %20.40f", h,y,k,hp,s, K_temp, VF_temp, coh_temp, cmin[hp], COH[inxCH((length_c-1),hp)], PHI[inxKHY(k,h,y)]);getchar();
                    
                } // end s
            } // end hp
            
            // 4.1 COMPUTE THE SMOOTHED VALUE FUNCTION //
            //if(smoothpar*log(PHI[inxKHY(k,h,y)]) < -1000000){VFnew[inxKHY(k,h,y)] = VFmax[inxKHY(k,h,y)];}else{VFnew[inxKHY(k,h,y)] = smoothpar*log(PHI[inxKHY(k,h,y)]);}
            
            VFnew[inxKHY(k,h,y)] = VFmax[inxKHY(k,h,y)]; // TEMPORAIRE.
            
            //if(VFnew[inxKHY(k,h,y)] < -1000000){printf("SMOOTHED VF= %f %d %d %d", VFnew[inxKHY(k,h,y)], k , h, y);getchar();}
            
            
        } // end y    (should be OK)
    }
}

#if OMP == 1
}
#endif



/** Find the best choice of S and H **/
for(y = 0; y<length_y; y++){
    for(k=0;k<length_k;k++){ // previous k
        for(h=0;h<length_h;h++){ // previous h
        
            // compute criterion (on VFI) //
            critereVF = max(fabs(VFnew[inxKHY(k,h,y)] - VF[inxKHY(k,h,y)]),critereVF);
            
            VF[inxKHY(k,h,y)] = VFnew[inxKHY(k,h,y)];
            
            // COMPUTE EXPECTATIONS //
            eVF[inxKHY(k,h,y)] = 0.0;
            
            for(ynext = 0; ynext<length_y; ynext++){
                if(y == 6 && ynext == 0){
                    eVF[inxKHY(k,h,y)] += PI_tot[y][ynext]*betapar*0.0;     // take into account the death
                }else{
                    eVF[inxKHY(k,h,y)] += PI_tot[y][ynext]*betapar*VF[inxKHY(k,h,ynext)];
                }
            } // end ynext
        
        
            // COMPUTE DERIVATIVES: with respect to k //
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
            // restrict s <= h' if h > 0 // CAUTION :: here "h" is the end of period house value.
            if(h == 0){maxS = length_s;} // case where h=0, car choose any s
            if(h > 0 && h < 6){maxS = h + 1;} // case  where Sgrid is identical to hgrid
            if(h >= 6){maxS = length_s;} // have so much housing capital, so can be landlord.

            for(s = 0; s<maxS; s++){ // s current
            
                if(k == 1){COHendo_min[inxHYS(h,y,s)] = K[inxKH(0,h)] + invMU(deVF[inxKHY(0,h,y)],S[s]);}
                if(k >= 2){COHendo[inxKHYS((k-1),h,y,s)] = K[inxKH((k-1),h)] + invMU(deVF[inxKHY((k-1),h,y)],S[s]);}  // here "k" is the end of period capital value.
                if(k == (length_k-1)){COHendo[inxKHYS(k,h,y,s)] = K[inxKH(k,h)] + invMU(deVF[inxKHY(k,h,y)],S[s]);}
                if(k == 1){COHendo[inxKHYS(0,h,y,s)] = K[inxKH(0,h)] + invMU(deVF[inxKHY(0,h,y)],S[s]);}
                

            } // end s
        } // end y
    } // end h
} // end k


if(critereVF > epsilonValue){
    UEC(deVF, eVF, COHendo);
}





iter++;

//printf("Rule CNVG %d\t%20.15f\n",iter,critereVF);




}//while (critereVF > epsilonValue)





//FILE *valfileout;
//
//valfileout=fopen(valfile, "w");
//
//setbuf ( valfileout, NULL );
//
//fprintf(valfileout,"%5s\t%5s\t%5s\t%20s\t%20s\t%5s\t%5s\t%20s\t%20s\t%20s\t%20s\n","PROD", "HOUSE","ASSET", "ASSETVAL","K'", "H'", "S", "VF", "eVF", "deVF", "phi");
//
//for(y=0;y<length_y;y++)
//{
//    for(h=0; h<length_h; h++)
//    {
//        for(k=0; k<length_k; k++)
//        {
//            fprintf(valfileout,"%5d\t%5d\t%5d\t%20.15f\t%20.15f\t%5d\t%5d\t%20.15f\t%20.15f\t%20.15f\t%20.15f\n", y, h, k, K[inxKH(k,h)], Knext[inxKHY(k,h,y)],Hnow[inxKHY(k,h,y)],Snow[inxKHY(k,h,y)],VFmax[inxKHY(k,h,y)],eVF[inxKHY(k,h,y)],deVF[inxKHY(k,h,y)], smoothpar*log(PHI[inxKHY(k,h,y)]));
//        }
//    }
//}
//fclose(valfileout);




    


// FREE calloc //
free(VFtilde);
free(Ktilde);
free(VFnew);

}


