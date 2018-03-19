
///***************************************/
///** UPPER ENVELOPE CONDITION (DC-EGM) **/
///***************************************/


void UEC(double *deVF, double *eVF, double *COHendo){

double slope, extrapole, remove, Kcase, slope2;
int y, h, s, maxS, klim, k, kk, kkk;

for(y = 0; y<length_y; y++){
    for(h=0;h<length_h;h++){ // chosen h
    
        // restrict s <= h' if h > 0 // CAUTION :: here "h" is the end of period house value.
        if(h == 0){maxS = length_s;} // case where h=0, car choose any s
        if(h > 0 && h < 6){maxS = h + 1;} // case  where Sgrid is identical to hgrid
        if(h >= 6){maxS = length_s;} // have so much housing capital, so can be landlord.

        for(s = 0; s<maxS; s++){ // s current
            
            klim = 0;

            for(k=0;k<(length_k);k++){ // previous k

                if(k > 0 && k < (length_k-1)){

                if(COHendo[inxKHYS(k,h,y,s)] < COHendo[inxKHYS(k-1,h,y,s)] && k >= klim){ // then means that there is a discontinuity here.

                    //printf("UEC :: %f %f", COHendo[inxKHYS(k,h,y,s)], COHendo[inxKHYS(k-1,h,y,s)]);getchar();
                    
                    // 1. find the point where COHendo(k) < COHendo(k+1)
                    kk = k;

                    // Bound the set J2, where there are several value of COH around the same place.
                    while((COHendo[inxKHYS(kk,h,y,s)] > COHendo[inxKHYS((kk+1),h,y,s)]) && (kk < length_k)){
                        if(kk >= (length_k-2)){Kcase = 2;
                        //printf("stop problem in upper envelope condition :: %f %f %f %f %d %d %d", deVF[inxKHY(k,h,y)], deVF[inxKHY(k-1,h,y)], COHendo[inxKHYS(k,h,y,s)], COHendo[inxKHYS(k-1,h,y,s)], kk, k, kk+1);getchar();
                        
                        }
                        kk++;
                    }

                    // AT THIS STAGE: we get the point at which we recover something increasing (kk), and we have the point at which it falls (k-1)
                    
                    
                    if(Kcase < 2){
                    
                    // NOW TEST THE VALUE of COH AT (k-1) EVALUATED AT THE value after the kink where coh is increasing //
                    slope = (eVF[inxKHY((kk+1),h,y)] - eVF[inxKHY((kk),h,y)])/(COHendo[inxKHYS((kk+1),h,y,s)] - COHendo[inxKHYS((kk),h,y,s)]);
                    extrapole = slope2*(COHendo[inxKHYS((k-1),h,y,s)] - COHendo[inxKHYS(kk,h,y,s)]) + eVF[inxKHY((kk),h,y)];
                    remove = 0;
                    kkk = (k-1);

                    // remove the bad solution values
                    // 1. if extrapole(k-1) > eVF[k-1]:: then remove COH(k-1), since not optimal to choose this. Compute new extrapolation.
                    // 2. if extrapole(k-1) < eVF[k-1]:: then remove COH(kk), and try the next point kk+1. Compute the new extrapolation.
                    while(remove == 0){

                        if(extrapole >= eVF[inxKHY(kkk,h,y)] && (COHendo[inxKHYS(kkk,h,y,s)] > COHendo[inxKHYS(kk,h,y,s)])){
                            COHendo[inxKHYS(kkk,h,y,s)] = -10000; // we don't care about this saving point, since we have a better point.
                            remove = 0;
                            kkk = kkk-1;
                            slope = (eVF[inxKHY((kk+1),h,y)] - eVF[inxKHY((kk),h,y)])/(COHendo[inxKHYS((kk+1),h,y,s)] - COHendo[inxKHYS((kk),h,y,s)]);
                            extrapole = slope*(COHendo[inxKHYS(kkk,h,y,s)] - COHendo[inxKHYS(kk,h,y,s)]) + eVF[inxKHY((kk),h,y)]; // may be before kk (in that case controlled).
                            if(kkk < 0){remove = 1;}
                        }

                        if(extrapole < eVF[inxKHY(kkk,h,y)] && (COHendo[inxKHYS(kkk,h,y,s)] > COHendo[inxKHYS(kk,h,y,s)])){
                            COHendo[inxKHYS(kk,h,y,s)] = -10000;
                            kk = kk+1;
                            remove = 0;
                            if(kk == length_k-1){
                                remove = 1;
                            }else{
                            slope = (eVF[inxKHY((kk+1),h,y)] - eVF[inxKHY((kk),h,y)])/(COHendo[inxKHYS((kk+1),h,y,s)] - COHendo[inxKHYS((kk),h,y,s)]);
                            extrapole = slope*(COHendo[inxKHYS(kkk,h,y,s)] - COHendo[inxKHYS(kk,h,y,s)]) + eVF[inxKHY((kk),h,y)];
                            
                            //printf("kk too big, extrapolate = %f, VF = %f, VF-10 = %f, VF+1= %f , slope = %f, kk = %d, kkk = %d, COHendoNew = %f, COHendoNew = %f, COHendoNew = %f, COHendoNew = %f, COHendocompare = %f", extrapole, eVF[inxKHY((kk),h,y)], eVF[inxKHY((kk-10),h,y)], eVF[inxKHY((kk+1),h,y)], slope, kk, kkk, COHendo[inxKHYS(kk,h,y,s)],COHendo[inxKHYS(kk-1,h,y,s)],COHendo[inxKHYS(kk-2,h,y,s)],COHendo[inxKHYS(kk-3,h,y,s)], COHendo[inxKHYS(kkk,h,y,s)]);
                                
                            }
                        }

                        if((COHendo[inxKHYS(kkk,h,y,s)] < COHendo[inxKHYS(kk,h,y,s)])){
                            remove = 1;
                        }

                    } // end remove some COH values.


                    // THEN SHOULD HAVE A FINE GRID WITH ONLY ONE SOLUTION FOR COH VALUE.
                    klim = kk;
                    
                    } // end Kcase == 2

                    if(k >= 2){
//                        if(COHendo[inxKHYS((k-1),h,y,s)] < COHendo[inxKHYS((k-2),h,y,s)] && COHendo[inxKHYS((k-1),h,y,s)] > -1000){printf("BUG COHk < COHk-1: %f %f", COHendo[inxKHYS(k-1,h,y,s)], COHendo[inxKHYS((k-2),h,y,s)]);getchar();
//                        }
                    }
                    
                    //printf("%f %f %f %f %f %f %f %f %f %f", COHendo[inxKHYS(k+4,h,y,s)], COHendo[inxKHYS(k+3,h,y,s)], COHendo[inxKHYS(k+2,h,y,s)], COHendo[inxKHYS(k+1,h,y,s)], COHendo[inxKHYS(k,h,y,s)], COHendo[inxKHYS(k-1,h,y,s)], COHendo[inxKHYS(k-2,h,y,s)], COHendo[inxKHYS(k-3,h,y,s)], COHendo[inxKHYS(k-4,h,y,s)], COHendo[inxKHYS(k-5,h,y,s)]);getchar();
                    
                
                } // end work around the kinks.
                } // end of k> 0.

            } // end k
        } // end s
    } // end y
} // end h


} // end UEC

