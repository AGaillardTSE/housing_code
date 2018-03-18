/* **************************************************** */
//      REPLICATION SOMMER ET SULLIVAN (AER, 2017)      //
/*                  February 2018                       */
//               Alexandre Gaillard                     //
/* **************************************************** */



/** GET RESIDUALS OF THE HISTOGRAM **/
void getresid(double x, double *residout, int *ixgrid, int h)
{
    double xgrid;
    
    // PUT THESE VALUE ON A GRID //
    if(h > 0){
        xgrid = invlinspace((-(1-theta)*price*H[h]),kmax,length_k,x);
    }else{
        xgrid = invlinspace(0,kmax,length_k,x);
    }
    
    *ixgrid=min((length_k-1),(int)(floor(xgrid)));
    *ixgrid=max(0,*ixgrid);
    
    if (*ixgrid>=(length_k-1))
    {
        //printf("getresid: (ixgrid>(dimweight-1))");getchar();
        *ixgrid=min((length_k-1),(int)(floor(xgrid)));
        *ixgrid=max(0,*ixgrid);
    }
    
    if (*ixgrid<=0)
    {
        //printf("getresid: (ixgrid<0) %d", *ixgrid);getchar();
        *ixgrid=min((length_k-1),(int)(floor(xgrid)));
        *ixgrid=max(0,*ixgrid);
    }
    
    *residout=(x-K[inxKH(*ixgrid,h)])/(K[inxKH((*ixgrid+1),h)]-K[inxKH(*ixgrid,h)]);
}




/** SIMULATION **/
void SIMULATION(double *save, int *Hopti, int *Sopti, double *dist, double *aggR, double *aggH, double *currentstock){
   
   
int k, y, h, i, *isave, ynext;
double *distold, *saveres, critdist, *distasset, housesold;

// Allocate the old and new dist //
distasset = (double *) calloc((dimKH), sizeof(double));
distold = (double *) calloc((dimKHY), sizeof(double));
saveres = (double *) calloc((dimKHY), sizeof(double));
isave = (int *) calloc((dimKHY), sizeof(int));

// Compute RESIDUAL distribution (FOR HISTOGRAM METHOD) //
for(k = 0; k < length_k; k++){
    for(h = 0; h < length_h; h++) {
        for(y = 0; y < length_y; y++) {
            if(save[inxKHY(k,h,y)]>kmax){save[inxKHY(k,h,y)]=kmax-0.00000001;} // printf("BUG up %f %f",save[inxKHY(k,h,y)],kmax);getchar();
            if(save[inxKHY(k,h,y)]<K[inxKH(0,Hopti[inxKHY(k,h,y)])]){printf("BUG %f %f",save[inxKHY(k,h,y)],K[inxKH(0,h)]);getchar();save[inxKHY(k,h,y)]=K[inxKH(0,Hopti[inxKHY(k,h,y)])]+0.00000001;}
            getresid(save[inxKHY(k,h,y)],&saveres[inxKHY(k,h,y)],&isave[inxKHY(k,h,y)],Hopti[inxKHY(k,h,y)]);
        }
    }
}




// STATIONARY DISTRIBUTION //
critdist=1.0;


while(critdist>epsdist)
{

// First bascule last iteration dist to the old one //
bascule(dist, distold, dimKHY);

for(i=0;i<dimKHY;i++)
{
    dist[i]=0.0;
}

housesold = 0.0;
for(k = 0; k < length_k; k++){
    for(h = 0; h < length_h; h++) {
        for(y = 0; y < length_y; y++) {
            for(ynext = 0; ynext < length_y; ynext++) {
            
                if(y == 6 && ynext == 0){ // individual die, so becomes renter with 0 net worth
                    
                    // first histogram, (1 - saveres) jump to isave-1 value
                    dist[0] += PI_tot[y][ynext]*(distold[inxKHY(k,h,y)]);
                    housesold += H[h]*PI_tot[y][ynext]*(distold[inxKHY(k,h,y)]);
                    
                }else{
                
                    // first histogram, (1 - saveres) jump to isave-1 value
                    dist[inxKHY(isave[inxKHY(k,h,y)],Hopti[inxKHY(k,h,y)],ynext)] += PI_tot[y][ynext]*(1.0-saveres[inxKHY(k,h,y)])*(distold[inxKHY(k,h,y)]);
                    
                    // second histogram saveres fraction of individual switch to isave value
                    dist[inxKHY((min((isave[inxKHY(k,h,y)]+1),(length_k-1))),Hopti[inxKHY(k,h,y)],ynext)] += PI_tot[y][ynext]*saveres[inxKHY(k,h,y)]*(distold[inxKHY(k,h,y)]);
                    
                } // end if
            } // end y
        } // end ygrid
    } // end hgrid
} // end igrid


// convergence criterion //
critdist=0.0;
for(k = 0; k < length_k; k++){
    for(h = 0; h < length_h; h++) {
        for(y = 0; y < length_y; y++) {
            critdist=(max(critdist,fabs(dist[inxKHY(k,h,y)]-(distold[inxKHY(k,h,y)]))));
        } // end ygrid
    } // end hgrid
} // end igrid

//printf("critdist = %f", critdist); getchar();

} // end while




/**************************************/
/** START TO COMPUTE SOME STATISTICS **/
/**************************************/

double aggK, fracL, fracH, fracR;
    
aggK = 0.0; fracL = 0.0; fracH = 0.0; *aggH =0.0; *aggR = 0.0; *currentstock = 0.0; fracR = 0.0;

for(k = 0; k < length_k; k++){
    for(h = 0; h < length_h; h++) {
        for(y = 0; y < length_y; y++) {
            
            distasset[inxKH(k,h)] += dist[inxKHY(k,h,y)];
            
            aggK += (K[inxKH(k,h)] + H[h])*dist[inxKHY(k,h,y)];
            
            *aggR += (S[Sopti[inxKHY(k,h,y)]] - H[Hopti[inxKHY(k,h,y)]])*dist[inxKHY(k,h,y)];
            
            *aggH += (H[Hopti[inxKHY(k,h,y)]])*dist[inxKHY(k,h,y)];
            
            *currentstock += (H[h])*dist[inxKHY(k,h,y)];
            
            // HOUSING TENURE //
            if(H[Hopti[inxKHY(k,h,y)]] == 0){
                fracR += dist[inxKHY(k,h,y)];
            } else {
                if(S[Sopti[inxKHY(k,h,y)]] < H[Hopti[inxKHY(k,h,y)]]){
                    fracL += dist[inxKHY(k,h,y)];
                }
                fracH += dist[inxKHY(k,h,y)];
            }
        } // end ygrid
    } // end hgrid
} // end igrid

*aggH = *aggH - housesold;

//printf("Homeowner = %f, Renter = %f, Landlord = %f, aggK = %f \n", fracH, fracR, fracL, aggK);
    

    
///** PRINT DISTRIBUTION **/
//
//// initialize files
//FILE *distfileout, *distfileasset;
//
//distfileout=fopen(distfile, "w");
//
//setbuf ( distfileout, NULL );
//
//fprintf(distfileout,"%5s\t%20s\t%5s\t%5s\t%5s\t%20s\n","ASSET","ASSETVAL", "PROD", "HOUSE", "SERVICE", "dist");
//for(k = 0; k < length_k; k++){
//    for(h = 0; h < length_h; h++) {
//        for(y = 0; y < length_y; y++) {
//            fprintf(distfileout,"%5d\t%20.15f\t%5d\t%5d\t%5d\t%20.15f\t%20.15f\t%5d\n",k,  save[inxKHY(k,h,y)], y, Hopti[inxKHY(k, h, y)], Sopti[inxKHY(k, h, y)], dist[inxKHY(k, h, y)], saveres[inxKHY(k, h, y)], isave[inxKHY(k, h, y)]);
//        }
//    }
//}
//fclose(distfileout);
//
//
//distfileasset=fopen(distassetfile, "w");
//
//setbuf ( distfileasset, NULL );
//
//fprintf(distfileasset,"%5s\t%5s\t%20s\t%20s\n","ASSET","HOUSE", "ASSETVAL", "dist");
//for(k = 0; k < length_k; k++){
//    for(h = 0; h < length_h; h++) {
//        fprintf(distfileasset,"%5d\t%5d\t%20.15f\t%20.15f\n",k, h, K[inxKH(k,h)], distasset[inxKH(k,h)]);
//    }
//}
//fclose(distfileasset);
//
    
free(distold);
free(saveres);
free(isave);

}

