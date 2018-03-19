/////////////////////
// Useful function //
/////////////////////

/******************/
//  MATH FUNCTION //
/******************/

#define max(a,b) (((a)>(b))?(a):(b))
#define min(a,b) (((a)<(b))?(a):(b))
#define interpol(x,y,z) (y+(x-floor(x))*(z-y))
#define inter1d(x1,y1,y2) ((1.0-(x1))*(y1)+(x1)*(y2))
#define inter2d(x1,x2,y11,y21,y12,y22) ((1.0-(x2))*((1.0-(x1))*(y11)+(x1)*(y21))+(x2)*((1.0-(x1))*(y12)+(x1)*(y22)))

// TO DEFINE GRID SPACING //
#define linspace(x0,xmax,n,i) ((i)*(((xmax)-(x0))/(n))+(x0)) // transform a grid into a value
#define invlinspace(x0,xmax,n,x) ((((x)-(x0))/((xmax)-(x0)))*(n)) // transform a value into a grid
#define expspace(i,xmin,xmax,echelle,n) ( echelle*(exp((log((xmax/echelle)-((xmin/echelle)-1.0))/(n-1))*(i))+((xmin/echelle)-1.0)) ) // transform a grid into a value
#define invexpspace(x,xmin,xmax,echelle,n) ( log((x)/echelle - ((xmin/echelle)-1.0))/(log((xmax/echelle)-((xmin/echelle)-1.0))/(n-1)) ) // transform a value into a grid



double timer_alex(timeval t1, timeval t2){
    double elapsedTime;
    elapsedTime = (t2.tv_sec - t1.tv_sec) * 1000.0;
    elapsedTime += (t2.tv_usec - t1.tv_usec) / 1000.0;
    
    return (elapsedTime/1000);
}


// FIND THE WEIGHT FOR INTERPOLATION //
double weightinter(double x, double xgrid, int *ixgrid, double vector[], int dimweight)
{
    double dxgrid;
    
    // PUT THESE VALUE ON A GRID //
    *ixgrid=min((dimweight-1),(int)(floor(xgrid)));
    *ixgrid=max(0,*ixgrid);
    
    if (*ixgrid>=(dimweight-1))
    {
        xgrid=(dimweight-1)-0.000000001;
        *ixgrid=min((dimweight-1),(int)(floor(xgrid)));
        *ixgrid=max(0,*ixgrid);
    }
    
    if (*ixgrid<=0)
    {
        xgrid=0.000000001;
        *ixgrid=min((dimweight-1),(int)(floor(xgrid)));
        *ixgrid=max(0,*ixgrid);
    }
    
    dxgrid=(x-vector[*ixgrid])/(vector[(*ixgrid+1)]-vector[*ixgrid]);
    
    return dxgrid;
    
}




// FOR QSORT
int comparefun2(const void* a, const void* b) 
{ 
	double* da = (double*)a; 
	double* db = (double*)b; 
	int diff1 = (da[0] > db[0]) - (da[0] < db[0]); 
	if (diff1 != 0) return diff1; 
	return (da[1] > db[1]) - (da[1] < db[1]); 
}




// SHIFTER //

namespace {
	inline void shft2(double &a, double &b, const double c)
	{
		a=b;
		b=c;
	}
    
	inline void shft3(double &a, double &b, double &c, const double d)
	{
		a=b;
		b=c;
		c=d;
	}
}

using namespace std;


void bascule(double *VectorIN, double *VectorOUT, int dim)
{
    int i;
    for(i=0; i<dim; i++){VectorOUT[i]=VectorIN[i];}
}



//// GET RESIDUALS OF THE HISTOGRAM //
//void getresid(double amylevel,double *residout,int *ixgrid, double vector[], int dimweight, double Gridmax)
//{
//    assert(residout);
//    assert(ixgrid);
//    
//    *ixgrid=min((dimweight-1),(int)(floor(phiinv((amylevel)))));
//    *ixgrid=max(0,*ixgrid);
//    
//    if (*ixgrid>(dimweight-1))
//    {
//        printf("getresid: (ixgrid>(maxgrid-1))");
//        getchar();
//    }
//    
//    if (*ixgrid<0)
//    {
//        printf("getresid: (ixgrid<0)");
//        getchar();
//    }
//    
//    
//    if ((amylevel>(Gridmax))){printf("bascule: (xval>(Gridmax)) %20.15f\n",amylevel);getchar();}
//    if ((amylevel<(0.0))){printf("bascule: (xval<(0.0)) %20.15f\n",amylevel);getchar();}
//    
//    
//    *residout=(amylevel-vector[*ixgrid])/(vector[(*ixgrid+1)]-vector[*ixgrid]);
//}




// INVARIANT DISTRIBUTION //
template<size_t dim, size_t dim2>
void inv_distri(double (&invdist)[dim], double (&prob)[dim][dim2])
{
    double tempdist[dim], critdist, sumdist;
    int i, j;
    
    invdist[1] = 1.0;
    
    critdist = 1.0;
    while(critdist > 0.00000001) {
        
        for(i=0;i<dim;i++){
            tempdist[i] = invdist[i];
        }
        
        // compute the invdist //
        for(i = 0; i<dim; i++){
            invdist[i] = 0.0;
        }
        
        for(i = 0; i<dim; i++){
            for(j = 0; j<dim2; j++) {
                invdist[i] += tempdist[j]*prob[j][i];
            }
        }
        
        critdist = 0.0;
        for(i=0;i<dim;i++) {
            critdist = max(abs(invdist[i] - tempdist[i]), critdist);
        }
        
       //printf("%f %f %f %f %f %f %f", invdist[0], invdist[1], invdist[2], tempdist[0], tempdist[1], tempdist[2], critdist); getchar();
        
    }
    
    // renormalize invdist //
    sumdist = 0.0;
    for(i = 0; i<dim; i++) {
        sumdist += invdist[i];
    }
    for(i=0;i<dim;i++) {
        invdist[i] = invdist[i]/sumdist;
    }
    
}
