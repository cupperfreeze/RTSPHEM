/* Programmer: Stephan Gärttner                          */

/* calculate (P*(R\(R'\(P'*x)))) for permutation matrix P, upper triangular R and 2 rhs (x vector double length)*/

#include "mex.h"
#include "matrix.h"

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
/* read input */
    mwSize m, n, nrow;
    double *Mpr1, *Vpr, *Cpr;
    mwIndex *Mir1, *Mjc1,*Mir2;

    Mpr1 = mxGetPr(prhs[0]);	
    Mir1 = mxGetIr(prhs[0]);
    Mjc1 = mxGetJc(prhs[0]);   	
    Mir2 = mxGetIr(prhs[1]);
    Vpr = mxGetPr(prhs[2]);	

    m = mxGetM(prhs[0]);
    n = mxGetN(prhs[0]);
  
    
/* Create output */
    plhs[0] = mxCreateDoubleMatrix(2*n, 1,mxREAL);
    Cpr = mxGetPr(plhs[0]);

  	
  	
/* Calculate result */

	int j;
	int i;
	
	double temp[2*n]; 
	int nnrow[n];

	/* sort according to permutation */
	for(j=0; j<n; j++ ) {
		temp[j] = Vpr[Mir2[j]];
		temp[j+n] = Vpr[Mir2[j]+n];
        }


	/* first substitution (forward) */
	for(j=0; j<n; j++ ) {nnrow[j]= Mjc1[j+1] - Mjc1[j];}
		
	temp[0] = temp[0]/ Mpr1[0];
	temp[n] = temp[n]/ Mpr1[0];
    		for(j=1; j<n; j++ ) {
         		for(i=0; i<nnrow[j]-1; i++ )  {
            			temp[j] -= Mpr1[Mjc1[j]+i] * temp[Mir1[Mjc1[j]+i]];  
				temp[j+n] -= Mpr1[Mjc1[j]+i] * temp[Mir1[Mjc1[j]+i]+n];  
				}
			temp[j] = temp[j]/ Mpr1[Mjc1[j+1]-1];
			temp[j+n] = temp[j+n]/ Mpr1[Mjc1[j+1]-1];
			
	 }	
	
	
	/* second substitution (backward) */
	int offset2 = Mjc1[n]-2;
	temp[n-1] = temp[n-1]/ Mpr1[offset2+1];
	temp[2*n-1] = temp[2*n-1]/ Mpr1[offset2+1];
	for(j=0; j<n; j++ ) {nnrow[j]= Mjc1[j+1] - Mjc1[j];}

	

    		for(j=n-1; j>0; j-- ) {
		        		/* Number of row elements for this column */
         		for(i=0; i<nnrow[j]-1; i++ )  {
            			temp[Mir1[offset2-i]] -= Mpr1[offset2-i] * temp[j]; 
 				temp[Mir1[offset2-i]+n] -= Mpr1[offset2-i] * temp[j+n]; 
			}
			offset2=Mjc1[j]-2;
			temp[j-1] = temp[j-1]/ Mpr1[offset2+1];
			temp[j-1+n] = temp[j-1+n]/ Mpr1[offset2+1];
 
	       }	

	
	/* sort according to permutation */
	for(j=0; j<n; j++ ) {
		Cpr[Mir2[j]] = temp[j];
		Cpr[Mir2[j]+n] = temp[j+n];
        }
           
		

	


}