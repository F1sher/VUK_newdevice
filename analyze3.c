#include <stdio.h>
#include <stdlib.h>
#include <fcntl.h>
#include <math.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <string.h>
#include <unistd.h>
#include <time.h>
#include <errno.h>

#define	X_POINTS	256
#define SIZEOF_DATA		X_POINTS 
#define CRIT_MAX	0.2
#define STEP_SQR	.01
#define EPS .000001

#define sqr(x)	((x)*(x))

double *square_osci[2], *time_osci[2];
float ucalc[9];
int num_arr_el = 0;

int DenisCalc(const int N[], const unsigned int F[], int LCH, int RCH, char flag, float time_N, float time_F, float u[])
{
// Ucalc[0] = square
// Ucalc[1] = Ucalc[0] error
// Ucalc[2] = pick position
// Ucalc[3] = Ucalc[2] error
// Ucalc[4] = FWHM
// Ucalc[5] = k (background approximation)
// Ucalc[6] = b (background approximation)
// Ucalc[7] = Ucalc[4] error
// Ucalc[8] = clear square (with background)
	
/* IF flag == 0 BEGIN */
	if (flag == 0) {
		/* IF no data */
		if (RCH <= LCH) return -13;
		
		static double *Q;
		double time_ratio1 = 1.0, time_ratio2 = 1.0;
		double avlb, avrb, k, b;
		int i, j;
		double s = 0.0, s_err = 0.0;
		double p = 0.0, p_err = 0.0;
		int c = 0, e = 0 ;
		static double *G;
		double max, hm;
		double t;
		
		/* Allocate memory for Q array + test */
		Q = (double *)malloc(sizeof(double) * (RCH+3));
		if (Q == NULL) return -1;
	
		/* Calculate time ratio */
		if (time_N > time_F) 
		{  
 	     // printf("flag=%d time_N=%e  time_F=%e \n",(int)flag,time_N , time_F );
			if (time_F < EPS) return -2;
			time_ratio1 = time_N/time_F;
		}
		
		if (time_N < time_F) {
			if (time_N < EPS) return -3; 
			time_ratio2 = time_F/time_N;
		}
	
		for (i = 0; i <= RCH+2; i++) {
			double a;
			a = ((double) N[i])/time_ratio1 - ((double) F[i])/time_ratio2;
			if (a > EPS) 
				Q[i] = a;
			else Q[i] = 0.0;
		}
	
		/* Calculate b, k of the background line (y=k*x+b) */
		avlb = (Q[LCH-2]+Q[LCH-1]+Q[LCH])/3.0;
		avrb = (Q[RCH]+Q[RCH+1]+Q[RCH+2])/3.0;
		k = (avrb-avlb)/(RCH-LCH+2);
		b = Q[LCH] - k*(LCH);
	
		/* Calculate S, S error, p */
		for(i = LCH; i <= RCH; i++) {
			if ((t = Q[i] - (k*i + b)) < EPS ) t = 0.0;
			s = s + t;
			if (t > EPS)
				s_err = s_err + (Q[i] + k*i+b);
			p = p + t*i;
		}
	
		/* IF s<=0 or p<=0 only */
		if (s < EPS || p < EPS) 
		{
			printf("s=%e  p=%e  EPS=%f \n",s ,p ,(float)EPS  );
			return -14;
		}
	
		u[0] = (float) s; 
		u[1] = (float) sqrt(s_err);
		u[2] = (float) (p = p/s);
	
		/* Calculate p error */
		for(i = LCH; i <= RCH; i++) {
			  if((t = Q[i] - (k*i+b)) < EPS) t = 0.0;
		      p_err = p_err + t*(p - i)*(p - i);
 			//p_err = p_err + (float) ((Q[i] - (k*i+b))*sqr(p - i));
 		}
		
		/* IF p_err <= 0 */	
		if (p_err < EPS) return -4;	
		
		u[3] = (float) ( p_err = sqrt(p_err/((RCH-LCH)*s)) );
	
		/* Allocate memory for G array + test */		
		G = (double *)malloc(sizeof(double) * (RCH-LCH+1));
		if (G == NULL) { free(Q); return -5; }
	
		/* Calculate G[K] = Q[K] - B[K] */
		for(i = LCH; i <= RCH; i++) {
			if ((t = Q[i] - (k*i+b)) < EPS) t = 0.0;
			G[i-LCH] = t;
		}
		
		/* Find max element in G[] array */
		for(i = 0, j = 0, max = G[0]; i <= RCH-LCH; i++)
			if (G[i] >= max) {
				max = G[i];
				j = i;
			}
		
		/* Half maximum */
		hm = G[j]/2;
	
		/* Calculate c and e */
		for(i = 0; i <= RCH-LCH; i++) {
			if (G[i] < hm && G[i+1] > hm && i < j)
				c = i;
			if (G[i] > hm && G[i+1] < hm && i > j)
				e = i;
		}
		
		/* IF if G[c+1]-G[c] == 0 or G[e]-G[e+1] == 0 */
		if ((G[c+1]-G[c]) < EPS || (G[e]-G[e+1]) < EPS) return -6; 
		
		u[4] = (float) ((e-c-1)+(G[c+1]-hm)/(G[c+1]-G[c])+(G[e]-hm)/(G[e]-G[e+1]));
	
		free(Q);
		free(G);
		
	return 0;
/* IF flag == 0 END */
	}
/* IF flag == 1 BEGIN */
	if (flag !=0 ) {
		double avlb, avrb, k, b;
		int i, j;
		double s = 0.0, s_err = 0.0;
		double s_en = 0.0, s_err_en = 0.0;
		double p = 0.0, p_err = 0.0;
		int c = 0, e = 0 ;
		static double *G;
		double max, hm;
		double t;
		double sqr_clear = 0.0;
	
		/* IF no data */
		if (RCH <= LCH) return -7;
	
		for(i=LCH; i<=RCH; i++)
			sqr_clear = sqr_clear + (double)N[i];
		
		u[8] = (float)sqr_clear;
		printf("sqr clear = %.2f u[7] = %.2f\n", sqr_clear, u[7]);
	
		/* Calculate b, k of the background line (y=k*x+b) */
		avlb = (N[LCH-2]+N[LCH-1]+N[LCH])/3.0;
		avrb = (N[RCH]+N[RCH+1]+N[RCH+2])/3.0;
		k = (avrb-avlb)/(RCH-LCH+2);
		b = N[LCH] - k*(LCH);
	
		u[5] = k;
		u[6] = b;
	
		//k = b = 0.0;
		/* Calculate S, S error, p */
		for(i = LCH; i <= RCH; i++) {
			if ((t = N[i] - (k*i + b)) < EPS ) t = 0.0;
			s = s + t;
			if (t > EPS)
				s_err = s_err + (N[i] + k*i+b);
			p = p + t*i;
		}
	
		/* IF s<=0 or p<=0 only */
		if (s < EPS || p < EPS) return -8;
	
		u[0] = (float) s; 
		u[1] = (float) sqrt(s_err);
		u[2] = (float) (p = p/s);
	
		/* Calculate p error */
		for(i = LCH; i <= RCH; i++) {
			if((t = N[i] - (k*i+b)) < EPS) t = 0.0;
 		    p_err = p_err + t*(p - i)*(p - i);
//			p_err = p_err + (float) ((N[i] - (k*i+b))*sqr(p - i));
		}
			
		/* IF p_err <= 0 */	
		if (p_err < EPS) return -9;
			
		u[3] = (float) ( p_err = sqrt(p_err/((RCH-LCH)*s)) );
		
		/* Allocate memory for G array + test */		
		G = (double *)malloc(sizeof(double) * (RCH-LCH+1));
		if (G == NULL) return -10;
	
		/* Calculate G[K] = N[K] - B[K] */
		for(i = LCH; i <= RCH; i++) {
			if ((t = N[i] - (k*i+b)) < EPS) t = 0.0;
			G[i-LCH] = t;
		}
		
		/* Find max element in G[] array */
		for(i = 0, j = 0, max = G[0]; i <= RCH-LCH; i++)
			if (G[i] >= max) {
				max = G[i];
				j = i;
			}
		
		/* Half maximum */
		hm = G[j]/2;
	
		/* Calculate c and e */
		for(i = 0; i <= RCH-LCH; i++) {
			if (G[i] < hm && G[i+1] > hm && i < j)
				c = i;
			if (G[i] > hm && G[i+1] < hm && i > j)
				e = i;
		}
		
		/* IF if G[c+1]-G[c] == 0 or G[e]-G[e+1] == 0 */
		if ((G[c+1]-G[c]) < EPS || (G[e]-G[e+1]) < EPS) return -11; 
		
		u[4] = (float) ((e-c-1)+(G[c+1]-hm)/(G[c+1]-G[c])+(G[e]-hm)/(G[e]-G[e+1]));
		
		free(G);
		
		u[7] = (float) (( (G[c+1]-hm)/(G[c+1]-G[c])+(G[e]-hm)/(G[e]-G[e+1]) )*sqrt( 1.0/fabs(G[c+1]-hm) + 1.0/fabs(G[c+1]-G[c]) + 1.0/fabs(G[e]-hm) + 1.0/fabs(G[e]-G[e+1]) ));
		
	return 0;
/* IF flag == 1 END */	
	}	

return -12;	
}

char isNum(char *s)
{
	while(*s++) {
		if(*s == '-')
			return 0;
	}
			
	return 1;
}

double max_bubble(double *x, int nums)
{
	int i;
	double z;
	
	for(i=1, z = fabs(x[0]); i<=nums-1; i++)
		if(fabs(x[i]) > z)
			z = fabs(x[i]);
	
	return z;
}

int min_bubble_num(int *x, int nums)
{
	int i, z;
	int max_num;
	
	for(i=1, z = x[0]; i<=nums-1; i++)
		if(x[i] < z) {
			z = x[i];
			max_num = i;
		}
	
	return max_num;
}

int det3(int a[3][3])
{
	return (a[0][0]*a[1][1]*a[2][2])-(a[0][0]*a[1][2]*a[2][1])
		+(a[0][1]*a[1][2]*a[2][0])-(a[0][1]*a[1][0]*a[2][2])
		+(a[0][2]*a[1][0]*a[2][1])-(a[0][2]*a[1][1]*a[2][0]);
}

double find_start_pick_lsm(int *x)
{
	const int start_n = 47, end_n = 51; // number of points to aproximation
	int i, j;
	int A[3][3], B[3]; // A - matix for Least Square Method
	double X[3];
	double z;
	double start_pos;
	
	for(i=0; i<3; i++) {
		B[i] = 0;
		for(j=0; j<3; j++)
			A[i][j] = 0;
	}
	
	for(i=start_n; i<=end_n; i++) {
		A[0][0] += i*i;
		A[0][1] += i;
		A[1][0] += i*i*i;
		A[2][0] += i*i*i*i;
		
		B[0] += x[i];
		B[1] += i*x[i];
		B[2] += i*i*x[i];
	}

	A[1][1] = A[2][2] = A[0][0];
	A[0][2] = (end_n-start_n+1);
	A[2][1] = A[1][0];
	A[1][2] = A[0][1];
/*
	for(i=0; i<3; i++) {
		printf("B[%d]=%d\n", i, B[i]);
		for(j=0; j<3; j++)
			printf("A[%d][%d]=%d ", i, j, A[i][j]);
		printf("\n");
	}*/

	int detx[3][3] = {{B[0],A[0][1],A[0][2]},{B[1],A[1][1],A[1][2]},
						{B[2],A[2][1],A[2][2]}};
	int dety[3][3] = {{A[0][0],B[0],A[0][2]},{A[1][0],B[1],A[1][2]},
						{A[2][0],B[2],A[2][2]}};
	int detz[3][3] = {{A[0][0],A[0][1],B[0]},{A[1][0],A[1][1],B[1]},
						{A[2][0],A[2][1],B[2]}};
	
	if(det3(A) != 0) {
		X[0] = (double)det3(detx)/det3(A);
		X[1] = (double)det3(dety)/det3(A);
		X[2] = (double)det3(detz)/det3(A);
	}
	/*
	for(i=0; i<3; i++)
		printf("X[%d] = %.2f ", i, X[i]);
	printf("\n");
	*/
	for(i=1, z = x[0]; i<=SIZEOF_DATA-2; i++)
		if(x[i] > z && x[i-1] >= 0.2*x[i] && x[i+1] >= 0.2*x[i]) {z = x[i]; j = i;}
	z = ((x[j]+x[j+1]+x[j+2])/3.0-(x[0]+x[1]+x[2])/3.0)*CRIT_MAX + (x[0]+x[1]+x[2])/3.0;
	
	double step = 0.001;
	double xi, xi_next;
	int numOFiter = (end_n-start_n+1)*1000;
	for(i=0; i<=numOFiter; i+=1) {
		xi = i*step+start_n;
		if( ((X[0]*xi*xi + X[1]*xi + X[2] - z) >= EPS)) return xi;
	}
	
	return 0;
}

double second_derivative_zero(int *x)
{
	int i;
	double *derivative = (double *) calloc(X_POINTS, sizeof(double));
	double start_pos=0;
	FILE *in;
	
	for(i=1; i<=X_POINTS-2; i++) {
		derivative[i] = (double)(x[i+1]-2.0*x[i]+x[i-1]);
		printf("%d %e %d %d %d\n", i, derivative[i], x[i+1], x[i], x[i-1]);
	}
		
	for(i=1; i<=X_POINTS-2; i++)
		if(derivative[i] > EPS && derivative[i+1] < EPS)
			start_pos = -1.0*(derivative[i]-i*(derivative[i+1]-derivative[i]))/(derivative[i+1]-derivative[i]);
	
	free(derivative);
	
	return start_pos;
}

double trap_area(int *a, int max_num, int min_num)
{
	int i;
	double area = 0.0;
	double baseline = (a[0] + a[1] + a[2] + a[3] + a[4]) / 5.0;
	
	for(i=max_num; i<=min_num-1; i++)
		area += 0.5*(a[i]+a[i+1]-2.0*baseline);
		
	return fabs(area);
}

double simpson_area(long int *a, int max_num, int min_num)
{
	int i;
	double area = 0.0;
	
	area = a[max_num] + a[min_num];
	for(i=max_num/2+1; i<=min_num/2; i++)
		area += 4.0*a[2*i-1];
	for(i=max_num/2+1; i<=min_num/2-1; i++)
		area += 2.0*a[2*i];
	
	return (area/3.0);
}

double monte_area(long int *a, int max_num, int min_num)
{
	int i;
	const int itt = 400000;
	double t, area = 0.0;
	
	for(i=0; i<=itt-1; i++) {
		t = (double)rand()/RAND_MAX;
		area += (double)a[(int)((min_num-max_num)*t+max_num)];
	}
	
	area = area/itt;
	
	return ((min_num-max_num)*area);
}

double monte_area2(long int *a, int max_num, int min_num)
{	
	int i;
	const int itt = 40000;
	double t, area = 0.0;

	srandom(a[0]);
	
	for(i=0; i<=itt-1; i++) {
		t = (double)random()/RAND_MAX;
		area += (double)a[(int)((min_num-max_num)*t+max_num)];
	}
	
	area = area/itt;
	
	return ((min_num-max_num)*area);
}

double monte_area3(long int *a, int max_num, int min_num)
{	
	int i;
	const int itt = 400000;
	double t, area = 0.0;

	srand48(a[0]);
	
	for(i=0; i<=itt-1; i++) {
		t = drand48();
		area += (double)a[(int)((min_num-max_num)*t+max_num)];
	}
	
	area = area/itt;
	
	return ((min_num-max_num)*area);
}

double newton_cotes_3(long int *a, int max_num, int min_num)
{
	int i;
	double area = 0.0;
	
	area = a[max_num] + a[min_num];
	
	for(i=max_num/3; i<=min_num/3-1; i++)
		area += 3.0*(a[3*i+1]+a[3*i+2]);
	
	for(i=max_num/3; i<=min_num/3-1; i++)
		area += 2.0*(a[3*i+3]);
		
	return (3.0*area/8.0);
}

double newton_cotes_4(long int *a, int max_num, int min_num)
{
	int i;
	double area = 0.0;
	
	area = 7.0*(a[max_num] + a[min_num]);
	
	for(i=max_num/2; i<=min_num/2-1; i++)
		area += 32.0*(a[2*i+1]);
	
	for(i=max_num/4; i<=min_num/4-1; i++)
		area += 12.0*(a[4*i+2]);
		
	for(i=max_num/4; i<=min_num/4-1; i++)
		area += 14.0*(a[4*i+4]);
		
	return (2.0*area/45.0);
}

int find_pick_start_stop(int *a, int *max_num, int *min_num)
{
	int i, j;
    int baseline = (a[0] + a[1] + a[2] + a[3]) / 4;
     
    for(i=2, j=1; i<=SIZEOF_DATA-1; i++) {
        if(a[i] <= 0.95*a[j]) {
		//	printf("i = %d, a[i] = %d, a[j] = %d\n", i, a[i], a[j]);
			if( (a[i+1] <= a[j]) && (a[i+2] <= a[j]) && (a[i] <= baseline) )
                j = i;
                break;
        }
    }
    *max_num = j-7;
    
    for(i=j; i <= SIZEOF_DATA-4; i++) {
		if((a[i] >= a[j-7]) &&	(a[i+1]>=a[j-7]) && (a[i]>=baseline)) {
			*min_num = i; 
			return 0;
		}
	}
    
    //*min_num = 150;
    
    return -1;
}

double stats(double *s_trap, double *avg)
{
	int i;
	double sqr_sums = 0.0, sums = 0.0;
	double disp;
		
	for(i=0; i<=num_arr_el-1; i++) {
		sqr_sums += sqr(s_trap[i]);
		sums += s_trap[i];
	}
	*avg = sums;
	sums = sqr(sums);
		
	disp = (sqr_sums - sums/num_arr_el)/(num_arr_el-1);
		
	return sqrt(disp);
}

int energy_histo(int **a, const char *filename_to_save, int ADC_NUM)
{
	printf("\nENERGY HISTO IN WORK\n");
	
	int i, j, m;
	double **delta = (double **)malloc(num_arr_el*sizeof(double *));
	if(delta == NULL) {
		perror("Error in delta malloc");
		return -1;
	}
	
	int max, min;
	double *s_trap = (double *) malloc(num_arr_el*sizeof(double));
	if(s_trap == NULL) {
		perror("Error in s_trap malloc");
		return -1;
	}
	
	int *shiftted_array = (int *)calloc(X_POINTS, sizeof(int));
	if(shiftted_array == NULL) {
		perror("Error in shiftted_array calloc");
		return -1;
	}
	
	double *k = (double *)calloc(X_POINTS, sizeof(double));
	if(k == NULL) {
		perror("Error in k calloc");
		return -1;
	}
	
	double *b = (double *)calloc(X_POINTS, sizeof(double));
	if(b == NULL) {
		perror("Error in b calloc");
		return -1;
	}
	
	FILE *in_to_spk, *in_to_shift;
	char shiftname[40];
	printf("0.2\n");
    double max_trap;
    
    int *s_trap_int = (int *)calloc(num_arr_el, sizeof(int));
	if(s_trap_int == NULL) {
		perror("Error in s_trap_int calloc");
		return -1;
	}
    
    int *hiztogram_trap = (int *)calloc(4096, sizeof(int));
    if(hiztogram_trap == NULL) {
		perror("Error in hiztogram_trap calloc");
		return -1;
	}
	
	for(i=0; i<=num_arr_el-1; i++) {
		delta[i] = (double *)malloc((int)(sizeof(double)*1.0/STEP_SQR));
		if(delta[i] == NULL) {
			perror("Error in delta[i] malloc");
			return -1;
		}
	}
	
	printf("0 num_arr_el = %d\n", num_arr_el);
	for(i=0; i<=num_arr_el-1; i++) {
		find_pick_start_stop(a[i], &max, &min);
		printf("max:%d, min:%d\n", max, min);

		s_trap[i] = trap_area(a[i], max, min);
		//s_trap[i] = trap_area(a[i], 1, X_POINTS-2);
	
		delta[i][0] = fabs(s_trap[0]-s_trap[i]);
	//	printf("delta[%d][0] = %.2f\n", i, delta[i][0]);
	}

	double avg = 0.0;
	
	printf("%.2f +- %.2f\n", avg/(num_arr_el), stats(s_trap, &avg));
	
	/*
	sprintf(shiftname, "./shift/a0");
	in_to_shift = fopen(shiftname, "w+");
	for(i=0; i<=X_POINTS-2; i++)
		fprintf(in_to_shift, "%d %d\n", i, a[0][i]);
	fclose(in_to_shift);
	
	
	for(m=1; m<=num_arr_el-1; m++) {
		find_pick_start_stop(a[m], &max, &min);

		
		for(i=0; i<=X_POINTS-2; i++) {
			k[i] = (double)(a[m][i+1]-a[m][i]);
			b[i] = a[m][i]-k[i]*(double)i;
			shiftted_array[i] = k[i]*((double)i+STEP_SQR)+b[i];
		}
		
		j=0;
		sprintf(shiftname, "./shift/shift-m%d-j%d", m, j);
		in_to_shift = fopen(shiftname, "w+");
		for(i=0; i<=X_POINTS-2; i++)
			fprintf(in_to_shift, "%d %d\n", i, a[m][i]);
		fclose(in_to_shift);
		
		j=1;
		sprintf(shiftname, "./shift/shift-m%d-j%d", m, j);
		in_to_shift = fopen(shiftname, "w+");
		for(i=0; i<=X_POINTS-2; i++)
			fprintf(in_to_shift, "%d %d\n", i, shiftted_array[i]);
		fclose(in_to_shift);
		
		find_pick_start_stop(shiftted_array, &max, &min);
		
		//printf("max: (%d, %d), min: (%d, %d) for m = %d, s = %.1f\n", max, a[m][max], min, a[m][min], m, trap_area(a[m], max, min));
		
		//if(...) -> SHIFT_TO_RIGHT else -> SHIFT_TO_LEFT
		delta[m][1] = fabs(s_trap[0]-trap_area(shiftted_array, max, min));
		printf("delta[%d][0] = %.2f; %.2f (max:%d, min:%d)\n", m, delta[m][0], delta[m][1], max, min);
		if( delta[m][1] < delta[m][0] ) {
			printf("IF#\n");
			for(j=2; j<=(int)(1.0/STEP_SQR)-1; j++) {
				for(i=0; i<=X_POINTS-2; i++)
					shiftted_array[i] = k[i]*((double)i+(double)j*STEP_SQR)+b[i];
				
				
				sprintf(shiftname, "./shift/shift-m%d-j%d", m, j);
				in_to_shift = fopen(shiftname, "w+");
				for(i=0; i<=X_POINTS-2; i++)
					fprintf(in_to_shift, "%d %d\n", i, shiftted_array[i]);
				fclose(in_to_shift);
				
				find_pick_start_stop(shiftted_array, &max, &min);
			
				delta[m][j] = fabs(s_trap[0]-trap_area(shiftted_array, max, min));
				//delta[m][j] = fabs(s_trap[0]-trap_area(shiftted_array, 1, X_POINTS-2));
			
				if(delta[m][j] <= delta[m][j-1]) {
					continue; 
				}
				else {
					//printf("END delta[%d] = %.4f delta[%d] = %.4f, delta[0] = %.4f\n", j, delta[m][j], j-1, delta[m][j-1], delta[m][0]);
				
					for(i=0; i<=X_POINTS-2; i++)
						shiftted_array[i] = k[i]*((double)i+(double)(j-1)*STEP_SQR)+b[i];
					
					find_pick_start_stop(shiftted_array, &max, &min);

					s_trap[m] = trap_area(shiftted_array, max, min);
					//s_trap[m] = trap_area(shiftted_array, 1, X_POINTS-2);
					
					//printf("PLUS SQR0 = %.4f; SQR1 = %.4f m = %d j = %d\n", s_trap[0], s_trap[m], m, j-1);
					//printf("--------------------------------------------\n");
				
					break;
				}
			}
		}
		else {
			printf("ELSE#\n");
			
			for(i=0; i<=X_POINTS-2; i++)
				shiftted_array[i] = k[i]*((double)i-(double)STEP_SQR)+b[i];
			
			//if( (delta[m][1] = fabs(s_trap[0]-trap_area(shiftted_array, max, min))) < delta[m][0] ) {
			if( (delta[m][1] = fabs(s_trap[0]-trap_area(shiftted_array, 1, X_POINTS-2))) < delta[m][0] ) {
			
			for(j=2; j<=(int)(1.0/STEP_SQR)-1; j++) {
				for(i=0; i<=X_POINTS-2; i++)
					shiftted_array[i] = k[i]*((double)i-(double)j*STEP_SQR)+b[i];
				
			
				find_pick_start_stop(shiftted_array, &max, &min);
			
				delta[m][j] = fabs(s_trap[0]-trap_area(shiftted_array, max, min));
				//delta[m][j] = fabs(s_trap[0]-trap_area(shiftted_array, 1, X_POINTS-2));
			
				if(delta[m][j] <= delta[m][j-1]) {printf("delta[%d][%d] = %.2f delta[%d][%d] = %.2f\n", m, j, delta[m][j], m, j-1, delta[m][j-1]); continue; }
				else {
					printf("m = %d END delta[%d] = %.4f delta[%d] = %.4f, delta[0] = %.4f\n", m, j, delta[m][j], j-1, delta[m][j-1], delta[m][0]);
				
					for(i=0; i<=X_POINTS-2; i++) 
						shiftted_array[i] = k[i]*((double)i-(double)(j-1)*STEP_SQR)+b[i];
					
					s_trap[m] = trap_area(shiftted_array, max, min);
					//s_trap[m] = trap_area(shiftted_array, 1, X_POINTS-2);
					
					printf("MINUS SQR0 = %.4f; SQR1 = %.4f m = %d j = %d\n", s_trap[2], s_trap[m], m, j-1);
					printf("--------------------------------------------\n");
				
					break;
				}
			}
			}
		}
	}*/
	
	//max_trap = max_bubble(s_trap, num_arr_el);
	max_trap = 2*4095.0*2800/2800.0;
	printf("max_trap = %.2f\n", max_trap);
	for(i=0; i<=num_arr_el-1; i++) {
		s_trap_int[i] = (int) (4095.0*s_trap[i]/max_trap);
		if(s_trap_int[i] >= 4095) {s_trap_int[i] = 4095;}
	//	if(s_trap_int[i] < 4000) {find_pick_start_stop(a[i], &max, &min); printf("i= %d; max = %d; min = %d; s_trap = %.2f, max_trap = %.2f\n", i, max, min, s_trap[i], max_trap);}
	}

	//find_pick_start_stop(a[0], &max, &min); printf("i= %d; max = %d; min = %d; s_trap = %.2f\n", i, max, min, s_trap[1]);

	int xnum=0;
	for(j=1; j<=4096-1; j++) {
		xnum += hiztogram_trap[j];
	}
	printf("x = %d num_arr_el = %d, max_trap = %.2f\n", xnum, num_arr_el, max_trap);

	for(i=0; i<=num_arr_el-1; i++)
		for(j=0; j<=4096-2; j++)
			if( (s_trap_int[i] >= j) && (s_trap_int[i] <= (j+1)) ) {
				hiztogram_trap[j+1]++;
				break;
			}

	for(j=1; j<=4096-1; j++)
		xnum += hiztogram_trap[j];
	printf("x = %d num_arr_el = %d\n", xnum, num_arr_el);
	
	FILE *in = fopen(filename_to_save, "w+");

	if(in == NULL) {
		printf("Open error\n");
		return -1;
	}
	
	for(j=1; j<=4096-1; j++)
		fprintf(in, "%d %d\n", j, hiztogram_trap[j]);
	
	fclose(in);
	
	DenisCalc(hiztogram_trap, NULL, 2000, 2100, 1, 0, 0, ucalc);
	for(i=0; i<sizeof(ucalc)/sizeof(float); i++)
		printf("Ucalc[%d] = %.3f\n", i, ucalc[i]); //pick pos = u[2] +- u[3]; resolution = ucalc[4] +- u[7]
	
	if(ADC_NUM == 0)
		for(i=0; i<=num_arr_el-1; i++)
			square_osci[0][i] = s_trap[i];
	else if(ADC_NUM == 1)
		for(i=0; i<=num_arr_el-1; i++)
			square_osci[1][i] = s_trap[i];
	
	free(delta);
	free(s_trap);
	free(shiftted_array);
	free(k); free(b);
	free(s_trap_int);
	free(hiztogram_trap);
}

#define f(x, x0) ( -2.49991703469447e-09*(exp(-((x)-(x0))/115069.384046802) - exp(-((x)-(x0))/1.45088424300257)) - 517.224857987907*exp(-((x)-(x0))/24.9381408089486) + 2506.78894844553 )
#define NORMAL_SQUARE	(-27702.6)

double find_start_pick(int *x, char flag)
{
//	#define f(x, x0) ( -0.000725124002756408*(exp(-((x)-(x0))/27461468.6) - exp(-((x)-(x0))/1.31784519694857)) -1855.99996745509*exp(-((x)-(x0))/34.3120374841288) + 3298.00001427814 )
	#define STEP	0.001 
	#define	EPSILON	0.000001
	double x0INIT = 52.0;
	
	double S_optim, S_plus, S_minus, S_test;
	double S_plus_prev, S_minus_prev;
	double x0new = x0INIT;
	int i, j;
	int start_approx = 13;
	int end_approx = 20;
	
	int x_min_num = min_bubble_num(x, SIZEOF_DATA);
	
	if(flag == 0) {
		double baseline = (x[0] + x[1] + x[2] + x[3] + x[4]) / 5.0;
		int max, min;
		find_pick_start_stop(x, &max, &min);
		
		double normalization = -1.0*NORMAL_SQUARE/trap_area(x, max, min);
		printf("normalization = %.4f\n", normalization);
		
		for(i=0; i<=SIZEOF_DATA-1; i++) {
			x[i] = (int)( baseline - normalization*(baseline-(double)x[i]) );
		} 
	}
	
	for(i=0; i<x_min_num; i++) {
		if((x[i] <= 0.97*x[i-1]) && (x[i] <= 0.97*x[i-1])) {
			start_approx = i;
			break;
		}
	}
 
	end_approx = start_approx + 7;
	printf("start = %d; end = %d approx\n", start_approx, end_approx);
	/*
	for(i = x_min_num; i<SIZEOF_DATA(FILETYPE)-1; i++) {
		if( (x[i] >= 0.7*x[x_min_num]) && (x[i+1] >= 0.7*x[x_min_num-4]) ) {
			x0INIT = (double)(i);
			break;
		}
	}
	printf("x_min_num = %d, x0INIT = %.4f\n", x_min_num, x0INIT);
	*/
	S_optim = 0.0;
	S_test = 0.0;
	
	for(i=start_approx; i<=end_approx-1; i++) {
		S_optim += fabs(x[i] - f((double)i, x0INIT));
		S_test += fabs(x[i] - f((double)i, 53.7495));
	}
	
	S_plus = S_minus = 0.0;
	for(j=1; j <= (int)(3.0/STEP); j++) {
		S_plus_prev = S_plus;
		S_minus_prev = S_minus;
		S_plus = S_minus = 0.0;
		
		for(i=start_approx; i<=end_approx-1; i++) {
			S_plus += fabs(x[i] - f((double)i, x0INIT + j*STEP));
			S_minus += fabs(x[i] - f((double)i, x0INIT - j*STEP));
		}
		
		if(S_plus <= S_minus) {
			if(S_plus < S_optim) {
				S_optim = S_plus;
				x0new = x0INIT + STEP*j;
				if(fabs(S_plus - S_plus_prev) <= EPSILON) {
					break;
				}
			}
		}
		else {
			if(S_minus < S_optim) {
		//		printf("j=%d Soptim = %.2f, S+ = %.2f, S- = %.2f, S_test = %.2f; x0new = %.2f\n", j, S_optim, S_plus, S_minus, S_test, x0new);
				S_optim = S_minus;
				x0new = x0INIT - STEP*j;
				if(fabs(S_minus - S_minus_prev) <= EPSILON) {
					break;
				}
			}
		}
	}
	
	printf("j=%d Soptim = %.2f, S+ = %.2f, S- = %.2f, S_test = %.2f; x0new = %.2f\n", j, S_optim, S_plus, S_minus, S_test, x0new);
	
	return x0new;
}

double find_start_cftrace(int *x)
{
	int i;
	for(i=1; i<=SIZEOF_DATA-1; i++) {
		x[i] = 3000 - x[i];
	}
	
	int *CFTrace = (int *)calloc(SIZEOF_DATA, sizeof(int));
	for(i=1; i<=SIZEOF_DATA-1; i++) {
		CFTrace[i] = 0.5*x[i-1] - x[i-1-2];
	}
	
	double background = (CFTrace[5] + CFTrace[6] + CFTrace[7] + CFTrace[8] + CFTrace[9])/5.0;
	double a, b, x0 = -1.0;
	
	for(i=10; i<SIZEOF_DATA; i++) {
		if((CFTrace[i-1] >= background) && (CFTrace[i] <= background) && (CFTrace[i-1] >= background+5)) {
			printf("i = %d, background = %.2f, data[5] = %d\n", i, background, CFTrace[5]);
			a = (double)(CFTrace[i]-CFTrace[i-1]);
			b = (double)CFTrace[i] - a*i;
			
			x0 = (background - b)/a;
			break;
		}
	}
	
	free(CFTrace); 

	return x0;
}

int time_histo(int **a, int **b, const char *filename_to_save)
{
	printf("\nTIME HISTO IN WORK\n");
	
	int i, j;
	FILE *in, *in_to_time;
	double *start_a = (double *)calloc(num_arr_el, sizeof(double));
	double *start_b = (double *)calloc(num_arr_el, sizeof(double));
	double *delta_start = (double *)calloc(num_arr_el, sizeof(double));
	int *delta_start_int = (int *)calloc(num_arr_el, sizeof(int));
	double min_delta, max_delta;
	double avg_start = 0.0;
    double disp_start = 0.0;
    int *hiztogram_start = (int *)calloc(4096, sizeof(int));
    char timename[60];
	
	clock_t start_time, end_time;
	double cpu_time_used;
	
	start_time = clock();
	for(i=0; i<=num_arr_el-1; i++) {
		//start_a[i] = find_start_pick(a[i], 0);
		//start_b[i] = find_start_pick(b[i], 1);
		
		start_a[i] = find_start_cftrace(a[i]);
		start_b[i] = find_start_cftrace(b[i]);
		delta_start[i] = start_a[i]-start_b[i];
	//	printf("starta[%d]=%.4f startb[%d]=%.4f ds[%d]=%.4f\t", i, start_a[i], i, start_b[i], i, delta_start[i]);
	}
	end_time = clock();
	cpu_time_used = ((double) (end_time - start_time)) / CLOCKS_PER_SEC;
	printf("\n------------------- cpu_time_used = %e; end=%ld, start=%ld\n", cpu_time_used, end_time, start_time);
	
	for(i=1, max_delta = fabs(delta_start[0]), min_delta = delta_start[0]; i<=num_arr_el-1; i++) {
		if(delta_start[i]<min_delta) min_delta = delta_start[i];
		if(fabs(delta_start[i])>max_delta) max_delta = fabs(delta_start[i]);
	}
	printf("min start = %.2f, max_start = %.2f\n", min_delta, max_delta);
	
	/*
	for(i=0; i<num_arr_el-1; i++) {
		sprintf(timename, "./times/start%d-%.3f-%.3f", i, start_a[i], start_b[i]);
		in_to_time = fopen(timename, "w+");
		for(j=0; j<=X_POINTS-1; j++)
			fprintf(in_to_time, "%d %d %d\n", j, a[i][j], b[i][j]);
		fclose(in_to_time);
	}*/
	
	for(i=0; i<=num_arr_el-1; i++) {
		/*if((int)(100.0*fabs(delta_start[i]-min_delta))>200) hiztogram_start[200]++;
		else hiztogram_start[(int)(100.0*fabs(delta_start[i]-min_delta))]++;*/
		delta_start_int[i] = 100.0*(delta_start[i]+1.0); //(int) (4095.0*fabs(delta_start[i])/max_delta);
	//	printf("dsint[%d]=%d\t", i, delta_start_int[i]);
		if((delta_start_int[i] > 4095) || (delta_start_int[i] < 0)) { printf("d = %d, i = %d\n", delta_start_int[i], i); }
	}
	printf("\n-------------------\n");
	
	for(i=0; i<=num_arr_el-1; i++) {
		for(j=0; j<=4096-2; j++)
			if(delta_start_int[i] >= j && delta_start_int[i] <= (j+1)) {
				hiztogram_start[j+1]++;
				break;
			}
		if(delta_start_int[i]<0) {
			printf("Delta start[%d]=%d\t", i, delta_start_int[i]);
			hiztogram_start[4095]++;
		}
	}
	
	printf("\n-------------------\n");
//	for(i=0; i<=101-1; i++)
//		printf("%d %d\n", i, hiztogram_start[i]);
	
	in = fopen(filename_to_save, "w+");
	for(i=0; i<=4096-1; i++) 
		fprintf(in, "%d %d\n", i, hiztogram_start[i]);
	fclose(in);
	
	int sum_hist_start = 0;
	for(i=0; i<=4096-1; i++)
		sum_hist_start += hiztogram_start[i];
	printf("sum_hist_Start = %d\n", sum_hist_start);
	
	for(i=0; i<=num_arr_el-1; i++) {
		avg_start += delta_start[i];
		disp_start += sqr(delta_start[i]);
	}	
	
	//disp_start = sqrt( (disp_start - sqr(avg_start)/num_arr_el)/(num_arr_el-1.0) );
	avg_start = avg_start/num_arr_el;
	disp_start = sqrt( disp_start/(num_arr_el-1.0) - sqr(avg_start) );
	printf("Time start = %.2f +- %.4f\n", avg_start, disp_start);
	
	for(i=0; i<=num_arr_el-1; i++) {
		time_osci[0][i] = start_a[i];
		time_osci[1][i] = start_b[i];
	}
	
	DenisCalc(hiztogram_start, NULL, 155, 185, 1, 0, 0, ucalc);
	for(i=0; i<sizeof(ucalc)/sizeof(float); i++)
		printf("Ucalc[%d] = %.3f\n", i, ucalc[i]); //pick pos = u[2] +- u[3]; resolution = ucalc[4] +- u[7]
	
	printf("Trap pick = %.2f, trap other = %.2f\n", trap_area(hiztogram_start, 100, 190), trap_area(hiztogram_start, 210, 300));
	
	free(start_a);
	free(start_b);
	free(delta_start);
	free(delta_start_int);
	free(hiztogram_start);
}

int max_bubble_int(int *x)
{
	int i;
	int z;

	for(i=1, z = (int)fabs(x[0]); i<=X_POINTS-1; i++)
		if( ((int)fabs(x[i])) > z) {
			z = (int)fabs(x[i]);
		//	printf("z = %d\n", z);
		}
	
	return z;
}

int ampl_stats(int **a)
{
	int i;
	int *ampl = (int *) malloc(num_arr_el*sizeof(int));
	if(ampl == NULL) {
		perror("Error in ampl malloc");
		return -1;
	}
	double *ampl_dbl = (double *) malloc(num_arr_el*sizeof(double));
	if(ampl_dbl == NULL) {
		perror("Error in ampl_dbl malloc");
		return -1;
	}
	
	for(i=0; i<=num_arr_el-1; i++) {
		ampl[i] = max_bubble_int(a[i]);
		ampl_dbl[i] = (double)ampl[i];
	//	printf("Ai = %d\n", ampl[i]);
	}
	
	double avg = 0.0;
	
	printf("%.2f +- %.2f\n", avg/(num_arr_el), stats(ampl_dbl, &avg));
	
	free(ampl);
	free(ampl_dbl);
	
	return 1;
}


int read_data(char *filename, int num_start, int num_stop, int **reg1, int **reg2) 
{
	printf("start=%d, stop = %d; filename = %s\n", num_start, num_stop, filename);
	
	int fd, i;
	FILE *fd_gp_out;
	int read_size = 0;
	
	if((fd = open(filename, O_RDONLY)) < -1)
		return -1;
	if(lseek(fd, X_POINTS*sizeof(int)*num_start, SEEK_SET) < 0)
		return -2;
	 
	for(i=0; i<=num_arr_el-1; i++) {
		if(read(fd, reg1[i], X_POINTS*sizeof(int)) == -1) break;
		if(read(fd, reg2[i], X_POINTS*sizeof(int)) == -1) break;
	}
		
	close(fd);

	return 1;
}
 

int main(int argc, char **argv)
{
	int opt = 0, i, j, k;
	char *in_fname = NULL;
	char filename_to_read[256];
	int number_start = 0, number_finish = 1, ret = 0;
	FILE *out;
	int **reg1, **reg2;

	while ((opt = getopt(argc, argv, "i:n:")) != -1) {
		switch(opt) {
			case 'i':
			in_fname = optarg;
			
			if(strlen(in_fname) > 256) {
				perror("Please decrease filename\n");
				break;
			}
			
			printf("Input option value=%s\n", in_fname);
			strcpy(filename_to_read, in_fname);
			
			break;
		case 'n':
			if(isNum(optarg)) {
				number_start = atoi(optarg);
				number_finish = number_start;
			}
			else {
				sscanf(optarg, "%d-%d", &number_start, &number_finish);
			}
			
			break;
		case '?':
			printf("Missing mandatory input option\n");
			break;
		}
	}
	
	if(number_finish != 1) {
		num_arr_el = (number_finish-number_start+1)/2;
	}
	else {
		int fd;
		if((fd = open(in_fname, O_RDONLY)) < -1) {
			perror("Open file error");
			return -1;
		}
		
		num_arr_el = lseek(fd, 0, SEEK_END);
		num_arr_el = num_arr_el/(2*X_POINTS*sizeof(int));
		printf("num_arr_el = %d\n", num_arr_el);
		
		close(fd);
	}
	
	
	printf("number_start = %d, number_finish = %d , %d\n", number_start, number_finish, num_arr_el);
	
	reg1 = (int **) malloc(num_arr_el*sizeof(int *));
	reg2 = (int **) malloc(num_arr_el*sizeof(int *));
	for(i=0; i<num_arr_el; i++) {
		reg1[i] = (int *) malloc(X_POINTS*sizeof(int));
		reg2[i] = (int *) malloc(X_POINTS*sizeof(int));
	}
	
	square_osci[0] = (double *) malloc(num_arr_el*sizeof(double));
	square_osci[1] = (double *) malloc(num_arr_el*sizeof(double));
	time_osci[0] = (double *) malloc(num_arr_el*sizeof(double));
	time_osci[1] = (double *) malloc(num_arr_el*sizeof(double));
	
	
	if((ret = read_data(filename_to_read, number_start, number_finish, reg1, reg2)) != 1) {
		printf("Error in read_data func ret = %d; errno = %d\n", ret, errno);
		return 1;
	}

	
    
    out = fopen("/home/das/job/miniPC/last_spk.txt", "w+");
    if(out == NULL) {
		printf("Open out error\n");
		return -1;
	}
	
	for(i=0; i<=num_arr_el-1; i++) {
		for(j=0; j<=X_POINTS-1; j++)
			fprintf(out, "%d %d\n", j+2*i*X_POINTS, reg1[i][j]);
		for(j=0; j<=X_POINTS-1; j++)
			fprintf(out, "%d %d\n", j+(2*i+1)*X_POINTS, reg2[i][j]);
	}

    fclose(out);
	
	printf("0.1\n");
	/*int max, min;
	for(i=0; i<=num_arr_el-1; i++) {
		find_pick_start_stop(reg1[i], &max, &min);
		printf("\nnum pick for reg1[%d]= %d, %d\n", i, max, min);
	}*/
    
	energy_histo(reg1, "./histo_en/en1", 0);
//	ampl_stats(reg1);
//	printf("0\n");
	energy_histo(reg2, "./histo_en/en2", 1);
//	printf("1\n");
	time_histo(reg1, reg2, "./histo_t/t");
	/*
	out = fopen("/home/das/job/miniPC/calc_results", "w+");
    if(out == NULL) {
		printf("Open out error\n");
		return -1;
	}
    for(i=0; i<=ARR_NUMS-1; i++)
        fprintf(out, "%.2f\t%.2f\n%8.2f\t%8.2f\n", square_osci[0][i], square_osci[1][i], time_osci[0][i], time_osci[1][i]);
    close(out);
	*/
	
	/*for(i=0; i<=num_arr_el-1; i++)
		printf("%.3f %.3f diff = %e\n", second_derivative_zero(reg1[i]), second_derivative_zero(reg2[i]), fabs(second_derivative_zero(reg1[i])-second_derivative_zero(reg2[i])));*/
	//second_derivative_zero(reg1[0]);
	
	return 0;
}
