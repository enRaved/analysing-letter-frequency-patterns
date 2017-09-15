#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <float.h>

#define N 2
#define M 27
#define T 50000

static double A[2][2] = {
					{0.47468, 0.52532},
					{0.51656, 0.48344}
};

static double B[N][M] = {
					{
						0.03735, 0.03408, 0.03455, 0.03828, 0.03782, 
						0.03922, 0.03688, 0.03408, 0.03875, 0.04062, 
						0.03735, 0.03968, 0.03548, 0.03735, 0.04062, 
						0.03595, 0.03641, 0.03408, 0.04062, 0.03548, 
						0.03922, 0.04062, 0.03455, 0.03595, 0.03408, 
						0.03408, 0.03688
					},
					{
						0.03909, 0.03537, 0.03537, 0.03909, 0.03583, 
						0.03630, 0.04048, 0.03537, 0.03816, 0.03909, 
						0.03490, 0.03723, 0.03537, 0.03909, 0.03397, 
						0.03397, 0.03816, 0.03676, 0.04048, 0.03443, 
						0.03537, 0.03955, 0.03816, 0.03723, 0.03769, 
						0.03955, 0.03397
					}
};



static double pi[2] = {0.51316, 0.48684};
static int maxIters;

static int observationSequence[T];
static double C[T];
static double alpha[T][N];
static double beta[T][N];
static double mygamma[T][N];
static double digamma[T][N][N];

int iters = 0;
double oldLogProb = -DBL_MAX, logProb;

void computeAlphaPass();
void computeBetaPass();
void computeGammaDigamma();
void reestimate();

void printGreekLetter(double [][N]);

int convertLettersToNumbers(char c){
	printf("%c\n", c);
	if(c != ' ')
		return c - 'a';
	else
		return 26;
}



int main(int argc, char *argv[]){

FILE *fr;
fr = fopen("brownfile.txt", "r");

	int k = 0;
	int c;

	
	//T = (long) atoi(argv[2]);
	maxIters = atoi(argv[3]);
	int seed = atoi(argv[4]);



//read from file to observation sequence
	if(fr){
		while((c=fgetc(fr))!=EOF && k < T){

 			observationSequence[k] = convertLettersToNumbers(c);
 			k++;
 		}
	}

	fclose(fr);
	//print observation sequence
	for (int i = 0; i < T; ++i)
	{
		printf("%d\t%d\n", i, observationSequence[i]);
	}


while(iters < maxIters || logProb > oldLogProb){

	computeAlphaPass();
	computeBetaPass();
	computeGammaDigamma();

	reestimate();

	// 6. Compute Compute log[P (O | λ)]
	logProb = 0;
	for(int i = 0; i < T; i++){
		logProb = logProb + log (C[i]);
	}
	logProb = -logProb;

	//To iterate or not to iterate, that is the question. . .
	

	iters = iters + 1;

	printf("Interation %d: %f\n", iters, logProb);

	//if(iters < maxIters && logProb > oldLogProb){
	if(iters < maxIters){
		oldLogProb = logProb;
		
	}
	else{
		//print pi
		// printf("\n Final Pi");
		// for(int i = 0; i < N; i++){
		// 	printf("\t%f", pi[i]);
		// }
		// //print A
		// printf("\nFinal A:");
		// for(int i = 0; i < N; i++){
		// 	printf("\n");
		// 	for(int j = 0; j < N; j++){
		// 		printf("%f\t", A[i][j]);
		// 	}
		// }
		printf("\n");
		//print B
		printf("\nFinal B:");
		for(int i = 0; i < N; i++){
			printf("\n");
			for(int j = 0; j < M; j++){
				printf("%f ", B[i][j]);
			}
		}
		printf("\n");
		break;
	}


} //while ends here



	//printGreekLetter(alpha);
	//printGreekLetter(beta);
	//printGreekLetter(mygamma);

}

void computeAlphaPass(){

	
	// compute α0(i)
		C[0] = 0;
		for(int i = 0; i < N; i++){
		alpha[0][i] = pi[i] * B[i][observationSequence[0]];
		C[0] = C[0] + alpha[0][i];
	}

	//scale the α0(i)
	C[0] = 1/C[0];
	for(int i = 0; i < N; i++){
		alpha[0][i] = C[0] * alpha[0][i];
	}

	//compute αt(i)
	for(int t = 1; t < T; t++){
		C[t] = 0;
		for(int i = 0; i < N; i++){
			alpha[t][i] = 0;
			for(int j = 0; j < N; j++){
				alpha[t][i] = alpha[t][i] + alpha[t-1][j] * A[j][i];
			}
			alpha[t][i] = alpha[t][i] * B[i][observationSequence[t]];
			C[t] = C[t] + alpha[t][i];
		}
		// scale αt(i)
		C[t] = 1/C[t];
		for(int i = 0; i < N; i++){
			alpha[t][i] = C[t] * alpha[t][i];
		}
	}

}

void computeBetaPass(){
	for(int i = 0; i < N; i++){
		beta[T-1][i] = C[T-1];
	}

	// β-pass
	for(int t = (T - 2); t >= 0; t--){
		for(int i = 0; i < N; i++){
			beta[t][i] = 0;
			for(int j = 0; j < N; j++){
				beta[t][i] = beta[t][i] + A[i][j] * B[j][observationSequence[t+1]] * beta[t+1][j];
			}
			// scale βt(i) with same scale factor as αt(i)
			beta[t][i] = C[t] * beta[t][i];
		}
	}
}

void computeGammaDigamma(){
	double numer, denom;
		for(int t = 0; t < T-1; t++){
		denom = 0;
		for(int i = 0; i < N; i++){
			for(int j = 0; j < N; j++){
				denom = denom + alpha[t][i] * A[i][j] * B[j][observationSequence[t+1]] * beta[t+1][j];
			}
		}
		for(int i = 0; i < N; i++){
			mygamma[t][i] = 0;

			for(int j = 0; j < N; j++){

				digamma[t][i][j] = ( alpha[t][i] * A[i][j] * B[j][observationSequence[t+1]] * beta[t+1][j] ) / denom;
				mygamma[t][i] = mygamma[t][i] + digamma[t][i][j];

			}
		}
	}

	// Special case for γT−1(i)
	denom = 0;
	for(int i = 0; i < N; i++){
		denom = denom + alpha[T-1][i];
	}
	for(int i = 0; i < N; i++){
		mygamma[T-1][i] = alpha[T-1][i]/denom;
	}
}

void reestimate(){
	// 5. Re-estimate A, B and π

	//re-estimate π
	for(int i = 0; i < N; i++){
		pi[i] = mygamma[0][i];
	}

	//re-estimate A
	double numer, denom;
	for(int i = 0; i < N; i++){
		for(int j = 0; j < N; j++){
			numer = 0;
			denom = 0;
			for(int t  = 0; t < T-1; t++){
				numer = numer + digamma[t][i][j];
				denom = denom + mygamma[t][i];
			}
			A[i][j] = numer/denom;
		}
	}

	//re-estimate B
	for(int i = 0; i < N; i++){
		for(int j = 0; j < M; j++){
			numer = 0;
			denom = 0;
			for(int t = 0; t < T; t++){
				if(observationSequence[t]==j){
					numer = numer + mygamma[t][i];
				}
				denom = denom + mygamma[t][i];
			}
			B[i][j] = numer/denom;
		}
	}
}

void printGreekLetter(double greekletter[][N]){
//PRINT STATEMENT REMOVE
//print alpha, beta and gamma pass
printf("Pass\n");
for(int t = 0; t < T; t++){
	for(int i = 0; i < N; i++){
			printf("%lf\t", greekletter[t][i]);
		}
		printf("\n");
}
}
