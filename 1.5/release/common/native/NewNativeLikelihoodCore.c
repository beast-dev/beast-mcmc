#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "NewNativeLikelihoodCore.h"

#define MATRIX_SIZE (STATE_COUNT + 1) * STATE_COUNT
#if (STATE_COUNT==4)
#define IS_NUCLEOTIDES
#endif

int nodeCount;
int stateTipCount;
int patternCount;
int partialsSize;
int matrixCount;

double* cMatrix;
double* storedCMatrix;
double* eigenValues;
double* storedEigenValues;

double* frequencies;
double* storedFrequencies;
double* categoryProportions;
double* storedCategoryProportions;
double* categoryRates;
double* storedCategoryRates;

double* integrationTmp;

double*** partials;
int** states;
double*** matrices;

int* currentMatricesIndices;
int* storedMatricesIndices;
int* currentPartialsIndices;
int* storedPartialsIndices;

void updateStatesStates(int nodeIndex1, int nodeIndex2, int nodeIndex3);
void updateStatesPartials(int nodeIndex1, int nodeIndex2, int nodeIndex3);
void updatePartialsPartials(int nodeIndex1, int nodeIndex2, int nodeIndex3);

void printArray(char* name, double *array, int length) {
	fprintf(stdout, "%s:", name);
	for (int i = 0; i < length; i++) {
		fprintf(stdout, " %f", array[i]);
	}
	fprintf(stdout, "\n");
}

/*
 * Class:     dr_evomodel_newtreelikelihood_NativeLikelihoodCore
 * Method:    initialize
 * Signature: (III)V
 */
JNIEXPORT void JNICALL Java_dr_evomodel_newtreelikelihood_NativeLikelihoodCore_initialize
  (JNIEnv *env, jobject obj, jint inNodeCount, jint inStateTipCount, jint inPatternCount, jint inMatrixCount)
{	
	nodeCount = inNodeCount;
	stateTipCount = inStateTipCount;
	patternCount = inPatternCount;
	matrixCount = inMatrixCount;
	
	cMatrix = (double *)malloc(sizeof(double) * STATE_COUNT * STATE_COUNT * STATE_COUNT);
	storedCMatrix = (double *)malloc(sizeof(double) * STATE_COUNT * STATE_COUNT * STATE_COUNT);

	eigenValues = (double *)malloc(sizeof(double) * STATE_COUNT);
	storedEigenValues = (double *)malloc(sizeof(double) * STATE_COUNT);

	frequencies = (double *)malloc(sizeof(double) * STATE_COUNT);
	storedFrequencies = (double *)malloc(sizeof(double) * STATE_COUNT);

	categoryRates = (double *)malloc(sizeof(double) * matrixCount);
	storedCategoryRates = (double *)malloc(sizeof(double) * matrixCount);

	categoryProportions = (double *)malloc(sizeof(double) * matrixCount);
	storedCategoryProportions = (double *)malloc(sizeof(double) * matrixCount);

	// a temporary array used in calculating log likelihoods
	integrationTmp = (double *)malloc(sizeof(double) * patternCount * STATE_COUNT);

	partialsSize = patternCount * STATE_COUNT * matrixCount;

	partials = (double ***)malloc(sizeof(double**) * 2);
	partials[0] = (double **)malloc(sizeof(double*) * nodeCount);
	partials[1] = (double **)malloc(sizeof(double*) * nodeCount);

	states = (int **)malloc(sizeof(int*) * nodeCount);
	
	for (int i = 0; i < nodeCount; i++) {
		partials[0][i] = (double *)malloc(sizeof(double) * partialsSize);
		partials[1][i] = (double *)malloc(sizeof(double) * partialsSize);
		states[i] = (int *)malloc(sizeof(int) * patternCount * matrixCount);
	}

  	currentMatricesIndices = (int *)malloc(sizeof(int) * nodeCount);
  	memset(currentMatricesIndices, 0, sizeof(int) * nodeCount);
  	storedMatricesIndices = (int *)malloc(sizeof(int) * nodeCount);

  	currentPartialsIndices = (int *)malloc(sizeof(int) * nodeCount);
  	memset(currentPartialsIndices, 0, sizeof(int) * nodeCount);
  	storedPartialsIndices = (int *)malloc(sizeof(int) * nodeCount);
       
	matrices = (double ***)malloc(sizeof(double**) * 2);
	matrices[0] = (double **)malloc(sizeof(double*) * nodeCount);
	matrices[1] = (double **)malloc(sizeof(double*) * nodeCount);
	for (int i = 0; i < nodeCount; i++) {
		matrices[0][i] = (double *)malloc(sizeof(double) * matrixCount * MATRIX_SIZE);
		matrices[1][i] = (double *)malloc(sizeof(double) * matrixCount * MATRIX_SIZE);
	}
	
}

/*
 * Class:     dr_evomodel_newtreelikelihood_NativeLikelihoodCore
 * Method:    freeNativeMemory
 * Signature: ()V
 */
JNIEXPORT void JNICALL Java_dr_evomodel_newtreelikelihood_NativeLikelihoodCore_freeNativeMemory
  (JNIEnv *env, jobject obj) 
{
	// free all that stuff...
}

/*
 * Class:     dr_evomodel_newtreelikelihood_NativeLikelihoodCore
 * Method:    setTipPartials
 * Signature: (I[D)V
 */
JNIEXPORT void JNICALL Java_dr_evomodel_newtreelikelihood_NativeLikelihoodCore_setTipPartials
  (JNIEnv *env, jobject obj, jint tipIndex, jdoubleArray inTipPartials)
{
	jdouble *tipPartials = (jdouble*)(*env)->GetPrimitiveArrayCritical(env, inTipPartials, 0);

	int k = 0;
	for (int i = 0; i < matrixCount; i++) {
		memcpy(partials[0][tipIndex] + k, tipPartials, sizeof(double) * patternCount * STATE_COUNT);
		k += patternCount * STATE_COUNT;
	}
	
	(*env)->ReleasePrimitiveArrayCritical(env, inTipPartials, tipPartials, JNI_ABORT);
}

/*
 * Class:     dr_evomodel_newtreelikelihood_NativeLikelihoodCore
 * Method:    setTipStates
 * Signature: (I[I)V
 */
JNIEXPORT void JNICALL Java_dr_evomodel_newtreelikelihood_NativeLikelihoodCore_setTipStates
  (JNIEnv *env, jobject obj, jint tipIndex, jintArray inTipStates)
{
	jint *tipStates = (jint*)(*env)->GetPrimitiveArrayCritical(env, inTipStates, 0);

	int k = 0;
	for (int i = 0; i < matrixCount; i++) {
		for (int j = 0; j < patternCount; j++) {
			states[tipIndex][k] = (tipStates[j] < STATE_COUNT ? tipStates[j] : STATE_COUNT);
			k++;
		}
	}
	
	(*env)->ReleasePrimitiveArrayCritical(env, inTipStates, tipStates, JNI_ABORT);
}

/*
 * Class:     dr_evomodel_newtreelikelihood_NativeLikelihoodCore
 * Method:    updateRootFrequencies
 * Signature: ([D)V
 */
JNIEXPORT void JNICALL Java_dr_evomodel_newtreelikelihood_NativeLikelihoodCore_updateRootFrequencies
  (JNIEnv *env, jobject obj, jdoubleArray inFrequencies) 
{
    (*env)->GetDoubleArrayRegion(env, inFrequencies, 0, STATE_COUNT, frequencies);
}

/*
 * Class:     dr_evomodel_newtreelikelihood_NativeLikelihoodCore
 * Method:    updateEigenDecomposition
 * Signature: ([[D[[D[D)V
 */
JNIEXPORT void JNICALL Java_dr_evomodel_newtreelikelihood_NativeLikelihoodCore_updateEigenDecomposition
  (JNIEnv *env, jobject obj, jobjectArray inEigenVectors, jobjectArray inInvEigenVectors, jdoubleArray inEigenValues)
{
	double Evec[STATE_COUNT][STATE_COUNT];
	double Ievc[STATE_COUNT][STATE_COUNT];
	
	for (int i = 0; i < STATE_COUNT; i++) {
		jdoubleArray oneDim = (jdoubleArray)(*env)->GetObjectArrayElement(env, inEigenVectors, i);
		jdouble *element = (*env)->GetDoubleArrayElements(env, oneDim, 0);
		for(int j = 0; j < STATE_COUNT; j++) {
			Evec[i][j]= element[j];
		}
	}

	for (int i = 0; i < STATE_COUNT; i++) {
		jdoubleArray oneDim = (jdoubleArray)(*env)->GetObjectArrayElement(env, inInvEigenVectors, i);
		jdouble *element = (*env)->GetDoubleArrayElements(env, oneDim, 0);
		for(int j = 0; j < STATE_COUNT; j++) {
			Ievc[i][j]= element[j];
		}
	}
	
	int l =0;
	for (int i = 0; i < STATE_COUNT; i++) {
		for (int j = 0; j < STATE_COUNT; j++) {
			for (int k = 0; k < STATE_COUNT; k++) {
				cMatrix[l] = Evec[i][k] * Ievc[k][j];
				l++;
			}
		}
	}

    (*env)->GetDoubleArrayRegion(env, inEigenValues, 0, STATE_COUNT, eigenValues);
}

/*
 * Class:     dr_evomodel_newtreelikelihood_NativeLikelihoodCore
 * Method:    updateCategoryRates
 * Signature: ([D)V
 */
JNIEXPORT void JNICALL Java_dr_evomodel_newtreelikelihood_NativeLikelihoodCore_updateCategoryRates
  (JNIEnv *env, jobject obj, jdoubleArray inRates)
{
    (*env)->GetDoubleArrayRegion(env, inRates, 0, matrixCount, categoryRates);
}

/*
 * Class:     dr_evomodel_newtreelikelihood_NativeLikelihoodCore
 * Method:    updateCategoryProportions
 * Signature: ([D)V
 */
JNIEXPORT void JNICALL Java_dr_evomodel_newtreelikelihood_NativeLikelihoodCore_updateCategoryProportions
  (JNIEnv *env, jobject obj, jdoubleArray inProportions)
{
    (*env)->GetDoubleArrayRegion(env, inProportions, 0, matrixCount, categoryProportions);
}

/*
 * Class:     dr_evomodel_newtreelikelihood_NativeLikelihoodCore
 * Method:    updateMatrices
 * Signature: ([I[DI)V
 */
JNIEXPORT void JNICALL Java_dr_evomodel_newtreelikelihood_NativeLikelihoodCore_updateMatrices
  (JNIEnv *env, jobject obj, jintArray inNodeIndices, jdoubleArray inBranchLengths, jint count)
{
	jint *nodeIndices = (jint*)(*env)->GetPrimitiveArrayCritical(env, inNodeIndices, 0);
	jdouble *branchLengths = (jdouble*)(*env)->GetPrimitiveArrayCritical(env, inBranchLengths, 0);
	        	
	static double tmp[STATE_COUNT];
	
	for (int x = 0; x < count; x++) {
		int nodeIndex = nodeIndices[x];
		double branchLength = branchLengths[x];
		
		currentMatricesIndices[nodeIndex] = 1 - currentMatricesIndices[nodeIndex];

        int n = 0;
        for (int l = 0; l < matrixCount; l++) {
            for (int i = 0; i < STATE_COUNT; i++) {
                tmp[i] =  exp(eigenValues[i] * branchLength * categoryRates[l]);
            }

            int m = 0;
            for (int i = 0; i < STATE_COUNT; i++) {
                for (int j = 0; j < STATE_COUNT; j++) {
                    double sum = 0.0;
                    for (int k = 0; k < STATE_COUNT; k++) {
                        sum += cMatrix[m] * tmp[k];
                        m++;
                    }
                    matrices[currentMatricesIndices[nodeIndex]][nodeIndex][n] = sum;
                    n++;
                }
                matrices[currentMatricesIndices[nodeIndex]][nodeIndex][n] = 1.0;
                n++;
            }
        }
        

	}
	
	(*env)->ReleasePrimitiveArrayCritical(env, inBranchLengths, branchLengths, JNI_ABORT);
	(*env)->ReleasePrimitiveArrayCritical(env, inNodeIndices, nodeIndices, JNI_ABORT);
}

/*
 * Class:     dr_evomodel_newtreelikelihood_NativeLikelihoodCore
 * Method:    updatePartials
 * Signature: ([I[II)V
 */
JNIEXPORT void JNICALL Java_dr_evomodel_newtreelikelihood_NativeLikelihoodCore_updatePartials
  (JNIEnv *env, jobject obj, jintArray inOperations, jintArray inDependencies, jint operationCount)
{
	jint *operations = (jint*)(*env)->GetPrimitiveArrayCritical(env, inOperations, 0);
	/* not required:
	jint *dependencies = (jint*)(*env)->GetPrimitiveArrayCritical(env, inDependencies, 0);
	*/

    int x = 0;
	for (int op = 0; op < operationCount; op++) {
		int nodeIndex1 = operations[x];
		x++;
		int nodeIndex2 = operations[x];
		x++;
		int nodeIndex3 = operations[x];
		x++;
		currentPartialsIndices[nodeIndex3] = 1 - currentPartialsIndices[nodeIndex3];

		if (nodeIndex1 < stateTipCount) {
			if (nodeIndex2 < stateTipCount) {
				updateStatesStates(nodeIndex1, nodeIndex2, nodeIndex3);
			} else {
				updateStatesPartials(nodeIndex1, nodeIndex2, nodeIndex3);
			}
		} else {
			if (nodeIndex2 < stateTipCount) {
				updateStatesPartials(nodeIndex2, nodeIndex1, nodeIndex3);
			} else {
				updatePartialsPartials(nodeIndex1, nodeIndex2, nodeIndex3);
			}
		}
	}
	/* not required:
	(*env)->ReleasePrimitiveArrayCritical(env, inDependencies, dependencies, JNI_ABORT);
	*/
	(*env)->ReleasePrimitiveArrayCritical(env, inOperations, operations, JNI_ABORT);
}

/*
 * Calculates partial likelihoods at a node when both children have states.
 */
void updateStatesStates(int nodeIndex1, int nodeIndex2, int nodeIndex3)
{
	double* matrices1 = matrices[currentMatricesIndices[nodeIndex1]][nodeIndex1];
	double* matrices2 = matrices[currentMatricesIndices[nodeIndex2]][nodeIndex2];

	int* states1 = states[nodeIndex1];
	int* states2 = states[nodeIndex2];

	double* partials3 = partials[currentPartialsIndices[nodeIndex3]][nodeIndex3];

    #ifdef IS_NUCLEOTIDES

	int v = 0;
	for (int l = 0; l < matrixCount; l++) {

		for (int k = 0; k < patternCount; k++) {

			int state1 = states1[k];
			int state2 = states2[k];

			int w = l * MATRIX_SIZE;

			partials3[v] = matrices1[w + state1] * matrices2[w + state2];
			v++;	w += (STATE_COUNT + 1);
			partials3[v] = matrices1[w + state1] * matrices2[w + state2];
			v++;	w += (STATE_COUNT + 1);
			partials3[v] = matrices1[w + state1] * matrices2[w + state2];
			v++;	w += (STATE_COUNT + 1);
			partials3[v] = matrices1[w + state1] * matrices2[w + state2];
			v++;	w += (STATE_COUNT + 1);

		}
	}

	#else

	int v = 0;
	for (int l = 0; l < matrixCount; l++) {

		for (int k = 0; k < patternCount; k++) {

			int state1 = states1[k];
			int state2 = states2[k];

			int w = l * MATRIX_SIZE;

			for (int i = 0; i < STATE_COUNT; i++) {

				partials3[v] = matrices1[w + state1] * matrices2[w + state2];

				v++;
				w += (STATE_COUNT + 1);
			}

		}
	}
	#endif
}

/*
 * Calculates partial likelihoods at a node when one child has states and one has partials.
 */
void updateStatesPartials(int nodeIndex1, int nodeIndex2, int nodeIndex3)
{
	double* matrices1 = matrices[currentMatricesIndices[nodeIndex1]][nodeIndex1];
	double* matrices2 = matrices[currentMatricesIndices[nodeIndex2]][nodeIndex2];

	int* states1 = states[nodeIndex1];
	double* partials2 = partials[currentPartialsIndices[nodeIndex2]][nodeIndex2];

	double* partials3 = partials[currentPartialsIndices[nodeIndex3]][nodeIndex3];

    #ifdef IS_NUCLEOTIDES

	int u = 0;
	int v = 0;

	for (int l = 0; l < matrixCount; l++) {
		for (int k = 0; k < patternCount; k++) {

			int state1 = states1[k];

			int w = l * MATRIX_SIZE;

			partials3[u] = matrices1[w + state1];

			double sum = matrices2[w] * partials2[v]; w++;
			sum +=	matrices2[w] * partials2[v + 1]; w++;
			sum +=	matrices2[w] * partials2[v + 2]; w++;
			sum +=	matrices2[w] * partials2[v + 3]; w++;
			w++; // increment for the extra column at the end
			partials3[u] *= sum;	u++;

			partials3[u] = matrices1[w + state1];

			sum = matrices2[w] * partials2[v]; w++;
			sum +=	matrices2[w] * partials2[v + 1]; w++;
			sum +=	matrices2[w] * partials2[v + 2]; w++;
			sum +=	matrices2[w] * partials2[v + 3]; w++;
			w++; // increment for the extra column at the end
			partials3[u] *= sum;	u++;

			partials3[u] = matrices1[w + state1];

			sum = matrices2[w] * partials2[v]; w++;
			sum +=	matrices2[w] * partials2[v + 1]; w++;
			sum +=	matrices2[w] * partials2[v + 2]; w++;
			sum +=	matrices2[w] * partials2[v + 3]; w++;
			w++; // increment for the extra column at the end
			partials3[u] *= sum;	u++;

			partials3[u] = matrices1[w + state1];

			sum = matrices2[w] * partials2[v]; w++;
			sum +=	matrices2[w] * partials2[v + 1]; w++;
			sum +=	matrices2[w] * partials2[v + 2]; w++;
			sum +=	matrices2[w] * partials2[v + 3]; w++;
			w++; // increment for the extra column at the end
			partials3[u] *= sum;	u++;

			v += 4;

		}
	}

	#else
	int u = 0;
	int v = 0;

	for (int l = 0; l < matrixCount; l++) {
		for (int k = 0; k < patternCount; k++) {

			int state1 = states1[k];

			int w = l * MATRIX_SIZE;

			for (int i = 0; i < STATE_COUNT; i++) {

				double tmp = matrices1[w + state1];

				double sum = 0.0;
				for (int j = 0; j < STATE_COUNT; j++) {
					sum += matrices2[w] * partials2[v + j];
					w++;
				}

				// increment for the extra column at the end
				w++;

				partials3[u] = tmp * sum;
				u++;
			}

			v += STATE_COUNT;
		}
	}
	#endif
}

void updatePartialsPartials(int nodeIndex1, int nodeIndex2, int nodeIndex3)
{
	double* matrices1 = matrices[currentMatricesIndices[nodeIndex1]][nodeIndex1];
	double* matrices2 = matrices[currentMatricesIndices[nodeIndex2]][nodeIndex2];

	double* partials1 = partials[currentPartialsIndices[nodeIndex1]][nodeIndex1];
	double* partials2 = partials[currentPartialsIndices[nodeIndex2]][nodeIndex2];

	double* partials3 = partials[currentPartialsIndices[nodeIndex3]][nodeIndex3];

	/* fprintf(stdout, "*** operation %d: %d, %d -> %d\n", op, nodeIndex1, nodeIndex2, nodeIndex3); */

	double sum1, sum2;

    #ifdef IS_NUCLEOTIDES

	int u = 0;
	int v = 0;

	for (int l = 0; l < matrixCount; l++) {
		for (int k = 0; k < patternCount; k++) {

			int w = l * MATRIX_SIZE;

			sum1 = matrices1[w] * partials1[v];
			sum2 = matrices2[w] * partials2[v]; w++;
			sum1 += matrices1[w] * partials1[v + 1];
			sum2 += matrices2[w] * partials2[v + 1]; w++;
			sum1 += matrices1[w] * partials1[v + 2];
			sum2 += matrices2[w] * partials2[v + 2]; w++;
			sum1 += matrices1[w] * partials1[v + 3];
			sum2 += matrices2[w] * partials2[v + 3]; w++;
			w++; // increment for the extra column at the end
			partials3[u] = sum1 * sum2; u++;

			sum1 = matrices1[w] * partials1[v];
			sum2 = matrices2[w] * partials2[v]; w++;
			sum1 += matrices1[w] * partials1[v + 1];
			sum2 += matrices2[w] * partials2[v + 1]; w++;
			sum1 += matrices1[w] * partials1[v + 2];
			sum2 += matrices2[w] * partials2[v + 2]; w++;
			sum1 += matrices1[w] * partials1[v + 3];
			sum2 += matrices2[w] * partials2[v + 3]; w++;
			w++; // increment for the extra column at the end
			partials3[u] = sum1 * sum2; u++;

			sum1 = matrices1[w] * partials1[v];
			sum2 = matrices2[w] * partials2[v]; w++;
			sum1 += matrices1[w] * partials1[v + 1];
			sum2 += matrices2[w] * partials2[v + 1]; w++;
			sum1 += matrices1[w] * partials1[v + 2];
			sum2 += matrices2[w] * partials2[v + 2]; w++;
			sum1 += matrices1[w] * partials1[v + 3];
			sum2 += matrices2[w] * partials2[v + 3]; w++;
			w++; // increment for the extra column at the end
			partials3[u] = sum1 * sum2; u++;

			sum1 = matrices1[w] * partials1[v];
			sum2 = matrices2[w] * partials2[v]; w++;
			sum1 += matrices1[w] * partials1[v + 1];
			sum2 += matrices2[w] * partials2[v + 1]; w++;
			sum1 += matrices1[w] * partials1[v + 2];
			sum2 += matrices2[w] * partials2[v + 2]; w++;
			sum1 += matrices1[w] * partials1[v + 3];
			sum2 += matrices2[w] * partials2[v + 3]; w++;
			w++; // increment for the extra column at the end
			partials3[u] = sum1 * sum2; u++;

			v += 4;

		}
	}

	#else

	int u = 0;
	int v = 0;

	for (int l = 0; l < matrixCount; l++) {

		for (int k = 0; k < patternCount; k++) {

			int w = l * MATRIX_SIZE;

			for (int i = 0; i < STATE_COUNT; i++) {

				sum1 = sum2 = 0.0;

				for (int j = 0; j < STATE_COUNT; j++) {
					sum1 += matrices1[w] * partials1[v + j];
					sum2 += matrices2[w] * partials2[v + j];
					w++;
				}

				// increment for the extra column at the end
				w++;

				partials3[u] = sum1 * sum2;
				
				u++;
			}		
			v += STATE_COUNT;
		}
	}

	#endif
}

/*
 * Class:     dr_evomodel_newtreelikelihood_NativeLikelihoodCore
 * Method:    calculateLogLikelihoods
 * Signature: (I[D)V
 */

JNIEXPORT void JNICALL Java_dr_evomodel_newtreelikelihood_NativeLikelihoodCore_calculateLogLikelihoods
  (JNIEnv *env, jobject obj, jint rootNodeIndex, jdoubleArray outLogLikelihoods)
{
	jdouble *logLikelihoods = (jdouble*)(*env)->GetPrimitiveArrayCritical(env, outLogLikelihoods, 0);

	double* rootPartials = partials[currentPartialsIndices[rootNodeIndex]][rootNodeIndex];
	
	int u = 0;
	int v = 0;
	for (int k = 0; k < patternCount; k++) {

		for (int i = 0; i < STATE_COUNT; i++) {

			integrationTmp[u] = rootPartials[v] * categoryProportions[0];
			u++;
			v++;
		}
	}


	for (int l = 1; l < matrixCount; l++) {
		u = 0;

		for (int k = 0; k < patternCount; k++) {

			for (int i = 0; i < STATE_COUNT; i++) {

				integrationTmp[u] += rootPartials[v] * categoryProportions[l];
				u++;
				v++;
			}
		}
	}

	u = 0;
	for (int k = 0; k < patternCount; k++) {

		double sum = 0.0;
		for (int i = 0; i < STATE_COUNT; i++) {

			sum += frequencies[i] * integrationTmp[u];
			u++;
		}
		logLikelihoods[k] = log(sum);
	}
	
	(*env)->ReleasePrimitiveArrayCritical(env, outLogLikelihoods, logLikelihoods, JNI_ABORT);
}

/*
 * Class:     dr_evomodel_newtreelikelihood_NativeLikelihoodCore
 * Method:    storeState
 * Signature: ()V
 */
JNIEXPORT void JNICALL Java_dr_evomodel_newtreelikelihood_NativeLikelihoodCore_storeState
  (JNIEnv *env, jobject obj)
{
	memcpy(storedCMatrix, cMatrix, sizeof(double) * STATE_COUNT * STATE_COUNT * STATE_COUNT);
	memcpy(storedEigenValues, eigenValues, sizeof(double) * STATE_COUNT);

	memcpy(storedFrequencies, frequencies, sizeof(double) * STATE_COUNT);
	memcpy(storedCategoryRates, categoryRates, sizeof(double) * matrixCount);
	memcpy(storedCategoryProportions, categoryProportions, sizeof(double) * matrixCount);
	
	memcpy(storedMatricesIndices, currentMatricesIndices, sizeof(int) * nodeCount);
	memcpy(storedPartialsIndices, currentPartialsIndices, sizeof(int) * nodeCount);
}

/*
 * Class:     dr_evomodel_newtreelikelihood_NativeLikelihoodCore
 * Method:    restoreState
 * Signature: ()V
 */
JNIEXPORT void JNICALL Java_dr_evomodel_newtreelikelihood_NativeLikelihoodCore_restoreState
  (JNIEnv *env, jobject obj)
{
	// Rather than copying the stored stuff back, just swap the pointers...
	double* tmp = cMatrix;
	cMatrix = storedCMatrix;
	storedCMatrix = tmp;

	tmp = eigenValues;
	eigenValues = storedEigenValues;
	storedEigenValues = tmp;

	tmp = frequencies;
	frequencies = storedFrequencies;
	storedFrequencies = tmp;

	tmp = categoryRates;
	categoryRates = storedCategoryRates;
	storedCategoryRates = tmp;

	tmp = categoryProportions;
	categoryProportions = storedCategoryProportions;
	storedCategoryProportions = tmp;

	int* tmp2 = currentMatricesIndices;
	currentMatricesIndices = storedMatricesIndices;
	storedMatricesIndices = tmp2;

	tmp2 = currentPartialsIndices;
	currentPartialsIndices = storedPartialsIndices;
	storedPartialsIndices = tmp2;
}

