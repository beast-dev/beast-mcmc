#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "NewNativeLikelihoodCore.h"

#define STATE_COUNT 4
#define MATRIX_SIZE STATE_COUNT * STATE_COUNT

int stateCount;
int nodeCount;
int patternCount;
int partialsSize;
int matrixSize;
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
double*** matrices;

int* currentMatricesIndices;
int* storedMatricesIndices;
int* currentPartialsIndices;
int* storedPartialsIndices;

/*
 * Class:     dr_evomodel_newtreelikelihood_NativeLikelihoodCore
 * Method:    initialize
 * Signature: (III)V
 */
JNIEXPORT void JNICALL Java_dr_evomodel_newtreelikelihood_NativeLikelihoodCore_initialize
  (JNIEnv *env, jobject obj, jint inNodeCount, jint inPatternCount, jint inMatrixCount)
{
	stateCount = STATE_COUNT;
	
	nodeCount = inNodeCount;
	patternCount = inPatternCount;
	matrixCount = inMatrixCount;
	
	cMatrix = (double *)malloc(sizeof(double) * stateCount * stateCount * stateCount);
	storedCMatrix = (double *)malloc(sizeof(double) * stateCount * stateCount * stateCount);

	eigenValues = (double *)malloc(sizeof(double) * stateCount);
	storedEigenValues = (double *)malloc(sizeof(double) * stateCount);

	frequencies = (double *)malloc(sizeof(double) * stateCount);
	storedFrequencies = (double *)malloc(sizeof(double) * stateCount);

	categoryRates = (double *)malloc(sizeof(double) * matrixCount);
	storedCategoryRates = (double *)malloc(sizeof(double) * matrixCount);

	categoryProportions = (double *)malloc(sizeof(double) * matrixCount);
	storedCategoryProportions = (double *)malloc(sizeof(double) * matrixCount);

	// a temporary array used in calculating log likelihoods
	integrationTmp = (double *)malloc(sizeof(double) * patternCount * STATE_COUNT);

	partialsSize = patternCount * stateCount * matrixCount;

	partials = (double ***)malloc(sizeof(double**) * 2);
	partials[0] = (double **)malloc(sizeof(double*) * nodeCount);
	partials[1] = (double **)malloc(sizeof(double*) * nodeCount);
	for (int i = 0; i < nodeCount; i++) {
		partials[0][i] = (double *)malloc(sizeof(double) * partialsSize);
		partials[1][i] = (double *)malloc(sizeof(double) * partialsSize);
	}

  	currentMatricesIndices = (int *)malloc(sizeof(int) * nodeCount);
  	storedMatricesIndices = (int *)malloc(sizeof(int) * nodeCount);

  	currentPartialsIndices = (int *)malloc(sizeof(int) * nodeCount);
  	storedPartialsIndices = (int *)malloc(sizeof(int) * nodeCount);
       
    matrixSize = stateCount * stateCount;

	matrices = (double ***)malloc(sizeof(double**) * 2);
	matrices[0] = (double **)malloc(sizeof(double*) * nodeCount);
	matrices[1] = (double **)malloc(sizeof(double*) * nodeCount);
	for (int i = 0; i < nodeCount; i++) {
		matrices[0][i] = (double *)malloc(sizeof(double) * matrixCount * matrixSize);
		matrices[1][i] = (double *)malloc(sizeof(double) * matrixCount * matrixSize);
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
 * Method:    updateRootFrequencies
 * Signature: ([D)V
 */
JNIEXPORT void JNICALL Java_dr_evomodel_newtreelikelihood_NativeLikelihoodCore_updateRootFrequencies
  (JNIEnv *env, jobject obj, jdoubleArray inFrequencies) 
{
	jdouble *tmpFrequencies = (jdouble*)(*env)->GetPrimitiveArrayCritical(env, inFrequencies, 0);
	
	memcpy(frequencies, tmpFrequencies, sizeof(double) * STATE_COUNT);
	
	(*env)->ReleasePrimitiveArrayCritical(env, inFrequencies, tmpFrequencies, JNI_ABORT);
}

/*
 * Class:     dr_evomodel_newtreelikelihood_NativeLikelihoodCore
 * Method:    updateEigenDecomposition
 * Signature: ([D[D)V
 */
JNIEXPORT void JNICALL Java_dr_evomodel_newtreelikelihood_NativeLikelihoodCore_updateEigenDecomposition
  (JNIEnv *env, jobject obj, jdoubleArray inCMatrix, jdoubleArray inEigenValues)
{
	jdouble *tmpCMatrix = (jdouble*)(*env)->GetPrimitiveArrayCritical(env, inCMatrix, 0);
	jdouble *tmpEigenValues = (jdouble*)(*env)->GetPrimitiveArrayCritical(env, inEigenValues, 0);
	
	memcpy(cMatrix, tmpCMatrix, sizeof(double) * STATE_COUNT * STATE_COUNT * STATE_COUNT);
	memcpy(eigenValues, tmpEigenValues, sizeof(double) * STATE_COUNT);
	
	(*env)->ReleasePrimitiveArrayCritical(env, inEigenValues, tmpEigenValues, JNI_ABORT);
	(*env)->ReleasePrimitiveArrayCritical(env, inCMatrix, tmpCMatrix, JNI_ABORT);
}

/*
 * Class:     dr_evomodel_newtreelikelihood_NativeLikelihoodCore
 * Method:    updateCategoryRates
 * Signature: ([D)V
 */
JNIEXPORT void JNICALL Java_dr_evomodel_newtreelikelihood_NativeLikelihoodCore_updateCategoryRates
  (JNIEnv *env, jobject obj, jdoubleArray inRates)
{
	jdouble *tmpRates = (jdouble*)(*env)->GetPrimitiveArrayCritical(env, inRates, 0);
	
	memcpy(categoryRates, tmpRates, sizeof(double) * matrixCount);
	
	(*env)->ReleasePrimitiveArrayCritical(env, inRates, tmpRates, JNI_ABORT);
}

/*
 * Class:     dr_evomodel_newtreelikelihood_NativeLikelihoodCore
 * Method:    updateCategoryProportions
 * Signature: ([D)V
 */
JNIEXPORT void JNICALL Java_dr_evomodel_newtreelikelihood_NativeLikelihoodCore_updateCategoryProportions
  (JNIEnv *env, jobject obj, jdoubleArray inProportions)
{
	jdouble *tmpProportions = (jdouble*)(*env)->GetPrimitiveArrayCritical(env, inProportions, 0);
	
	memcpy(categoryProportions, tmpProportions, sizeof(double) * matrixCount);
	
	(*env)->ReleasePrimitiveArrayCritical(env, inProportions, tmpProportions, JNI_ABORT);
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
	
	for (int x; x < count; x++) {
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


		double* matrices1 = matrices[currentMatricesIndices[nodeIndex1]][nodeIndex1];
		double* matrices2 = matrices[currentMatricesIndices[nodeIndex2]][nodeIndex2];

		double* partials1 = partials[currentPartialsIndices[nodeIndex1]][nodeIndex1];
		double* partials2 = partials[currentPartialsIndices[nodeIndex2]][nodeIndex2];

		currentPartialsIndices[nodeIndex3] = 1 - currentPartialsIndices[nodeIndex3];
		double* partials3 = partials[currentPartialsIndices[nodeIndex3]][nodeIndex3];

		double sum1, sum2;

		int u = 0;
		int v = 0;

		for (int l = 0; l < matrixCount; l++) {

			for (int k = 0; k < patternCount; k++) {

				int w = l * matrixSize;

				for (int i = 0; i < STATE_COUNT; i++) {

					sum1 = sum2 = 0.0;

					for (int j = 0; j < STATE_COUNT; j++) {
						sum1 += matrices1[w] * partials1[v + j];
						sum2 += matrices2[w] * partials2[v + j];

						w++;
					}

					partials3[u] = sum1 * sum2;
					u++;
				}
				v += STATE_COUNT;
			}
		}
	}
	/* not required:
	(*env)->ReleasePrimitiveArrayCritical(env, inDependencies, dependencies, JNI_ABORT);
	*/
	(*env)->ReleasePrimitiveArrayCritical(env, inOperations, operations, JNI_ABORT);

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

