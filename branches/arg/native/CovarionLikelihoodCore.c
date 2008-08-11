
#include "CovarionLikelihoodCore.h"

#define STATE_COUNT 8
#define MATRIX_SIZE STATE_COUNT * STATE_COUNT

/**
	 CovarionLikelihoodCore - An implementation of LikelihoodCore for Covarions 
	 that calls native methods for maximum speed. The native methods should be in
	 a shared library called "CovarionLikelihoodCore" but the exact name will be system
	 dependent (i.e. "libCovarionLikelihoodCore.so" or "CovarionLikelihoodCore.dll").
	 @version $Id: CovarionLikelihoodCore.c,v 1.4 2004/03/16 17:56:35 rambaut Exp $
	 @author Erik Bloomquist
	 
	 
	 June 16, 2008
	 
	 Note that the methods
	 1. Java_dr_evomodel_treelikelihood_NativeCovarionLikelihoodCore_nativeStatesPartialsPruning
	 2. Java_dr_evomodel_treelikelihood_NativeCovarionLikelihoodCore_nativePartialsPartialsPruning
	 
	 have been unrolled.  Seems to be slightly faster than the loop versions,
	 but only by a mild amount 3-5%.  The looped versions have been commented out
	 and are at the bottom of the code.
*/


/*
 * Class:     dr_evomodel_treelikelihood_NativeCovarionLikelihoodCore
 * Method:    nativeStatesStatesPruning
 * Signature: ([I[D[I[DII[D)V
 */

/**
	 Calculates partial likelihoods at a node when both children have states.
*/

JNIEXPORT void JNICALL Java_dr_evomodel_treelikelihood_NativeCovarionLikelihoodCore_nativeStatesStatesPruning
	(JNIEnv *env, jobject obj, jintArray inStates1, jdoubleArray inMatrices1, 
								jintArray inStates2, jdoubleArray inMatrices2, 
								jint patternCount, jint matrixCount,
								jdoubleArray outPartials)
{
	int j, k, u, v, w, state1, state2;

	jint *states1 = (jint*)(*env)->GetPrimitiveArrayCritical(env, inStates1, 0);
	jdouble *matrices1 = (jdouble*)(*env)->GetPrimitiveArrayCritical(env, inMatrices1, 0);
	jint *states2 = (jint*)(*env)->GetPrimitiveArrayCritical(env, inStates2, 0);
	jdouble *matrices2 = (jdouble*)(*env)->GetPrimitiveArrayCritical(env, inMatrices2, 0);
	jdouble *partials3 = (jdouble*)(*env)->GetPrimitiveArrayCritical(env, outPartials, 0);

	v = 0;
	u = 0;
	for (j = 0; j < matrixCount; j++) {

		for (k = 0; k < patternCount; k++) {
		
			w = u;
			
			state1 = states1[k];
			state2 = states2[k];
							
			if (state1 < 8 && state2 < 8) {
			
				partials3[v] = matrices1[w + state1] * matrices2[w + state2];
				v++;	w += 8;
				partials3[v] = matrices1[w + state1] * matrices2[w + state2];
				v++;	w += 8;
				partials3[v] = matrices1[w + state1] * matrices2[w + state2];
				v++;	w += 8;
				partials3[v] = matrices1[w + state1] * matrices2[w + state2];
				v++;	w += 8;
				partials3[v] = matrices1[w + state1] * matrices2[w + state2];
				v++;	w += 8;
				partials3[v] = matrices1[w + state1] * matrices2[w + state2];
				v++;	w += 8;
				partials3[v] = matrices1[w + state1] * matrices2[w + state2];
				v++;	w += 8;
				partials3[v] = matrices1[w + state1] * matrices2[w + state2];
				v++;	w += 8;
				
				
				
				
							
				
			} else if (state1 < 8) {
				// child 2 has a gap or unknown state so don't use it

				partials3[v] = matrices1[w + state1];
				v++;	w += 8;
				partials3[v] = matrices1[w + state1];
				v++;	w += 8;
				partials3[v] = matrices1[w + state1];
				v++;	w += 8;
				partials3[v] = matrices1[w + state1];
				v++;	w += 8;
				partials3[v] = matrices1[w + state1];
				v++;	w += 8;
				partials3[v] = matrices1[w + state1];
				v++;	w += 8;
				partials3[v] = matrices1[w + state1];
				v++;	w += 8;
				partials3[v] = matrices1[w + state1];
				v++;	w += 8;
				
				
				
				
				
			} else if (state2 < 8) {
				// child 2 has a gap or unknown state so don't use it
				partials3[v] = matrices2[w + state2];
				v++;	w += 8;
				partials3[v] = matrices2[w + state2];
				v++;	w += 8;
				partials3[v] = matrices2[w + state2];
				v++;	w += 8;
				partials3[v] = matrices2[w + state2];
				v++;	w += 8;
				partials3[v] = matrices2[w + state2];
				v++;	w += 8;
				partials3[v] = matrices2[w + state2];
				v++;	w += 8;
				partials3[v] = matrices2[w + state2];
				v++;	w += 8;
				partials3[v] = matrices2[w + state2];
				v++;	w += 8;
				

				
			} else {
				// both children have a gap or unknown state so set partials to 1
				partials3[v] = 1.0;
				v++;
				partials3[v] = 1.0;
				v++;
				partials3[v] = 1.0;
				v++;
				partials3[v] = 1.0;
				v++;
				partials3[v] = 1.0;
				v++;
				partials3[v] = 1.0;
				v++;
				partials3[v] = 1.0;
				v++;
				partials3[v] = 1.0;
				v++;
			}
		}
		
		u += MATRIX_SIZE;
	}
	
	
	(*env)->ReleasePrimitiveArrayCritical(env, outPartials, partials3, 0);
	(*env)->ReleasePrimitiveArrayCritical(env, inMatrices2, matrices2, JNI_ABORT);
	(*env)->ReleasePrimitiveArrayCritical(env, inStates2, states2, JNI_ABORT);
	(*env)->ReleasePrimitiveArrayCritical(env, inMatrices1, matrices1, JNI_ABORT);
	(*env)->ReleasePrimitiveArrayCritical(env, inStates1, states1, JNI_ABORT);
}


/*
 * Class:     dr_evomodel_treelikelihood_NativeCovarionLikelihoodCore
 * Method:    nativeStatesPartialsPruning
 * Signature: ([I[D[D[DII[D)V
 */
JNIEXPORT void JNICALL Java_dr_evomodel_treelikelihood_NativeCovarionLikelihoodCore_nativeStatesPartialsPruning
	(JNIEnv *env, jobject obj, jintArray inStates1, jdoubleArray inMatrices1, 
								jdoubleArray inPartials2, jdoubleArray inMatrices2, 
								jint patternCount, jint matrixCount,
								jdoubleArray outPartials)
{
	jint *states1 = (jint*)(*env)->GetPrimitiveArrayCritical(env, inStates1, 0);
	jdouble *matrices1 = (jdouble*)(*env)->GetPrimitiveArrayCritical(env, inMatrices1, 0);
	jdouble *partials2 = (jdouble*)(*env)->GetPrimitiveArrayCritical(env, inPartials2, 0);
	jdouble *matrices2 = (jdouble*)(*env)->GetPrimitiveArrayCritical(env, inMatrices2, 0);
	jdouble *partials3 = (jdouble*)(*env)->GetPrimitiveArrayCritical(env, outPartials, 0);
	
	int l, k, state1;
	int u = 0;
	int v = 0;
	int w = 0;
	double sum;
	

	for (l = 0; l < matrixCount; l++) {
		for (k = 0; k < patternCount; k++) {
			state1 = states1[k];

				if (state1 < 8) {

				 sum = 0;

sum = matrices2[w] * partials2[v];
sum += matrices2[w + 1] * partials2[v + 1];
sum += matrices2[w + 2] * partials2[v + 2];
sum += matrices2[w + 3] * partials2[v + 3];
sum += matrices2[w + 4] * partials2[v + 4];
sum += matrices2[w + 5] * partials2[v + 5];
sum += matrices2[w + 6] * partials2[v + 6];
sum += matrices2[w + 7] * partials2[v + 7];
partials3[u] = matrices1[w + state1] * sum;	u++;

sum += matrices2[w + 8] * partials2[v];
sum += matrices2[w + 9] * partials2[v + 1];
sum += matrices2[w + 10] * partials2[v + 2];
sum += matrices2[w + 11] * partials2[v + 3];
sum += matrices2[w + 12] * partials2[v + 4];
sum += matrices2[w + 13] * partials2[v + 5];
sum += matrices2[w + 14] * partials2[v + 6];
sum += matrices2[w + 15] * partials2[v + 7];
partials3[u] = matrices1[w + 8 + state1] * sum;	u++;

sum += matrices2[w + 16] * partials2[v];
sum += matrices2[w + 17] * partials2[v + 1];
sum += matrices2[w + 18] * partials2[v + 2];
sum += matrices2[w + 19] * partials2[v + 3];
sum += matrices2[w + 20] * partials2[v + 4];
sum += matrices2[w + 21] * partials2[v + 5];
sum += matrices2[w + 22] * partials2[v + 6];
sum += matrices2[w + 23] * partials2[v + 7];
partials3[u] = matrices1[w + 16 + state1] * sum;	u++;

sum += matrices2[w + 24] * partials2[v];
sum += matrices2[w + 25] * partials2[v + 1];
sum += matrices2[w + 26] * partials2[v + 2];
sum += matrices2[w + 27] * partials2[v + 3];
sum += matrices2[w + 28] * partials2[v + 4];
sum += matrices2[w + 29] * partials2[v + 5];
sum += matrices2[w + 30] * partials2[v + 6];
sum += matrices2[w + 31] * partials2[v + 7];
partials3[u] = matrices1[w + 24 + state1] * sum;	u++;

sum += matrices2[w + 32] * partials2[v];
sum += matrices2[w + 33] * partials2[v + 1];
sum += matrices2[w + 34] * partials2[v + 2];
sum += matrices2[w + 35] * partials2[v + 3];
sum += matrices2[w + 36] * partials2[v + 4];
sum += matrices2[w + 37] * partials2[v + 5];
sum += matrices2[w + 38] * partials2[v + 6];
sum += matrices2[w + 39] * partials2[v + 7];
partials3[u] = matrices1[w + 32 + state1] * sum;	u++;

sum += matrices2[w + 40] * partials2[v];
sum += matrices2[w + 41] * partials2[v + 1];
sum += matrices2[w + 42] * partials2[v + 2];
sum += matrices2[w + 43] * partials2[v + 3];
sum += matrices2[w + 44] * partials2[v + 4];
sum += matrices2[w + 45] * partials2[v + 5];
sum += matrices2[w + 46] * partials2[v + 6];
sum += matrices2[w + 47] * partials2[v + 7];
partials3[u] = matrices1[w + 40 + state1] * sum;	u++;

sum += matrices2[w + 48] * partials2[v];
sum += matrices2[w + 49] * partials2[v + 1];
sum += matrices2[w + 50] * partials2[v + 2];
sum += matrices2[w + 51] * partials2[v + 3];
sum += matrices2[w + 52] * partials2[v + 4];
sum += matrices2[w + 53] * partials2[v + 5];
sum += matrices2[w + 54] * partials2[v + 6];
sum += matrices2[w + 55] * partials2[v + 7];
partials3[u] = matrices1[w + 48 + state1] * sum;	u++;

sum += matrices2[w + 56] * partials2[v];
sum += matrices2[w + 57] * partials2[v + 1];
sum += matrices2[w + 58] * partials2[v + 2];
sum += matrices2[w + 59] * partials2[v + 3];
sum += matrices2[w + 60] * partials2[v + 4];
sum += matrices2[w + 61] * partials2[v + 5];
sum += matrices2[w + 62] * partials2[v + 6];
sum += matrices2[w + 63] * partials2[v + 7];
partials3[u] = matrices1[w + 56 + state1] * sum;	u++;



					v += 8;

				} else {
					// Child 1 has a gap or unknown state so don't use it

				sum = 0;

sum = matrices2[w] * partials2[v];
sum += matrices2[w + 1] * partials2[v + 1];
sum += matrices2[w + 2] * partials2[v + 2];
sum += matrices2[w + 3] * partials2[v + 3];
sum += matrices2[w + 4] * partials2[v + 4];
sum += matrices2[w + 5] * partials2[v + 5];
sum += matrices2[w + 6] * partials2[v + 6];
sum += matrices2[w + 7] * partials2[v + 7];
partials3[u] =  sum;	u++;

sum += matrices2[w + 8] * partials2[v];
sum += matrices2[w + 9] * partials2[v + 1];
sum += matrices2[w + 10] * partials2[v + 2];
sum += matrices2[w + 11] * partials2[v + 3];
sum += matrices2[w + 12] * partials2[v + 4];
sum += matrices2[w + 13] * partials2[v + 5];
sum += matrices2[w + 14] * partials2[v + 6];
sum += matrices2[w + 15] * partials2[v + 7];
partials3[u] =  sum;	u++;

sum += matrices2[w + 16] * partials2[v];
sum += matrices2[w + 17] * partials2[v + 1];
sum += matrices2[w + 18] * partials2[v + 2];
sum += matrices2[w + 19] * partials2[v + 3];
sum += matrices2[w + 20] * partials2[v + 4];
sum += matrices2[w + 21] * partials2[v + 5];
sum += matrices2[w + 22] * partials2[v + 6];
sum += matrices2[w + 23] * partials2[v + 7];
partials3[u] =  sum;	u++;

sum += matrices2[w + 24] * partials2[v];
sum += matrices2[w + 25] * partials2[v + 1];
sum += matrices2[w + 26] * partials2[v + 2];
sum += matrices2[w + 27] * partials2[v + 3];
sum += matrices2[w + 28] * partials2[v + 4];
sum += matrices2[w + 29] * partials2[v + 5];
sum += matrices2[w + 30] * partials2[v + 6];
sum += matrices2[w + 31] * partials2[v + 7];
partials3[u] =  sum;	u++;

sum += matrices2[w + 32] * partials2[v];
sum += matrices2[w + 33] * partials2[v + 1];
sum += matrices2[w + 34] * partials2[v + 2];
sum += matrices2[w + 35] * partials2[v + 3];
sum += matrices2[w + 36] * partials2[v + 4];
sum += matrices2[w + 37] * partials2[v + 5];
sum += matrices2[w + 38] * partials2[v + 6];
sum += matrices2[w + 39] * partials2[v + 7];
partials3[u] =  sum;	u++;

sum += matrices2[w + 40] * partials2[v];
sum += matrices2[w + 41] * partials2[v + 1];
sum += matrices2[w + 42] * partials2[v + 2];
sum += matrices2[w + 43] * partials2[v + 3];
sum += matrices2[w + 44] * partials2[v + 4];
sum += matrices2[w + 45] * partials2[v + 5];
sum += matrices2[w + 46] * partials2[v + 6];
sum += matrices2[w + 47] * partials2[v + 7];
partials3[u] =  sum;	u++;

sum += matrices2[w + 48] * partials2[v];
sum += matrices2[w + 49] * partials2[v + 1];
sum += matrices2[w + 50] * partials2[v + 2];
sum += matrices2[w + 51] * partials2[v + 3];
sum += matrices2[w + 52] * partials2[v + 4];
sum += matrices2[w + 53] * partials2[v + 5];
sum += matrices2[w + 54] * partials2[v + 6];
sum += matrices2[w + 55] * partials2[v + 7];
partials3[u] =  sum;	u++;

sum += matrices2[w + 56] * partials2[v];
sum += matrices2[w + 57] * partials2[v + 1];
sum += matrices2[w + 58] * partials2[v + 2];
sum += matrices2[w + 59] * partials2[v + 3];
sum += matrices2[w + 60] * partials2[v + 4];
sum += matrices2[w + 61] * partials2[v + 5];
sum += matrices2[w + 62] * partials2[v + 6];
sum += matrices2[w + 63] * partials2[v + 7];
partials3[u] =  sum;	u++;



					v += 8;
				}
			}

			w += MATRIX_SIZE;
		}
	(*env)->ReleasePrimitiveArrayCritical(env, outPartials, partials3, 0);
	(*env)->ReleasePrimitiveArrayCritical(env, inMatrices2, matrices2, JNI_ABORT);
	(*env)->ReleasePrimitiveArrayCritical(env, inPartials2, partials2, JNI_ABORT);
	(*env)->ReleasePrimitiveArrayCritical(env, inMatrices1, matrices1, JNI_ABORT);
	(*env)->ReleasePrimitiveArrayCritical(env, inStates1, states1, JNI_ABORT);
}

/*
 * Class:     dr_evomodel_treelikelihood_NativeCovarionLikelihoodCore
 * Method:    nativePartialsPartialsPruning
 * Signature: ([D[D[D[DII[D)V
 */
JNIEXPORT void JNICALL Java_dr_evomodel_treelikelihood_NativeCovarionLikelihoodCore_nativePartialsPartialsPruning
	(JNIEnv *env, jobject obj, jdoubleArray inPartials1, jdoubleArray inMatrices1, 
								jdoubleArray inPartials2, jdoubleArray inMatrices2, 
								jint patternCount, jint matrixCount,
								jdoubleArray outPartials)
{
	jdouble *partials1 = (jdouble*)(*env)->GetPrimitiveArrayCritical(env, inPartials1, 0);
	jdouble *matrices1 = (jdouble*)(*env)->GetPrimitiveArrayCritical(env, inMatrices1, 0);
	jdouble *partials2 = (jdouble*)(*env)->GetPrimitiveArrayCritical(env, inPartials2, 0);
	jdouble *matrices2 = (jdouble*)(*env)->GetPrimitiveArrayCritical(env, inMatrices2, 0);
	jdouble *partials3 = (jdouble*)(*env)->GetPrimitiveArrayCritical(env, outPartials, 0);

	int i, k, l;
	double sum1, sum2;
	
	int u = 0;
	int v = 0;
	int w = 0;

	
	for (l = 0; l < matrixCount; l++) {
			for (k = 0; k < patternCount; k++) {

sum1 = matrices1[w] * partials1[v];
sum2 = matrices2[w] * partials2[v];
sum1 += matrices1[w + 1] * partials1[v + 1];
sum2 += matrices2[w + 1] * partials2[v + 1];
sum1 += matrices1[w + 2] * partials1[v + 2];
sum2 += matrices2[w + 2] * partials2[v + 2];
sum1 += matrices1[w + 3] * partials1[v + 3];
sum2 += matrices2[w + 3] * partials2[v + 3];
sum1 += matrices1[w + 4] * partials1[v + 4];
sum2 += matrices2[w + 4] * partials2[v + 4];
sum1 += matrices1[w + 5] * partials1[v + 5];
sum2 += matrices2[w + 5] * partials2[v + 5];
sum1 += matrices1[w + 6] * partials1[v + 6];
sum2 += matrices2[w + 6] * partials2[v + 6];
sum1 += matrices1[w + 7] * partials1[v + 7];
sum2 += matrices2[w + 7] * partials2[v + 7];
partials3[u] = sum1 * sum2; u++;

sum1 = matrices1[w + 8] * partials1[v];
sum2 = matrices2[w + 8] * partials2[v];
sum1 += matrices1[w + 9] * partials1[v + 1];
sum2 += matrices2[w + 9] * partials2[v + 1];
sum1 += matrices1[w + 10] * partials1[v + 2];
sum2 += matrices2[w + 10] * partials2[v + 2];
sum1 += matrices1[w + 11] * partials1[v + 3];
sum2 += matrices2[w + 11] * partials2[v + 3];
sum1 += matrices1[w + 12] * partials1[v + 4];
sum2 += matrices2[w + 12] * partials2[v + 4];
sum1 += matrices1[w + 13] * partials1[v + 5];
sum2 += matrices2[w + 13] * partials2[v + 5];
sum1 += matrices1[w + 14] * partials1[v + 6];
sum2 += matrices2[w + 14] * partials2[v + 6];
sum1 += matrices1[w + 15] * partials1[v + 7];
sum2 += matrices2[w + 15] * partials2[v + 7];
partials3[u] = sum1 * sum2; u++;

sum1 = matrices1[w + 16] * partials1[v];
sum2 = matrices2[w + 16] * partials2[v];
sum1 += matrices1[w + 17] * partials1[v + 1];
sum2 += matrices2[w + 17] * partials2[v + 1];
sum1 += matrices1[w + 18] * partials1[v + 2];
sum2 += matrices2[w + 18] * partials2[v + 2];
sum1 += matrices1[w + 19] * partials1[v + 3];
sum2 += matrices2[w + 19] * partials2[v + 3];
sum1 += matrices1[w + 20] * partials1[v + 4];
sum2 += matrices2[w + 20] * partials2[v + 4];
sum1 += matrices1[w + 21] * partials1[v + 5];
sum2 += matrices2[w + 21] * partials2[v + 5];
sum1 += matrices1[w + 22] * partials1[v + 6];
sum2 += matrices2[w + 22] * partials2[v + 6];
sum1 += matrices1[w + 23] * partials1[v + 7];
sum2 += matrices2[w + 23] * partials2[v + 7];
partials3[u] = sum1 * sum2; u++;

sum1 = matrices1[w + 24] * partials1[v];
sum2 = matrices2[w + 24] * partials2[v];
sum1 += matrices1[w + 25] * partials1[v + 1];
sum2 += matrices2[w + 25] * partials2[v + 1];
sum1 += matrices1[w + 26] * partials1[v + 2];
sum2 += matrices2[w + 26] * partials2[v + 2];
sum1 += matrices1[w + 27] * partials1[v + 3];
sum2 += matrices2[w + 27] * partials2[v + 3];
sum1 += matrices1[w + 28] * partials1[v + 4];
sum2 += matrices2[w + 28] * partials2[v + 4];
sum1 += matrices1[w + 29] * partials1[v + 5];
sum2 += matrices2[w + 29] * partials2[v + 5];
sum1 += matrices1[w + 30] * partials1[v + 6];
sum2 += matrices2[w + 30] * partials2[v + 6];
sum1 += matrices1[w + 31] * partials1[v + 7];
sum2 += matrices2[w + 31] * partials2[v + 7];
partials3[u] = sum1 * sum2; u++;

sum1 = matrices1[w + 32] * partials1[v];
sum2 = matrices2[w + 32] * partials2[v];
sum1 += matrices1[w + 33] * partials1[v + 1];
sum2 += matrices2[w + 33] * partials2[v + 1];
sum1 += matrices1[w + 34] * partials1[v + 2];
sum2 += matrices2[w + 34] * partials2[v + 2];
sum1 += matrices1[w + 35] * partials1[v + 3];
sum2 += matrices2[w + 35] * partials2[v + 3];
sum1 += matrices1[w + 36] * partials1[v + 4];
sum2 += matrices2[w + 36] * partials2[v + 4];
sum1 += matrices1[w + 37] * partials1[v + 5];
sum2 += matrices2[w + 37] * partials2[v + 5];
sum1 += matrices1[w + 38] * partials1[v + 6];
sum2 += matrices2[w + 38] * partials2[v + 6];
sum1 += matrices1[w + 39] * partials1[v + 7];
sum2 += matrices2[w + 39] * partials2[v + 7];
partials3[u] = sum1 * sum2; u++;

sum1 = matrices1[w + 40] * partials1[v];
sum2 = matrices2[w + 40] * partials2[v];
sum1 += matrices1[w + 41] * partials1[v + 1];
sum2 += matrices2[w + 41] * partials2[v + 1];
sum1 += matrices1[w + 42] * partials1[v + 2];
sum2 += matrices2[w + 42] * partials2[v + 2];
sum1 += matrices1[w + 43] * partials1[v + 3];
sum2 += matrices2[w + 43] * partials2[v + 3];
sum1 += matrices1[w + 44] * partials1[v + 4];
sum2 += matrices2[w + 44] * partials2[v + 4];
sum1 += matrices1[w + 45] * partials1[v + 5];
sum2 += matrices2[w + 45] * partials2[v + 5];
sum1 += matrices1[w + 46] * partials1[v + 6];
sum2 += matrices2[w + 46] * partials2[v + 6];
sum1 += matrices1[w + 47] * partials1[v + 7];
sum2 += matrices2[w + 47] * partials2[v + 7];
partials3[u] = sum1 * sum2; u++;

sum1 = matrices1[w + 48] * partials1[v];
sum2 = matrices2[w + 48] * partials2[v];
sum1 += matrices1[w + 49] * partials1[v + 1];
sum2 += matrices2[w + 49] * partials2[v + 1];
sum1 += matrices1[w + 50] * partials1[v + 2];
sum2 += matrices2[w + 50] * partials2[v + 2];
sum1 += matrices1[w + 51] * partials1[v + 3];
sum2 += matrices2[w + 51] * partials2[v + 3];
sum1 += matrices1[w + 52] * partials1[v + 4];
sum2 += matrices2[w + 52] * partials2[v + 4];
sum1 += matrices1[w + 53] * partials1[v + 5];
sum2 += matrices2[w + 53] * partials2[v + 5];
sum1 += matrices1[w + 54] * partials1[v + 6];
sum2 += matrices2[w + 54] * partials2[v + 6];
sum1 += matrices1[w + 55] * partials1[v + 7];
sum2 += matrices2[w + 55] * partials2[v + 7];
partials3[u] = sum1 * sum2; u++;

sum1 = matrices1[w + 56] * partials1[v];
sum2 = matrices2[w + 56] * partials2[v];
sum1 += matrices1[w + 57] * partials1[v + 1];
sum2 += matrices2[w + 57] * partials2[v + 1];
sum1 += matrices1[w + 58] * partials1[v + 2];
sum2 += matrices2[w + 58] * partials2[v + 2];
sum1 += matrices1[w + 59] * partials1[v + 3];
sum2 += matrices2[w + 59] * partials2[v + 3];
sum1 += matrices1[w + 60] * partials1[v + 4];
sum2 += matrices2[w + 60] * partials2[v + 4];
sum1 += matrices1[w + 61] * partials1[v + 5];
sum2 += matrices2[w + 61] * partials2[v + 5];
sum1 += matrices1[w + 62] * partials1[v + 6];
sum2 += matrices2[w + 62] * partials2[v + 6];
sum1 += matrices1[w + 63] * partials1[v + 7];
sum2 += matrices2[w + 63] * partials2[v + 7];
partials3[u] = sum1 * sum2; u++;





				v += 8;
			}

			w += MATRIX_SIZE;
		}



	(*env)->ReleasePrimitiveArrayCritical(env, outPartials, partials3, 0);
	(*env)->ReleasePrimitiveArrayCritical(env, inMatrices2, matrices2, JNI_ABORT);
	(*env)->ReleasePrimitiveArrayCritical(env, inPartials2, partials2, JNI_ABORT);
	(*env)->ReleasePrimitiveArrayCritical(env, inMatrices1, matrices1, JNI_ABORT);
	(*env)->ReleasePrimitiveArrayCritical(env, inPartials1, partials1, JNI_ABORT);
}




/*
 * Class:     dr_evomodel_treelikelihood_NativeCovarionLikelihoodCore
 * Method:    nativeIntegratePartials
 * Signature: ([D[DII[D)V
 */
JNIEXPORT void JNICALL Java_dr_evomodel_treelikelihood_NativeCovarionLikelihoodCore_nativeIntegratePartials
	(JNIEnv *env, jobject obj, jdoubleArray inPartials, jdoubleArray inProportions, 
								jint patternCount, jint matrixCount,
								jdoubleArray outPartials)
{
	jdouble *partials1 = (jdouble*)(*env)->GetPrimitiveArrayCritical(env, inPartials, 0);
	jdouble *proportions = (jdouble*)(*env)->GetPrimitiveArrayCritical(env, inProportions, 0);
	jdouble *partials3 = (jdouble*)(*env)->GetPrimitiveArrayCritical(env, outPartials, 0);

	int j, k;
		
	int u = 0;	
	int v = 0;
	
	for (k = 0; k < patternCount; k++) {
	
		partials3[u] = partials1[v] * proportions[0]; u++; v++;
		partials3[u] = partials1[v] * proportions[0]; u++; v++;
		partials3[u] = partials1[v] * proportions[0]; u++; v++;
		partials3[u] = partials1[v] * proportions[0]; u++; v++;

		partials3[u] = partials1[v] * proportions[0]; u++; v++;
		partials3[u] = partials1[v] * proportions[0]; u++; v++;
		partials3[u] = partials1[v] * proportions[0]; u++; v++;
		partials3[u] = partials1[v] * proportions[0]; u++; v++;
				
	}
	
	
	for (j = 1; j < matrixCount; j++) {
		u = 0;	
		for (k = 0; k < patternCount; k++) {
			partials3[u] += partials1[v] * proportions[j]; u++; v++;
			partials3[u] += partials1[v] * proportions[j]; u++; v++;
			partials3[u] += partials1[v] * proportions[j]; u++; v++;
			partials3[u] += partials1[v] * proportions[j]; u++; v++;

			partials3[u] += partials1[v] * proportions[j]; u++; v++;
			partials3[u] += partials1[v] * proportions[j]; u++; v++;
			partials3[u] += partials1[v] * proportions[j]; u++; v++;
			partials3[u] += partials1[v] * proportions[j]; u++; v++;

			


		}
	}

	(*env)->ReleasePrimitiveArrayCritical(env, outPartials, partials3, 0);
	(*env)->ReleasePrimitiveArrayCritical(env, inProportions, proportions, JNI_ABORT);
	(*env)->ReleasePrimitiveArrayCritical(env, inPartials, partials1, JNI_ABORT);
}


/*
 * Class:     dr_evomodel_treelikelihood_NativeCovarionLikelihoodCore
 * Method:    nativeIntegratePartials
 * Signature: ([D[DII[D)V
 */
 /*
JNIEXPORT void JNICALL Java_dr_evomodel_treelikelihood_NativeCovarionLikelihoodCore_nativeIntegratePartials
	(JNIEnv *env, jobject obj, jdoubleArray inPartials, jdoubleArray inProportions, 
								jint patternCount, jint matrixCount,
								jdoubleArray outPartials)
{
	jdouble *partials1 = (jdouble*)(*env)->GetPrimitiveArrayCritical(env, inPartials, 0);
	jdouble *proportions = (jdouble*)(*env)->GetPrimitiveArrayCritical(env, inProportions, 0);
	jdouble *partials3 = (jdouble*)(*env)->GetPrimitiveArrayCritical(env, outPartials, 0);

	int j, k;
		
	int u = 0;	
	int v = 0;
	
	for (k = 0; k < patternCount; k++) {
		partials3[u] = partials1[v] * proportions[0]; u++; v++;
		partials3[u] = partials1[v] * proportions[0]; u++; v++;
		partials3[u] = partials1[v] * proportions[0]; u++; v++;
		partials3[u] = partials1[v] * proportions[0]; u++; v++;
		partials3[u] = partials1[v] * proportions[0]; u++; v++;
		partials3[u] = partials1[v] * proportions[0]; u++; v++;
		partials3[u] = partials1[v] * proportions[0]; u++; v++;
		partials3[u] = partials1[v] * proportions[0]; u++; v++;
		
		
	}
	
	
	for (j = 1; j < matrixCount; j++) {
		u = 0;	
		for (k = 0; k < patternCount; k++) {
			partials3[u] += partials1[v] * proportions[j]; u++; v++;
			partials3[u] += partials1[v] * proportions[j]; u++; v++;
			partials3[u] += partials1[v] * proportions[j]; u++; v++;
			partials3[u] += partials1[v] * proportions[j]; u++; v++;
			partials3[u] += partials1[v] * proportions[j]; u++; v++;
			partials3[u] += partials1[v] * proportions[j]; u++; v++;
			partials3[u] += partials1[v] * proportions[j]; u++; v++;
			partials3[u] += partials1[v] * proportions[j]; u++; v++;
	
		}
	}

	(*env)->ReleasePrimitiveArrayCritical(env, outPartials, partials3, 0);
	(*env)->ReleasePrimitiveArrayCritical(env, inProportions, proportions, JNI_ABORT);
	(*env)->ReleasePrimitiveArrayCritical(env, inPartials, partials1, JNI_ABORT);
}
				*/
							

