
#include "NucleotideLikelihoodCore.h"

#define STATE_COUNT 4
#define MATRIX_SIZE STATE_COUNT * STATE_COUNT

/**
	 NucleotideLikelihoodCore - An implementation of LikelihoodCore for nucleotides 
	 that calls native methods for maximum speed. The native methods should be in
	 a shared library called "NucleotideLikelihoodCore" but the exact name will be system
	 dependent (i.e. "libNucleotideLikelihoodCore.so" or "NucleotideLikelihoodCore.dll").
	 @version $Id: NucleotideLikelihoodCore.c,v 1.4 2004/03/16 17:56:35 rambaut Exp $
	 @author Andrew Rambaut
*/


/*
 * Class:     dr_evomodel_treelikelihood_NativeNucleotideLikelihoodCore
 * Method:    nativeStatesStatesPruning
 * Signature: ([I[D[I[DII[D)V
 */

/**
	 Calculates partial likelihoods at a node when both children have states.
*/

JNIEXPORT void JNICALL Java_dr_evomodel_treelikelihood_NativeNucleotideLikelihoodCore_nativeStatesStatesPruning
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
							
			if (state1 < 4 && state2 < 4) {
			
				partials3[v] = matrices1[w + state1] * matrices2[w + state2];
				v++;	w += 4;
				partials3[v] = matrices1[w + state1] * matrices2[w + state2];
				v++;	w += 4;
				partials3[v] = matrices1[w + state1] * matrices2[w + state2];
				v++;	w += 4;
				partials3[v] = matrices1[w + state1] * matrices2[w + state2];
				v++;	w += 4;
				
			} else if (state1 < 4) {
				// child 2 has a gap or unknown state so don't use it

				partials3[v] = matrices1[w + state1];
				v++;	w += 4;
				partials3[v] = matrices1[w + state1];
				v++;	w += 4;
				partials3[v] = matrices1[w + state1];
				v++;	w += 4;
				partials3[v] = matrices1[w + state1];
				v++;	w += 4;
				
			} else if (state2 < 4) {
				// child 2 has a gap or unknown state so don't use it
				partials3[v] = matrices2[w + state2];
				v++;	w += 4;
				partials3[v] = matrices2[w + state2];
				v++;	w += 4;
				partials3[v] = matrices2[w + state2];
				v++;	w += 4;
				partials3[v] = matrices2[w + state2];
				v++;	w += 4;
				
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
 * Class:     dr_evomodel_treelikelihood_NativeNucleotideLikelihoodCore
 * Method:    nativeStatesPartialsPruning
 * Signature: ([I[D[D[DII[D)V
 */
JNIEXPORT void JNICALL Java_dr_evomodel_treelikelihood_NativeNucleotideLikelihoodCore_nativeStatesPartialsPruning
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
	double sum;
	int u = 0;
	int v = 0;
	int w = 0;

	for (l = 0; l < matrixCount; l++) {
		for (k = 0; k < patternCount; k++) {
		
			state1 = states1[k];
		
			if (state1 < 4) {
						
				sum =	matrices2[w] * partials2[v];
				sum +=	matrices2[w + 1] * partials2[v + 1];
				sum +=	matrices2[w + 2] * partials2[v + 2];
				sum +=	matrices2[w + 3] * partials2[v + 3];
				partials3[u] = matrices1[w + state1] * sum;	u++;
				
				sum =	matrices2[w + 4] * partials2[v];
				sum +=	matrices2[w + 5] * partials2[v + 1];
				sum +=	matrices2[w + 6] * partials2[v + 2];
				sum +=	matrices2[w + 7] * partials2[v + 3];
				partials3[u] = matrices1[w + 4 + state1] * sum;	u++;
				
				sum =	matrices2[w + 8] * partials2[v];
				sum +=	matrices2[w + 9] * partials2[v + 1];
				sum +=	matrices2[w + 10] * partials2[v + 2];
				sum +=	matrices2[w + 11] * partials2[v + 3];
				partials3[u] = matrices1[w + 8 + state1] * sum;	u++;
				
				sum =	matrices2[w + 12] * partials2[v];
				sum +=	matrices2[w + 13] * partials2[v + 1];
				sum +=	matrices2[w + 14] * partials2[v + 2];
				sum +=	matrices2[w + 15] * partials2[v + 3];
				partials3[u] = matrices1[w + 12 + state1] * sum;	u++;
				
				v += 4;

			} else {
				// Child 1 has a gap or unknown state so don't use it
				
				sum =	matrices2[w] * partials2[v];
				sum +=	matrices2[w + 1] * partials2[v + 1];
				sum +=	matrices2[w + 2] * partials2[v + 2];
				sum +=	matrices2[w + 3] * partials2[v + 3];
				partials3[u] = sum;	u++;
				
				sum =	matrices2[w + 4] * partials2[v];
				sum +=	matrices2[w + 5] * partials2[v + 1];
				sum +=	matrices2[w + 6] * partials2[v + 2];
				sum +=	matrices2[w + 7] * partials2[v + 3];
				partials3[u] = sum;	u++;
				
				sum =	matrices2[w + 8] * partials2[v];
				sum +=	matrices2[w + 9] * partials2[v + 1];
				sum +=	matrices2[w + 10] * partials2[v + 2];
				sum +=	matrices2[w + 11] * partials2[v + 3];
				partials3[u] = sum;	u++;
				
				sum =	matrices2[w + 12] * partials2[v];
				sum +=	matrices2[w + 13] * partials2[v + 1];
				sum +=	matrices2[w + 14] * partials2[v + 2];
				sum +=	matrices2[w + 15] * partials2[v + 3];
				partials3[u] = sum;	u++;
				
				v += 4;
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
 * Class:     dr_evomodel_treelikelihood_NativeNucleotideLikelihoodCore
 * Method:    nativePartialsPartialsPruning
 * Signature: ([D[D[D[DII[D)V
 */
JNIEXPORT void JNICALL Java_dr_evomodel_treelikelihood_NativeNucleotideLikelihoodCore_nativePartialsPartialsPruning
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

	int k, l;
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
			partials3[u] = sum1 * sum2; u++;

			sum1 = matrices1[w + 4] * partials1[v];
			sum2 = matrices2[w + 4] * partials2[v];	
			sum1 += matrices1[w + 5] * partials1[v + 1];
			sum2 += matrices2[w + 5] * partials2[v + 1];
			sum1 += matrices1[w + 6] * partials1[v + 2];
			sum2 += matrices2[w + 6] * partials2[v + 2];
			sum1 += matrices1[w + 7] * partials1[v + 3];
			sum2 += matrices2[w + 7] * partials2[v + 3];
			partials3[u] = sum1 * sum2; u++;

			sum1 = matrices1[w + 8] * partials1[v];
			sum2 = matrices2[w + 8] * partials2[v];	
			sum1 += matrices1[w + 9] * partials1[v + 1];
			sum2 += matrices2[w + 9] * partials2[v + 1];
			sum1 += matrices1[w + 10] * partials1[v + 2];
			sum2 += matrices2[w + 10] * partials2[v + 2];
			sum1 += matrices1[w + 11] * partials1[v + 3];
			sum2 += matrices2[w + 11] * partials2[v + 3];
			partials3[u] = sum1 * sum2; u++;

			sum1 = matrices1[w + 12] * partials1[v];
			sum2 = matrices2[w + 12] * partials2[v];	
			sum1 += matrices1[w + 13] * partials1[v + 1];
			sum2 += matrices2[w + 13] * partials2[v + 1];
			sum1 += matrices1[w + 14] * partials1[v + 2];
			sum2 += matrices2[w + 14] * partials2[v + 2];
			sum1 += matrices1[w + 15] * partials1[v + 3];
			sum2 += matrices2[w + 15] * partials2[v + 3];
			partials3[u] = sum1 * sum2; u++;

			v += 4;

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
 * Class:     dr_evomodel_treelikelihood_NativeNucleotideLikelihoodCore
 * Method:    nativeIntegratePartials
 * Signature: ([D[DII[D)V
 */
JNIEXPORT void JNICALL Java_dr_evomodel_treelikelihood_NativeNucleotideLikelihoodCore_nativeIntegratePartials
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
	}
	
	
	for (j = 1; j < matrixCount; j++) {
		u = 0;	
		for (k = 0; k < patternCount; k++) {
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
							
