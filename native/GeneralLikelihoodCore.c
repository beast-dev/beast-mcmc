#include <stdlib.h>
#include "GeneralLikelihoodCore.h"

/**

	 @author Marc Suchard
*/


/*
 * Class:     dr_evomodel_treelikelihood_NativeGeneralLikelihoodCore
 * Method:    nativeStatesStatesPruning
 * Signature: ([I[D[I[DII[D)V
 */

/**
	 Calculates partial likelihoods at a node when both children have states.
*/

JNIEXPORT void JNICALL Java_dr_evomodel_treelikelihood_NativeGeneralLikelihoodCore_nativeStatesStatesPruning
	(JNIEnv *env, jobject obj, jintArray inStates1, jdoubleArray inMatrices1, 
								jintArray inStates2, jdoubleArray inMatrices2, 
								jint patternCount, jint matrixCount,
								jdoubleArray outPartials,
								jint stateCount)
{
	int j, k, u, v, w, state1, state2, m;

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
							
			if (state1 < stateCount && state2 < stateCount) {
			
				for (m = 0; m < stateCount; m++) {
					partials3[v] = matrices1[w + state1] * matrices2[w + state2];
					v++;	w += stateCount;	
				}

			} else if (state1 < stateCount) {
			
				for (m = 0; m < stateCount; m++) {
					partials3[v] = matrices1[w + state1];
					v++;	w += stateCount;	
				}
							
			} else if (state2 < stateCount) {
			
				for (m = 0; m < stateCount; m++) {
					partials3[v] =  matrices2[w + state2];
					v++;	w += stateCount;	
				}
				
			} else {
				// both children have a gap or unknown state so set partials to 1
				
				for (m = 0; m < stateCount; m++) {
					partials3[v] = 1.0;
					v++;	
				}				
			}
		}
		
		u += stateCount * stateCount;
	}
	
	(*env)->ReleasePrimitiveArrayCritical(env, outPartials, partials3, 0);
	(*env)->ReleasePrimitiveArrayCritical(env, inMatrices2, matrices2, JNI_ABORT);
	(*env)->ReleasePrimitiveArrayCritical(env, inStates2, states2, JNI_ABORT);
	(*env)->ReleasePrimitiveArrayCritical(env, inMatrices1, matrices1, JNI_ABORT);
	(*env)->ReleasePrimitiveArrayCritical(env, inStates1, states1, JNI_ABORT);
}


/*
 * Class:     dr_evomodel_treelikelihood_NativeGeneralLikelihoodCore
 * Method:    nativeIntegratePartials
 * Signature: ([D[DII[D)V
 */
JNIEXPORT void JNICALL Java_dr_evomodel_treelikelihood_NativeGeneralLikelihoodCore_nativeIntegratePartials
	(JNIEnv *env, jobject obj, jdoubleArray inPartials, jdoubleArray inProportions, 
								jint patternCount, jint matrixCount,
								jdoubleArray outPartials,
								jint stateCount)
{
	jdouble *partials1 = (jdouble*)(*env)->GetPrimitiveArrayCritical(env, inPartials, 0);
	jdouble *proportions = (jdouble*)(*env)->GetPrimitiveArrayCritical(env, inProportions, 0);
	jdouble *partials3 = (jdouble*)(*env)->GetPrimitiveArrayCritical(env, outPartials, 0);

	int j, k, m;
		
	int u = 0;	
	int v = 0;
	
	for (k = 0; k < patternCount; k++) {
	
		for (m = 0; m < stateCount; m++) {
			partials3[u] = partials1[v] * proportions[0]; u++; v++;
		}
	}
	
	
	for (j = 1; j < matrixCount; j++) {
		u = 0;	
		for (k = 0; k < patternCount; k++) {
		
			for (m = 0; m < stateCount; m++) {
					partials3[u] += partials1[v] * proportions[j]; u++; v++;
			}			
		}
	}

	(*env)->ReleasePrimitiveArrayCritical(env, outPartials, partials3, 0);
	(*env)->ReleasePrimitiveArrayCritical(env, inProportions, proportions, JNI_ABORT);
	(*env)->ReleasePrimitiveArrayCritical(env, inPartials, partials1, JNI_ABORT);
}

/*
 * Class:     dr_evomodel_treelikelihood_NativeGeneralLikelihoodCore
 * Method:    nativeStatesPartialsPruning
 * Signature: ([I[D[D[DII[D)V
 */

JNIEXPORT void JNICALL Java_dr_evomodel_treelikelihood_NativeGeneralLikelihoodCore_nativeStatesPartialsPruning
	(JNIEnv *env, jobject obj, jintArray inStates1, jdoubleArray inMatrices1, 
								jdoubleArray inPartials2, jdoubleArray inMatrices2, 
								jint patternCount, jint matrixCount,
								jdoubleArray outPartials,
								jint stateCount)
{
	jint *states1 = (jint*)(*env)->GetPrimitiveArrayCritical(env, inStates1, 0);
	jdouble *matrices1 = (jdouble*)(*env)->GetPrimitiveArrayCritical(env, inMatrices1, 0);
	jdouble *partials2 = (jdouble*)(*env)->GetPrimitiveArrayCritical(env, inPartials2, 0);
	jdouble *matrices2 = (jdouble*)(*env)->GetPrimitiveArrayCritical(env, inMatrices2, 0);
	jdouble *partials3 = (jdouble*)(*env)->GetPrimitiveArrayCritical(env, outPartials, 0);
	
	int i, l, k, state1, m;
	double sum;
	int u = 0;
	int v = 0;
	int w = 0;
	int x, y;

	for (l = 0; l < matrixCount; l++) {
		for (k = 0; k < patternCount; k++) {
		
			int state1 = states1[k];
			
				if (state1 < stateCount) {
			
					x = w;
					for (i = 0; i < stateCount; i++) {
					
						y = v;
						double value = matrices1[x + state1];	
						
						sum = 0;
						for (m = 0; m < stateCount; m++) {
							sum += matrices2[x] * partials2[y]; x++; y++;
						}
																					
						partials3[u] = value * sum;	
						u++;
					}
											
					v += stateCount;

				} else {
	
					x = w;
					for (i = 0; i < stateCount; i++) {
					
						y = v;
						sum = 0;
						for (m = 0; m < stateCount; m++) {
							sum += matrices2[x] * partials2[y]; x++; y++;
						}
																	
						partials3[u] = sum;	
						u++;
					}
											
					v += stateCount;
				}
			}
			
			w += stateCount * stateCount;
		}
	
	(*env)->ReleasePrimitiveArrayCritical(env, outPartials, partials3, 0);
	(*env)->ReleasePrimitiveArrayCritical(env, inMatrices2, matrices2, JNI_ABORT);
	(*env)->ReleasePrimitiveArrayCritical(env, inPartials2, partials2, JNI_ABORT);
	(*env)->ReleasePrimitiveArrayCritical(env, inMatrices1, matrices1, JNI_ABORT);
	(*env)->ReleasePrimitiveArrayCritical(env, inStates1, states1, JNI_ABORT);
}
							

/*
 * Class:     dr_evomodel_treelikelihood_NativeGeneralLikelihoodCore
 * Method:    nativePartialsPartialsPruning
 * Signature: ([D[D[D[DII[D)V
 */

JNIEXPORT void JNICALL Java_dr_evomodel_treelikelihood_NativeGeneralLikelihoodCore_nativePartialsPartialsPruning
	(JNIEnv *env, jobject obj, jdoubleArray inPartials1, jdoubleArray inMatrices1, 
								jdoubleArray inPartials2, jdoubleArray inMatrices2, 
								jint patternCount, jint matrixCount,
								jdoubleArray outPartials,
								jint stateCount)
{
	jdouble *partials1 = (jdouble*)(*env)->GetPrimitiveArrayCritical(env, inPartials1, 0);
	jdouble *matrices1 = (jdouble*)(*env)->GetPrimitiveArrayCritical(env, inMatrices1, 0);
	jdouble *partials2 = (jdouble*)(*env)->GetPrimitiveArrayCritical(env, inPartials2, 0);
	jdouble *matrices2 = (jdouble*)(*env)->GetPrimitiveArrayCritical(env, inMatrices2, 0);
	jdouble *partials3 = (jdouble*)(*env)->GetPrimitiveArrayCritical(env, outPartials, 0);

	int i, k, l, j;
	double sum1, sum2;
	
	int u = 0;
	int v = 0;
	int w = 0;
	int x,y;

	for (l = 0; l < matrixCount; l++) {
		for (k = 0; k < patternCount; k++) {
		
				x = w;
				for (i = 0; i < stateCount; i++) {
					y = v;
					sum1 = 0;
					sum2 = 0;
					for (j =0; j < stateCount; j++) {
						sum1 += matrices1[x] * partials1[y];
						sum2 += matrices2[x] * partials2[y]; x++; y++;
					}										
					partials3[u] = sum1 * sum2;	
					u++;
				}
										
				v += stateCount;
		}
			
		w += stateCount * stateCount;
	}


	(*env)->ReleasePrimitiveArrayCritical(env, outPartials, partials3, 0);
	(*env)->ReleasePrimitiveArrayCritical(env, inMatrices2, matrices2, JNI_ABORT);
	(*env)->ReleasePrimitiveArrayCritical(env, inPartials2, partials2, JNI_ABORT);
	(*env)->ReleasePrimitiveArrayCritical(env, inMatrices1, matrices1, JNI_ABORT);
	(*env)->ReleasePrimitiveArrayCritical(env, inPartials1, partials1, JNI_ABORT);
}
