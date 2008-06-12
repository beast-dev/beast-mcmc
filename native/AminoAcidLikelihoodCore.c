
#include "AminoAcidLikelihoodCore.h"

#define STATE_COUNT 20
#define MATRIX_SIZE STATE_COUNT * STATE_COUNT

/**
	 AminoAcidLikelihoodCore - An implementation of LikelihoodCore for AminoAcids 
	 that calls native methods for maximum speed. The native methods should be in
	 a shared library called "AminoAcidLikelihoodCore" but the exact name will be system
	 dependent (i.e. "libAminoAcidLikelihoodCore.so" or "AminoAcidLikelihoodCore.dll").
	 @version $Id: AminoAcidLikelihoodCore.c,v 1.4 2004/03/16 17:56:35 rambaut Exp $
	 @author Andrew Rambaut
*/


/*
 * Class:     dr_evomodel_treelikelihood_NativeAminoAcidLikelihoodCore
 * Method:    nativeStatesStatesPruning
 * Signature: ([I[D[I[DII[D)V
 */

/**
	 Calculates partial likelihoods at a node when both children have states.
*/

JNIEXPORT void JNICALL Java_dr_evomodel_treelikelihood_NativeAminoAcidLikelihoodCore_nativeStatesStatesPruning
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
							
			if (state1 < 20 && state2 < 20) {
				partials3[v] = matrices1[w + state1] * matrices2[w + state2];
				v++;	w += 20;
				partials3[v] = matrices1[w + state1] * matrices2[w + state2];
				v++;	w += 20;
				partials3[v] = matrices1[w + state1] * matrices2[w + state2];
				v++;	w += 20;
				partials3[v] = matrices1[w + state1] * matrices2[w + state2];
				v++;	w += 20;

				partials3[v] = matrices1[w + state1] * matrices2[w + state2];
				v++;	w += 20;
				partials3[v] = matrices1[w + state1] * matrices2[w + state2];
				v++;	w += 20;
				partials3[v] = matrices1[w + state1] * matrices2[w + state2];
				v++;	w += 20;
				partials3[v] = matrices1[w + state1] * matrices2[w + state2];
				v++;	w += 20;

				partials3[v] = matrices1[w + state1] * matrices2[w + state2];
				v++;	w += 20;
				partials3[v] = matrices1[w + state1] * matrices2[w + state2];
				v++;	w += 20;
				partials3[v] = matrices1[w + state1] * matrices2[w + state2];
				v++;	w += 20;
				partials3[v] = matrices1[w + state1] * matrices2[w + state2];
				v++;	w += 20;

				partials3[v] = matrices1[w + state1] * matrices2[w + state2];
				v++;	w += 20;
				partials3[v] = matrices1[w + state1] * matrices2[w + state2];
				v++;	w += 20;
				partials3[v] = matrices1[w + state1] * matrices2[w + state2];
				v++;	w += 20;
				partials3[v] = matrices1[w + state1] * matrices2[w + state2];
				v++;	w += 20;

				partials3[v] = matrices1[w + state1] * matrices2[w + state2];
				v++;	w += 20;
				partials3[v] = matrices1[w + state1] * matrices2[w + state2];
				v++;	w += 20;
				partials3[v] = matrices1[w + state1] * matrices2[w + state2];
				v++;	w += 20;
				partials3[v] = matrices1[w + state1] * matrices2[w + state2];
				v++;	w += 20;			
				
			} else if (state1 < 20) {
				partials3[v] = matrices1[w + state1];
				v++;	w += 20;
				partials3[v] = matrices1[w + state1];
				v++;	w += 20;
				partials3[v] = matrices1[w + state1];
				v++;	w += 20;
				partials3[v] = matrices1[w + state1];
				v++;	w += 20;
				
				partials3[v] = matrices1[w + state1];
				v++;	w += 20;
				partials3[v] = matrices1[w + state1];
				v++;	w += 20;
				partials3[v] = matrices1[w + state1];
				v++;	w += 20;
				partials3[v] = matrices1[w + state1];
				v++;	w += 20;
	
				partials3[v] = matrices1[w + state1];
				v++;	w += 20;
				partials3[v] = matrices1[w + state1];
				v++;	w += 20;
				partials3[v] = matrices1[w + state1];
				v++;	w += 20;
				partials3[v] = matrices1[w + state1];
				v++;	w += 20;

				partials3[v] = matrices1[w + state1];
				v++;	w += 20;
				partials3[v] = matrices1[w + state1];
				v++;	w += 20;
				partials3[v] = matrices1[w + state1];
				v++;	w += 20;
				partials3[v] = matrices1[w + state1];
				v++;	w += 20;

				partials3[v] = matrices1[w + state1];
				v++;	w += 20;
				partials3[v] = matrices1[w + state1];
				v++;	w += 20;
				partials3[v] = matrices1[w + state1];
				v++;	w += 20;
				partials3[v] = matrices1[w + state1];
				v++;	w += 20;
				
			} else if (state2 < 20) {
				partials3[v] = matrices2[w + state2];
				v++;	w += 20;
				partials3[v] = matrices2[w + state2];
				v++;	w += 20;
				partials3[v] = matrices2[w + state2];
				v++;	w += 20;
				partials3[v] = matrices2[w + state2];
				v++;	w += 20;
				
				partials3[v] = matrices2[w + state2];
				v++;	w += 20;
				partials3[v] = matrices2[w + state2];
				v++;	w += 20;
				partials3[v] = matrices2[w + state2];
				v++;	w += 20;
				partials3[v] = matrices2[w + state2];
				v++;	w += 20;
				
				partials3[v] = matrices2[w + state2];
				v++;	w += 20;
				partials3[v] = matrices2[w + state2];
				v++;	w += 20;
				partials3[v] = matrices2[w + state2];
				v++;	w += 20;
				partials3[v] = matrices2[w + state2];
				v++;	w += 20;
				
				partials3[v] = matrices2[w + state2];
				v++;	w += 20;
				partials3[v] = matrices2[w + state2];
				v++;	w += 20;
				partials3[v] = matrices2[w + state2];
				v++;	w += 20;
				partials3[v] = matrices2[w + state2];
				v++;	w += 20;
				
				partials3[v] = matrices2[w + state2];
				v++;	w += 20;
				partials3[v] = matrices2[w + state2];
				v++;	w += 20;
				partials3[v] = matrices2[w + state2];
				v++;	w += 20;
				partials3[v] = matrices2[w + state2];
				v++;	w += 20;
				
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
	
//	for(v = 0; v < 10; v++)
//		printf("%i\n",MATRIX_SIZE);

	(*env)->ReleasePrimitiveArrayCritical(env, outPartials, partials3, 0);
	(*env)->ReleasePrimitiveArrayCritical(env, inMatrices2, matrices2, JNI_ABORT);
	(*env)->ReleasePrimitiveArrayCritical(env, inStates2, states2, JNI_ABORT);
	(*env)->ReleasePrimitiveArrayCritical(env, inMatrices1, matrices1, JNI_ABORT);
	(*env)->ReleasePrimitiveArrayCritical(env, inStates1, states1, JNI_ABORT);
}

/*
 * Class:     dr_evomodel_treelikelihood_NativeAminoAcidLikelihoodCore
 * Method:    nativeStatesPartialsPruning
 * Signature: ([I[D[D[DII[D)V
 */
JNIEXPORT void JNICALL Java_dr_evomodel_treelikelihood_NativeAminoAcidLikelihoodCore_nativeStatesPartialsPruning
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
	
	int i, l, k, state1;
	double sum;
	int u = 0;
	int v = 0;
	int w = 0;
	int x, y;

	for (l = 0; l < matrixCount; l++) {
		for (k = 0; k < patternCount; k++) {
		
			int state1 = states1[k];
			
				if (state1 < 20) {
			
					double sum;
					
					x = w;
					for (i = 0; i < 20; i++) {
					
						y = v;
						double value = matrices1[x + state1];
						sum =	matrices2[x] * partials2[y]; x++; y++;
						sum +=	matrices2[x] * partials2[y]; x++; y++;
						sum +=	matrices2[x] * partials2[y]; x++; y++;
						sum +=	matrices2[x] * partials2[y]; x++; y++;
						sum +=	matrices2[x] * partials2[y]; x++; y++;
						sum +=	matrices2[x] * partials2[y]; x++; y++;
						sum +=	matrices2[x] * partials2[y]; x++; y++;
						sum +=	matrices2[x] * partials2[y]; x++; y++;
						sum +=	matrices2[x] * partials2[y]; x++; y++;
						sum +=	matrices2[x] * partials2[y]; x++; y++;
						sum +=	matrices2[x] * partials2[y]; x++; y++;
						sum +=	matrices2[x] * partials2[y]; x++; y++;
						sum +=	matrices2[x] * partials2[y]; x++; y++;
						sum +=	matrices2[x] * partials2[y]; x++; y++;
						sum +=	matrices2[x] * partials2[y]; x++; y++;
						sum +=	matrices2[x] * partials2[y]; x++; y++;
						sum +=	matrices2[x] * partials2[y]; x++; y++;
						sum +=	matrices2[x] * partials2[y]; x++; y++;
						sum +=	matrices2[x] * partials2[y]; x++; y++;
						sum +=	matrices2[x] * partials2[y]; x++; y++;
						partials3[u] = value * sum;	
						u++;
					}
											
					v += 20;

				} else {
					// Child 1 has a gap or unknown state so don't use it

					double sum;
					
					x = w;
					for (i = 0; i < 20; i++) {
					
						y = v;
						sum =	matrices2[x] * partials2[y]; x++; y++;
						sum +=	matrices2[x] * partials2[y]; x++; y++;
						sum +=	matrices2[x] * partials2[y]; x++; y++;
						sum +=	matrices2[x] * partials2[y]; x++; y++;
						sum +=	matrices2[x] * partials2[y]; x++; y++;
						sum +=	matrices2[x] * partials2[y]; x++; y++;
						sum +=	matrices2[x] * partials2[y]; x++; y++;
						sum +=	matrices2[x] * partials2[y]; x++; y++;
						sum +=	matrices2[x] * partials2[y]; x++; y++;
						sum +=	matrices2[x] * partials2[y]; x++; y++;
						sum +=	matrices2[x] * partials2[y]; x++; y++;
						sum +=	matrices2[x] * partials2[y]; x++; y++;
						sum +=	matrices2[x] * partials2[y]; x++; y++;
						sum +=	matrices2[x] * partials2[y]; x++; y++;
						sum +=	matrices2[x] * partials2[y]; x++; y++;
						sum +=	matrices2[x] * partials2[y]; x++; y++;
						sum +=	matrices2[x] * partials2[y]; x++; y++;
						sum +=	matrices2[x] * partials2[y]; x++; y++;
						sum +=	matrices2[x] * partials2[y]; x++; y++;
						sum +=	matrices2[x] * partials2[y]; x++; y++;
						partials3[u] = sum;	
						u++;
					}
											
					v += 20;
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
 * Class:     dr_evomodel_treelikelihood_NativeAminoAcidLikelihoodCore
 * Method:    nativePartialsPartialsPruning
 * Signature: ([D[D[D[DII[D)V
 */
JNIEXPORT void JNICALL Java_dr_evomodel_treelikelihood_NativeAminoAcidLikelihoodCore_nativePartialsPartialsPruning
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
	int x,y;

	for (l = 0; l < matrixCount; l++) {
		for (k = 0; k < patternCount; k++) {
		
				x = w;
				for (i = 0; i < 20; i++) {
					y = v;
					sum1 =	matrices1[x] * partials1[y];
					sum2 =	matrices2[x] * partials2[y]; x++; y++;
					sum1 +=	matrices1[x] * partials1[y];
					sum2 +=	matrices2[x] * partials2[y]; x++; y++;
					sum1 +=	matrices1[x] * partials1[y];
					sum2 +=	matrices2[x] * partials2[y]; x++; y++;
					sum1 +=	matrices1[x] * partials1[y];
					sum2 +=	matrices2[x] * partials2[y]; x++; y++;
					sum1 +=	matrices1[x] * partials1[y];
					sum2 +=	matrices2[x] * partials2[y]; x++; y++;
					sum1 +=	matrices1[x] * partials1[y];
					sum2 +=	matrices2[x] * partials2[y]; x++; y++;
					sum1 +=	matrices1[x] * partials1[y];
					sum2 +=	matrices2[x] * partials2[y]; x++; y++;
					sum1 +=	matrices1[x] * partials1[y];
					sum2 +=	matrices2[x] * partials2[y]; x++; y++;
					sum1 +=	matrices1[x] * partials1[y];
					sum2 +=	matrices2[x] * partials2[y]; x++; y++;
					sum1 +=	matrices1[x] * partials1[y];
					sum2 +=	matrices2[x] * partials2[y]; x++; y++;
					sum1 +=	matrices1[x] * partials1[y];
					sum2 +=	matrices2[x] * partials2[y]; x++; y++;
					sum1 +=	matrices1[x] * partials1[y];
					sum2 +=	matrices2[x] * partials2[y]; x++; y++;
					sum1 +=	matrices1[x] * partials1[y];
					sum2 +=	matrices2[x] * partials2[y]; x++; y++;
					sum1 +=	matrices1[x] * partials1[y];
					sum2 +=	matrices2[x] * partials2[y]; x++; y++;
					sum1 +=	matrices1[x] * partials1[y];
					sum2 +=	matrices2[x] * partials2[y]; x++; y++;
					sum1 +=	matrices1[x] * partials1[y];
					sum2 +=	matrices2[x] * partials2[y]; x++; y++;
					sum1 +=	matrices1[x] * partials1[y];
					sum2 +=	matrices2[x] * partials2[y]; x++; y++;
					sum1 +=	matrices1[x] * partials1[y];
					sum2 +=	matrices2[x] * partials2[y]; x++; y++;
					sum1 +=	matrices1[x] * partials1[y];
					sum2 +=	matrices2[x] * partials2[y]; x++; y++;
					sum1 +=	matrices1[x] * partials1[y];
					sum2 +=	matrices2[x] * partials2[y]; x++; y++;

					partials3[u] = sum1 * sum2;	
					u++;
				}
										
				v += 20;
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
 * Class:     dr_evomodel_treelikelihood_NativeAminoAcidLikelihoodCore
 * Method:    nativeIntegratePartials
 * Signature: ([D[DII[D)V
 */
JNIEXPORT void JNICALL Java_dr_evomodel_treelikelihood_NativeAminoAcidLikelihoodCore_nativeIntegratePartials
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
		
		partials3[u] = partials1[v] * proportions[0]; u++; v++;
		partials3[u] = partials1[v] * proportions[0]; u++; v++;
		partials3[u] = partials1[v] * proportions[0]; u++; v++;
		partials3[u] = partials1[v] * proportions[0]; u++; v++;

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

			partials3[u] += partials1[v] * proportions[j]; u++; v++;
			partials3[u] += partials1[v] * proportions[j]; u++; v++;
			partials3[u] += partials1[v] * proportions[j]; u++; v++;
			partials3[u] += partials1[v] * proportions[j]; u++; v++;

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
							
