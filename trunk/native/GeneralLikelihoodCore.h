#ifndef _Included_dr_evomodel_treelikelihood_GeneralLikelihoodCore
#define _Included_dr_evomodel_treelikelihood_GeneralLikelihoodCore
#include <jni.h>
#ifdef __cplusplus
extern "C" {
#endif

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
	(JNIEnv *, jobject, jintArray, jdoubleArray, jintArray, jdoubleArray, jint, jint, jdoubleArray, jint);

/*
 * Class:     dr_evomodel_treelikelihood_NativeGeneralLikelihoodCore
 * Method:    nativeStatesPartialsPruning
 * Signature: ([I[D[D[DII[D)V
 */
JNIEXPORT void JNICALL Java_dr_evomodel_treelikelihood_NativeGeneralLikelihoodCore_nativeStatesPartialsPruning
	(JNIEnv *, jobject, jintArray, jdoubleArray, jdoubleArray, jdoubleArray, jint, jint, jdoubleArray, jint);

/*
 * Class:     dr_evomodel_treelikelihood_NativeGeneralLikelihoodCore
 * Method:    nativePartialsPartialsPruning
 * Signature: ([D[D[D[DII[D)V
 */
JNIEXPORT void JNICALL Java_dr_evomodel_treelikelihood_NativeGeneralLikelihoodCore_nativePartialsPartialsPruning
	(JNIEnv *, jobject, jdoubleArray, jdoubleArray, jdoubleArray, jdoubleArray, jint, jint, jdoubleArray, jint);

/*
 * Class:     dr_evomodel_treelikelihood_NativeGeneralLikelihoodCore
 * Method:    nativeIntegratePartials
 * Signature: ([D[DII[D)V
 */
JNIEXPORT void JNICALL Java_dr_evomodel_treelikelihood_NativeGeneralLikelihoodCore_nativeIntegratePartials
	(JNIEnv *, jobject, jdoubleArray, jdoubleArray, jint, jint, jdoubleArray, jint);
	
#ifdef __cplusplus
}
#endif
#endif
