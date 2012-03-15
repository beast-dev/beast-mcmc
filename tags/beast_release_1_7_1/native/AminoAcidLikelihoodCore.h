#ifndef _Included_dr_evomodel_treelikelihood_NativeAminoAcidLikelihoodCore
#define _Included_dr_evomodel_treelikelihood_NativeAminoAcidLikelihoodCore
#include <jni.h>
#ifdef __cplusplus
extern "C" {
#endif

/**
	 NativeAminoAcidLikelihoodCore - An implementation of LikelihoodCore for AminoAcids 
	 that calls native methods for maximum speed. The native methods should be in
	 a shared library called "AminoAcidLikelihoodCore" but the exact name will be system
	 dependent (i.e. "libAminoAcidLikelihoodCore.so" or "AminoAcidLikelihoodCore.dll").
	 @version $Id: AminoAcidLikelihoodCore.h,v 1.2 2003/10/28 16:26:31 rambaut Exp $
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
	(JNIEnv *, jobject, jintArray, jdoubleArray, jintArray, jdoubleArray, jint, jint, jdoubleArray);

/*
 * Class:     dr_evomodel_treelikelihood_NativeAminoAcidLikelihoodCore
 * Method:    nativeStatesPartialsPruning
 * Signature: ([I[D[D[DII[D)V
 */
JNIEXPORT void JNICALL Java_dr_evomodel_treelikelihood_NativeAminoAcidLikelihoodCore_nativeStatesPartialsPruning
	(JNIEnv *, jobject, jintArray, jdoubleArray, jdoubleArray, jdoubleArray, jint, jint, jdoubleArray);

/*
 * Class:     dr_evomodel_treelikelihood_NativeAminoAcidLikelihoodCore
 * Method:    nativePartialsPartialsPruning
 * Signature: ([D[D[D[DII[D)V
 */
JNIEXPORT void JNICALL Java_dr_evomodel_treelikelihood_NativeAminoAcidLikelihoodCore_nativePartialsPartialsPruning
	(JNIEnv *, jobject, jdoubleArray, jdoubleArray, jdoubleArray, jdoubleArray, jint, jint, jdoubleArray);

/*
 * Class:     dr_evomodel_treelikelihood_NativeAminoAcidLikelihoodCore
 * Method:    nativeIntegratePartials
 * Signature: ([D[DII[D)V
 */
JNIEXPORT void JNICALL Java_dr_evomodel_treelikelihood_NativeAminoAcidLikelihoodCore_nativeIntegratePartials
	(JNIEnv *, jobject, jdoubleArray, jdoubleArray, jint, jint, jdoubleArray);

#ifdef __cplusplus
}
#endif
#endif
