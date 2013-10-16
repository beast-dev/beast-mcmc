
#include "AminoAcidLikelihoodCore.h"

#define STATE_COUNT 20
#define MATRIX_SIZE STATE_COUNT * STATE_COUNT

/**
	 AminoAcidLikelihoodCore - An implementation of LikelihoodCore for AminoAcids 
	 that calls native methods for maximum speed. The native methods should be in
	 a shared library called "AminoAcidLikelihoodCore" but the exact name will be system
	 dependent (i.e. "libAminoAcidLikelihoodCore.so" or "AminoAcidLikelihoodCore.dll").
	 @version $Id: AminoAcidLikelihoodCore.c,v 1.4 2004/03/16 17:56:35 rambaut Exp $
	 @author Erik Bloomquist
	 
	 
	 June 16, 2008
	 
	 Note that the methods
	 1. Java_dr_evomodel_treelikelihood_NativeAminoAcidLikelihoodCore_nativeStatesPartialsPruning
	 2. Java_dr_evomodel_treelikelihood_NativeAminoAcidLikelihoodCore_nativePartialsPartialsPruning
	 
	 have been unrolled.  Seems to be slightly faster than the loop versions,
	 but only by a mild amount 3-5%.  The looped versions have been commented out
	 and are at the bottom of the code.
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
	int u = 0;
	int v = 0;
	int w = 0;
	double sum;
	

	for (l = 0; l < matrixCount; l++) {
		for (k = 0; k < patternCount; k++) {
			state1 = states1[k];

				if (state1 < 20) {

				 sum = 0;

sum = matrices2[w] * partials2[v];
sum += matrices2[w + 1] * partials2[v + 1];
sum += matrices2[w + 2] * partials2[v + 2];
sum += matrices2[w + 3] * partials2[v + 3];
sum += matrices2[w + 4] * partials2[v + 4];
sum += matrices2[w + 5] * partials2[v + 5];
sum += matrices2[w + 6] * partials2[v + 6];
sum += matrices2[w + 7] * partials2[v + 7];
sum += matrices2[w + 8] * partials2[v + 8];
sum += matrices2[w + 9] * partials2[v + 9];
sum += matrices2[w + 10] * partials2[v + 10];
sum += matrices2[w + 11] * partials2[v + 11];
sum += matrices2[w + 12] * partials2[v + 12];
sum += matrices2[w + 13] * partials2[v + 13];
sum += matrices2[w + 14] * partials2[v + 14];
sum += matrices2[w + 15] * partials2[v + 15];
sum += matrices2[w + 16] * partials2[v + 16];
sum += matrices2[w + 17] * partials2[v + 17];
sum += matrices2[w + 18] * partials2[v + 18];
sum += matrices2[w + 19] * partials2[v + 19];
partials3[u] = matrices1[w + state1] * sum;	u++;


sum = matrices2[w + 20] * partials2[v];
sum += matrices2[w + 21] * partials2[v + 1];
sum += matrices2[w + 22] * partials2[v + 2];
sum += matrices2[w + 23] * partials2[v + 3];
sum += matrices2[w + 24] * partials2[v + 4];
sum += matrices2[w + 25] * partials2[v + 5];
sum += matrices2[w + 26] * partials2[v + 6];
sum += matrices2[w + 27] * partials2[v + 7];
sum += matrices2[w + 28] * partials2[v + 8];
sum += matrices2[w + 29] * partials2[v + 9];
sum += matrices2[w + 30] * partials2[v + 10];
sum += matrices2[w + 31] * partials2[v + 11];
sum += matrices2[w + 32] * partials2[v + 12];
sum += matrices2[w + 33] * partials2[v + 13];
sum += matrices2[w + 34] * partials2[v + 14];
sum += matrices2[w + 35] * partials2[v + 15];
sum += matrices2[w + 36] * partials2[v + 16];
sum += matrices2[w + 37] * partials2[v + 17];
sum += matrices2[w + 38] * partials2[v + 18];
sum += matrices2[w + 39] * partials2[v + 19];
partials3[u] = matrices1[w + 20+ state1] * sum;	u++;


sum = matrices2[w + 40] * partials2[v];
sum += matrices2[w + 41] * partials2[v + 1];
sum += matrices2[w + 42] * partials2[v + 2];
sum += matrices2[w + 43] * partials2[v + 3];
sum += matrices2[w + 44] * partials2[v + 4];
sum += matrices2[w + 45] * partials2[v + 5];
sum += matrices2[w + 46] * partials2[v + 6];
sum += matrices2[w + 47] * partials2[v + 7];
sum += matrices2[w + 48] * partials2[v + 8];
sum += matrices2[w + 49] * partials2[v + 9];
sum += matrices2[w + 50] * partials2[v + 10];
sum += matrices2[w + 51] * partials2[v + 11];
sum += matrices2[w + 52] * partials2[v + 12];
sum += matrices2[w + 53] * partials2[v + 13];
sum += matrices2[w + 54] * partials2[v + 14];
sum += matrices2[w + 55] * partials2[v + 15];
sum += matrices2[w + 56] * partials2[v + 16];
sum += matrices2[w + 57] * partials2[v + 17];
sum += matrices2[w + 58] * partials2[v + 18];
sum += matrices2[w + 59] * partials2[v + 19];
partials3[u] = matrices1[w + 40+ state1] * sum;	u++;


sum = matrices2[w + 60] * partials2[v];
sum += matrices2[w + 61] * partials2[v + 1];
sum += matrices2[w + 62] * partials2[v + 2];
sum += matrices2[w + 63] * partials2[v + 3];
sum += matrices2[w + 64] * partials2[v + 4];
sum += matrices2[w + 65] * partials2[v + 5];
sum += matrices2[w + 66] * partials2[v + 6];
sum += matrices2[w + 67] * partials2[v + 7];
sum += matrices2[w + 68] * partials2[v + 8];
sum += matrices2[w + 69] * partials2[v + 9];
sum += matrices2[w + 70] * partials2[v + 10];
sum += matrices2[w + 71] * partials2[v + 11];
sum += matrices2[w + 72] * partials2[v + 12];
sum += matrices2[w + 73] * partials2[v + 13];
sum += matrices2[w + 74] * partials2[v + 14];
sum += matrices2[w + 75] * partials2[v + 15];
sum += matrices2[w + 76] * partials2[v + 16];
sum += matrices2[w + 77] * partials2[v + 17];
sum += matrices2[w + 78] * partials2[v + 18];
sum += matrices2[w + 79] * partials2[v + 19];
partials3[u] = matrices1[w + 60+ state1] * sum;	u++;


sum = matrices2[w + 80] * partials2[v];
sum += matrices2[w + 81] * partials2[v + 1];
sum += matrices2[w + 82] * partials2[v + 2];
sum += matrices2[w + 83] * partials2[v + 3];
sum += matrices2[w + 84] * partials2[v + 4];
sum += matrices2[w + 85] * partials2[v + 5];
sum += matrices2[w + 86] * partials2[v + 6];
sum += matrices2[w + 87] * partials2[v + 7];
sum += matrices2[w + 88] * partials2[v + 8];
sum += matrices2[w + 89] * partials2[v + 9];
sum += matrices2[w + 90] * partials2[v + 10];
sum += matrices2[w + 91] * partials2[v + 11];
sum += matrices2[w + 92] * partials2[v + 12];
sum += matrices2[w + 93] * partials2[v + 13];
sum += matrices2[w + 94] * partials2[v + 14];
sum += matrices2[w + 95] * partials2[v + 15];
sum += matrices2[w + 96] * partials2[v + 16];
sum += matrices2[w + 97] * partials2[v + 17];
sum += matrices2[w + 98] * partials2[v + 18];
sum += matrices2[w + 99] * partials2[v + 19];
partials3[u] = matrices1[w + 80+ state1] * sum;	u++;


sum = matrices2[w + 100] * partials2[v];
sum += matrices2[w + 101] * partials2[v + 1];
sum += matrices2[w + 102] * partials2[v + 2];
sum += matrices2[w + 103] * partials2[v + 3];
sum += matrices2[w + 104] * partials2[v + 4];
sum += matrices2[w + 105] * partials2[v + 5];
sum += matrices2[w + 106] * partials2[v + 6];
sum += matrices2[w + 107] * partials2[v + 7];
sum += matrices2[w + 108] * partials2[v + 8];
sum += matrices2[w + 109] * partials2[v + 9];
sum += matrices2[w + 110] * partials2[v + 10];
sum += matrices2[w + 111] * partials2[v + 11];
sum += matrices2[w + 112] * partials2[v + 12];
sum += matrices2[w + 113] * partials2[v + 13];
sum += matrices2[w + 114] * partials2[v + 14];
sum += matrices2[w + 115] * partials2[v + 15];
sum += matrices2[w + 116] * partials2[v + 16];
sum += matrices2[w + 117] * partials2[v + 17];
sum += matrices2[w + 118] * partials2[v + 18];
sum += matrices2[w + 119] * partials2[v + 19];
partials3[u] = matrices1[w + 100+ state1] * sum;	u++;


sum = matrices2[w + 120] * partials2[v];
sum += matrices2[w + 121] * partials2[v + 1];
sum += matrices2[w + 122] * partials2[v + 2];
sum += matrices2[w + 123] * partials2[v + 3];
sum += matrices2[w + 124] * partials2[v + 4];
sum += matrices2[w + 125] * partials2[v + 5];
sum += matrices2[w + 126] * partials2[v + 6];
sum += matrices2[w + 127] * partials2[v + 7];
sum += matrices2[w + 128] * partials2[v + 8];
sum += matrices2[w + 129] * partials2[v + 9];
sum += matrices2[w + 130] * partials2[v + 10];
sum += matrices2[w + 131] * partials2[v + 11];
sum += matrices2[w + 132] * partials2[v + 12];
sum += matrices2[w + 133] * partials2[v + 13];
sum += matrices2[w + 134] * partials2[v + 14];
sum += matrices2[w + 135] * partials2[v + 15];
sum += matrices2[w + 136] * partials2[v + 16];
sum += matrices2[w + 137] * partials2[v + 17];
sum += matrices2[w + 138] * partials2[v + 18];
sum += matrices2[w + 139] * partials2[v + 19];
partials3[u] = matrices1[w + 120+ state1] * sum;	u++;


sum = matrices2[w + 140] * partials2[v];
sum += matrices2[w + 141] * partials2[v + 1];
sum += matrices2[w + 142] * partials2[v + 2];
sum += matrices2[w + 143] * partials2[v + 3];
sum += matrices2[w + 144] * partials2[v + 4];
sum += matrices2[w + 145] * partials2[v + 5];
sum += matrices2[w + 146] * partials2[v + 6];
sum += matrices2[w + 147] * partials2[v + 7];
sum += matrices2[w + 148] * partials2[v + 8];
sum += matrices2[w + 149] * partials2[v + 9];
sum += matrices2[w + 150] * partials2[v + 10];
sum += matrices2[w + 151] * partials2[v + 11];
sum += matrices2[w + 152] * partials2[v + 12];
sum += matrices2[w + 153] * partials2[v + 13];
sum += matrices2[w + 154] * partials2[v + 14];
sum += matrices2[w + 155] * partials2[v + 15];
sum += matrices2[w + 156] * partials2[v + 16];
sum += matrices2[w + 157] * partials2[v + 17];
sum += matrices2[w + 158] * partials2[v + 18];
sum += matrices2[w + 159] * partials2[v + 19];
partials3[u] = matrices1[w + 140+ state1] * sum;	u++;


sum = matrices2[w + 160] * partials2[v];
sum += matrices2[w + 161] * partials2[v + 1];
sum += matrices2[w + 162] * partials2[v + 2];
sum += matrices2[w + 163] * partials2[v + 3];
sum += matrices2[w + 164] * partials2[v + 4];
sum += matrices2[w + 165] * partials2[v + 5];
sum += matrices2[w + 166] * partials2[v + 6];
sum += matrices2[w + 167] * partials2[v + 7];
sum += matrices2[w + 168] * partials2[v + 8];
sum += matrices2[w + 169] * partials2[v + 9];
sum += matrices2[w + 170] * partials2[v + 10];
sum += matrices2[w + 171] * partials2[v + 11];
sum += matrices2[w + 172] * partials2[v + 12];
sum += matrices2[w + 173] * partials2[v + 13];
sum += matrices2[w + 174] * partials2[v + 14];
sum += matrices2[w + 175] * partials2[v + 15];
sum += matrices2[w + 176] * partials2[v + 16];
sum += matrices2[w + 177] * partials2[v + 17];
sum += matrices2[w + 178] * partials2[v + 18];
sum += matrices2[w + 179] * partials2[v + 19];
partials3[u] = matrices1[w + 160+ state1] * sum;	u++;


sum = matrices2[w + 180] * partials2[v];
sum += matrices2[w + 181] * partials2[v + 1];
sum += matrices2[w + 182] * partials2[v + 2];
sum += matrices2[w + 183] * partials2[v + 3];
sum += matrices2[w + 184] * partials2[v + 4];
sum += matrices2[w + 185] * partials2[v + 5];
sum += matrices2[w + 186] * partials2[v + 6];
sum += matrices2[w + 187] * partials2[v + 7];
sum += matrices2[w + 188] * partials2[v + 8];
sum += matrices2[w + 189] * partials2[v + 9];
sum += matrices2[w + 190] * partials2[v + 10];
sum += matrices2[w + 191] * partials2[v + 11];
sum += matrices2[w + 192] * partials2[v + 12];
sum += matrices2[w + 193] * partials2[v + 13];
sum += matrices2[w + 194] * partials2[v + 14];
sum += matrices2[w + 195] * partials2[v + 15];
sum += matrices2[w + 196] * partials2[v + 16];
sum += matrices2[w + 197] * partials2[v + 17];
sum += matrices2[w + 198] * partials2[v + 18];
sum += matrices2[w + 199] * partials2[v + 19];
partials3[u] = matrices1[w + 180+ state1] * sum;	u++;


sum = matrices2[w + 200] * partials2[v];
sum += matrices2[w + 201] * partials2[v + 1];
sum += matrices2[w + 202] * partials2[v + 2];
sum += matrices2[w + 203] * partials2[v + 3];
sum += matrices2[w + 204] * partials2[v + 4];
sum += matrices2[w + 205] * partials2[v + 5];
sum += matrices2[w + 206] * partials2[v + 6];
sum += matrices2[w + 207] * partials2[v + 7];
sum += matrices2[w + 208] * partials2[v + 8];
sum += matrices2[w + 209] * partials2[v + 9];
sum += matrices2[w + 210] * partials2[v + 10];
sum += matrices2[w + 211] * partials2[v + 11];
sum += matrices2[w + 212] * partials2[v + 12];
sum += matrices2[w + 213] * partials2[v + 13];
sum += matrices2[w + 214] * partials2[v + 14];
sum += matrices2[w + 215] * partials2[v + 15];
sum += matrices2[w + 216] * partials2[v + 16];
sum += matrices2[w + 217] * partials2[v + 17];
sum += matrices2[w + 218] * partials2[v + 18];
sum += matrices2[w + 219] * partials2[v + 19];
partials3[u] = matrices1[w + 200+ state1] * sum;	u++;


sum = matrices2[w + 220] * partials2[v];
sum += matrices2[w + 221] * partials2[v + 1];
sum += matrices2[w + 222] * partials2[v + 2];
sum += matrices2[w + 223] * partials2[v + 3];
sum += matrices2[w + 224] * partials2[v + 4];
sum += matrices2[w + 225] * partials2[v + 5];
sum += matrices2[w + 226] * partials2[v + 6];
sum += matrices2[w + 227] * partials2[v + 7];
sum += matrices2[w + 228] * partials2[v + 8];
sum += matrices2[w + 229] * partials2[v + 9];
sum += matrices2[w + 230] * partials2[v + 10];
sum += matrices2[w + 231] * partials2[v + 11];
sum += matrices2[w + 232] * partials2[v + 12];
sum += matrices2[w + 233] * partials2[v + 13];
sum += matrices2[w + 234] * partials2[v + 14];
sum += matrices2[w + 235] * partials2[v + 15];
sum += matrices2[w + 236] * partials2[v + 16];
sum += matrices2[w + 237] * partials2[v + 17];
sum += matrices2[w + 238] * partials2[v + 18];
sum += matrices2[w + 239] * partials2[v + 19];
partials3[u] = matrices1[w + 220+ state1] * sum;	u++;


sum = matrices2[w + 240] * partials2[v];
sum += matrices2[w + 241] * partials2[v + 1];
sum += matrices2[w + 242] * partials2[v + 2];
sum += matrices2[w + 243] * partials2[v + 3];
sum += matrices2[w + 244] * partials2[v + 4];
sum += matrices2[w + 245] * partials2[v + 5];
sum += matrices2[w + 246] * partials2[v + 6];
sum += matrices2[w + 247] * partials2[v + 7];
sum += matrices2[w + 248] * partials2[v + 8];
sum += matrices2[w + 249] * partials2[v + 9];
sum += matrices2[w + 250] * partials2[v + 10];
sum += matrices2[w + 251] * partials2[v + 11];
sum += matrices2[w + 252] * partials2[v + 12];
sum += matrices2[w + 253] * partials2[v + 13];
sum += matrices2[w + 254] * partials2[v + 14];
sum += matrices2[w + 255] * partials2[v + 15];
sum += matrices2[w + 256] * partials2[v + 16];
sum += matrices2[w + 257] * partials2[v + 17];
sum += matrices2[w + 258] * partials2[v + 18];
sum += matrices2[w + 259] * partials2[v + 19];
partials3[u] = matrices1[w + 240+ state1] * sum;	u++;


sum = matrices2[w + 260] * partials2[v];
sum += matrices2[w + 261] * partials2[v + 1];
sum += matrices2[w + 262] * partials2[v + 2];
sum += matrices2[w + 263] * partials2[v + 3];
sum += matrices2[w + 264] * partials2[v + 4];
sum += matrices2[w + 265] * partials2[v + 5];
sum += matrices2[w + 266] * partials2[v + 6];
sum += matrices2[w + 267] * partials2[v + 7];
sum += matrices2[w + 268] * partials2[v + 8];
sum += matrices2[w + 269] * partials2[v + 9];
sum += matrices2[w + 270] * partials2[v + 10];
sum += matrices2[w + 271] * partials2[v + 11];
sum += matrices2[w + 272] * partials2[v + 12];
sum += matrices2[w + 273] * partials2[v + 13];
sum += matrices2[w + 274] * partials2[v + 14];
sum += matrices2[w + 275] * partials2[v + 15];
sum += matrices2[w + 276] * partials2[v + 16];
sum += matrices2[w + 277] * partials2[v + 17];
sum += matrices2[w + 278] * partials2[v + 18];
sum += matrices2[w + 279] * partials2[v + 19];
partials3[u] = matrices1[w + 260+ state1] * sum;	u++;


sum = matrices2[w + 280] * partials2[v];
sum += matrices2[w + 281] * partials2[v + 1];
sum += matrices2[w + 282] * partials2[v + 2];
sum += matrices2[w + 283] * partials2[v + 3];
sum += matrices2[w + 284] * partials2[v + 4];
sum += matrices2[w + 285] * partials2[v + 5];
sum += matrices2[w + 286] * partials2[v + 6];
sum += matrices2[w + 287] * partials2[v + 7];
sum += matrices2[w + 288] * partials2[v + 8];
sum += matrices2[w + 289] * partials2[v + 9];
sum += matrices2[w + 290] * partials2[v + 10];
sum += matrices2[w + 291] * partials2[v + 11];
sum += matrices2[w + 292] * partials2[v + 12];
sum += matrices2[w + 293] * partials2[v + 13];
sum += matrices2[w + 294] * partials2[v + 14];
sum += matrices2[w + 295] * partials2[v + 15];
sum += matrices2[w + 296] * partials2[v + 16];
sum += matrices2[w + 297] * partials2[v + 17];
sum += matrices2[w + 298] * partials2[v + 18];
sum += matrices2[w + 299] * partials2[v + 19];
partials3[u] = matrices1[w + 280+ state1] * sum;	u++;


sum = matrices2[w + 300] * partials2[v];
sum += matrices2[w + 301] * partials2[v + 1];
sum += matrices2[w + 302] * partials2[v + 2];
sum += matrices2[w + 303] * partials2[v + 3];
sum += matrices2[w + 304] * partials2[v + 4];
sum += matrices2[w + 305] * partials2[v + 5];
sum += matrices2[w + 306] * partials2[v + 6];
sum += matrices2[w + 307] * partials2[v + 7];
sum += matrices2[w + 308] * partials2[v + 8];
sum += matrices2[w + 309] * partials2[v + 9];
sum += matrices2[w + 310] * partials2[v + 10];
sum += matrices2[w + 311] * partials2[v + 11];
sum += matrices2[w + 312] * partials2[v + 12];
sum += matrices2[w + 313] * partials2[v + 13];
sum += matrices2[w + 314] * partials2[v + 14];
sum += matrices2[w + 315] * partials2[v + 15];
sum += matrices2[w + 316] * partials2[v + 16];
sum += matrices2[w + 317] * partials2[v + 17];
sum += matrices2[w + 318] * partials2[v + 18];
sum += matrices2[w + 319] * partials2[v + 19];
partials3[u] = matrices1[w + 300+ state1] * sum;	u++;


sum = matrices2[w + 320] * partials2[v];
sum += matrices2[w + 321] * partials2[v + 1];
sum += matrices2[w + 322] * partials2[v + 2];
sum += matrices2[w + 323] * partials2[v + 3];
sum += matrices2[w + 324] * partials2[v + 4];
sum += matrices2[w + 325] * partials2[v + 5];
sum += matrices2[w + 326] * partials2[v + 6];
sum += matrices2[w + 327] * partials2[v + 7];
sum += matrices2[w + 328] * partials2[v + 8];
sum += matrices2[w + 329] * partials2[v + 9];
sum += matrices2[w + 330] * partials2[v + 10];
sum += matrices2[w + 331] * partials2[v + 11];
sum += matrices2[w + 332] * partials2[v + 12];
sum += matrices2[w + 333] * partials2[v + 13];
sum += matrices2[w + 334] * partials2[v + 14];
sum += matrices2[w + 335] * partials2[v + 15];
sum += matrices2[w + 336] * partials2[v + 16];
sum += matrices2[w + 337] * partials2[v + 17];
sum += matrices2[w + 338] * partials2[v + 18];
sum += matrices2[w + 339] * partials2[v + 19];
partials3[u] = matrices1[w + 320+ state1] * sum;	u++;


sum = matrices2[w + 340] * partials2[v];
sum += matrices2[w + 341] * partials2[v + 1];
sum += matrices2[w + 342] * partials2[v + 2];
sum += matrices2[w + 343] * partials2[v + 3];
sum += matrices2[w + 344] * partials2[v + 4];
sum += matrices2[w + 345] * partials2[v + 5];
sum += matrices2[w + 346] * partials2[v + 6];
sum += matrices2[w + 347] * partials2[v + 7];
sum += matrices2[w + 348] * partials2[v + 8];
sum += matrices2[w + 349] * partials2[v + 9];
sum += matrices2[w + 350] * partials2[v + 10];
sum += matrices2[w + 351] * partials2[v + 11];
sum += matrices2[w + 352] * partials2[v + 12];
sum += matrices2[w + 353] * partials2[v + 13];
sum += matrices2[w + 354] * partials2[v + 14];
sum += matrices2[w + 355] * partials2[v + 15];
sum += matrices2[w + 356] * partials2[v + 16];
sum += matrices2[w + 357] * partials2[v + 17];
sum += matrices2[w + 358] * partials2[v + 18];
sum += matrices2[w + 359] * partials2[v + 19];
partials3[u] = matrices1[w + 340+ state1] * sum;	u++;


sum = matrices2[w + 360] * partials2[v];
sum += matrices2[w + 361] * partials2[v + 1];
sum += matrices2[w + 362] * partials2[v + 2];
sum += matrices2[w + 363] * partials2[v + 3];
sum += matrices2[w + 364] * partials2[v + 4];
sum += matrices2[w + 365] * partials2[v + 5];
sum += matrices2[w + 366] * partials2[v + 6];
sum += matrices2[w + 367] * partials2[v + 7];
sum += matrices2[w + 368] * partials2[v + 8];
sum += matrices2[w + 369] * partials2[v + 9];
sum += matrices2[w + 370] * partials2[v + 10];
sum += matrices2[w + 371] * partials2[v + 11];
sum += matrices2[w + 372] * partials2[v + 12];
sum += matrices2[w + 373] * partials2[v + 13];
sum += matrices2[w + 374] * partials2[v + 14];
sum += matrices2[w + 375] * partials2[v + 15];
sum += matrices2[w + 376] * partials2[v + 16];
sum += matrices2[w + 377] * partials2[v + 17];
sum += matrices2[w + 378] * partials2[v + 18];
sum += matrices2[w + 379] * partials2[v + 19];
partials3[u] = matrices1[w + 360+ state1] * sum;	u++;


sum = matrices2[w + 380] * partials2[v];
sum += matrices2[w + 381] * partials2[v + 1];
sum += matrices2[w + 382] * partials2[v + 2];
sum += matrices2[w + 383] * partials2[v + 3];
sum += matrices2[w + 384] * partials2[v + 4];
sum += matrices2[w + 385] * partials2[v + 5];
sum += matrices2[w + 386] * partials2[v + 6];
sum += matrices2[w + 387] * partials2[v + 7];
sum += matrices2[w + 388] * partials2[v + 8];
sum += matrices2[w + 389] * partials2[v + 9];
sum += matrices2[w + 390] * partials2[v + 10];
sum += matrices2[w + 391] * partials2[v + 11];
sum += matrices2[w + 392] * partials2[v + 12];
sum += matrices2[w + 393] * partials2[v + 13];
sum += matrices2[w + 394] * partials2[v + 14];
sum += matrices2[w + 395] * partials2[v + 15];
sum += matrices2[w + 396] * partials2[v + 16];
sum += matrices2[w + 397] * partials2[v + 17];
sum += matrices2[w + 398] * partials2[v + 18];
sum += matrices2[w + 399] * partials2[v + 19];
partials3[u] = matrices1[w + 380+ state1] * sum;	u++;

					v += 20;

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
sum += matrices2[w + 8] * partials2[v + 8];
sum += matrices2[w + 9] * partials2[v + 9];
sum += matrices2[w + 10] * partials2[v + 10];
sum += matrices2[w + 11] * partials2[v + 11];
sum += matrices2[w + 12] * partials2[v + 12];
sum += matrices2[w + 13] * partials2[v + 13];
sum += matrices2[w + 14] * partials2[v + 14];
sum += matrices2[w + 15] * partials2[v + 15];
sum += matrices2[w + 16] * partials2[v + 16];
sum += matrices2[w + 17] * partials2[v + 17];
sum += matrices2[w + 18] * partials2[v + 18];
sum += matrices2[w + 19] * partials2[v + 19];
partials3[u] =  sum;	u++;


sum = matrices2[w + 20] * partials2[v];
sum += matrices2[w + 21] * partials2[v + 1];
sum += matrices2[w + 22] * partials2[v + 2];
sum += matrices2[w + 23] * partials2[v + 3];
sum += matrices2[w + 24] * partials2[v + 4];
sum += matrices2[w + 25] * partials2[v + 5];
sum += matrices2[w + 26] * partials2[v + 6];
sum += matrices2[w + 27] * partials2[v + 7];
sum += matrices2[w + 28] * partials2[v + 8];
sum += matrices2[w + 29] * partials2[v + 9];
sum += matrices2[w + 30] * partials2[v + 10];
sum += matrices2[w + 31] * partials2[v + 11];
sum += matrices2[w + 32] * partials2[v + 12];
sum += matrices2[w + 33] * partials2[v + 13];
sum += matrices2[w + 34] * partials2[v + 14];
sum += matrices2[w + 35] * partials2[v + 15];
sum += matrices2[w + 36] * partials2[v + 16];
sum += matrices2[w + 37] * partials2[v + 17];
sum += matrices2[w + 38] * partials2[v + 18];
sum += matrices2[w + 39] * partials2[v + 19];
partials3[u] =  sum;	u++;


sum = matrices2[w + 40] * partials2[v];
sum += matrices2[w + 41] * partials2[v + 1];
sum += matrices2[w + 42] * partials2[v + 2];
sum += matrices2[w + 43] * partials2[v + 3];
sum += matrices2[w + 44] * partials2[v + 4];
sum += matrices2[w + 45] * partials2[v + 5];
sum += matrices2[w + 46] * partials2[v + 6];
sum += matrices2[w + 47] * partials2[v + 7];
sum += matrices2[w + 48] * partials2[v + 8];
sum += matrices2[w + 49] * partials2[v + 9];
sum += matrices2[w + 50] * partials2[v + 10];
sum += matrices2[w + 51] * partials2[v + 11];
sum += matrices2[w + 52] * partials2[v + 12];
sum += matrices2[w + 53] * partials2[v + 13];
sum += matrices2[w + 54] * partials2[v + 14];
sum += matrices2[w + 55] * partials2[v + 15];
sum += matrices2[w + 56] * partials2[v + 16];
sum += matrices2[w + 57] * partials2[v + 17];
sum += matrices2[w + 58] * partials2[v + 18];
sum += matrices2[w + 59] * partials2[v + 19];
partials3[u] =  sum;	u++;


sum = matrices2[w + 60] * partials2[v];
sum += matrices2[w + 61] * partials2[v + 1];
sum += matrices2[w + 62] * partials2[v + 2];
sum += matrices2[w + 63] * partials2[v + 3];
sum += matrices2[w + 64] * partials2[v + 4];
sum += matrices2[w + 65] * partials2[v + 5];
sum += matrices2[w + 66] * partials2[v + 6];
sum += matrices2[w + 67] * partials2[v + 7];
sum += matrices2[w + 68] * partials2[v + 8];
sum += matrices2[w + 69] * partials2[v + 9];
sum += matrices2[w + 70] * partials2[v + 10];
sum += matrices2[w + 71] * partials2[v + 11];
sum += matrices2[w + 72] * partials2[v + 12];
sum += matrices2[w + 73] * partials2[v + 13];
sum += matrices2[w + 74] * partials2[v + 14];
sum += matrices2[w + 75] * partials2[v + 15];
sum += matrices2[w + 76] * partials2[v + 16];
sum += matrices2[w + 77] * partials2[v + 17];
sum += matrices2[w + 78] * partials2[v + 18];
sum += matrices2[w + 79] * partials2[v + 19];
partials3[u] =  sum;	u++;


sum = matrices2[w + 80] * partials2[v];
sum += matrices2[w + 81] * partials2[v + 1];
sum += matrices2[w + 82] * partials2[v + 2];
sum += matrices2[w + 83] * partials2[v + 3];
sum += matrices2[w + 84] * partials2[v + 4];
sum += matrices2[w + 85] * partials2[v + 5];
sum += matrices2[w + 86] * partials2[v + 6];
sum += matrices2[w + 87] * partials2[v + 7];
sum += matrices2[w + 88] * partials2[v + 8];
sum += matrices2[w + 89] * partials2[v + 9];
sum += matrices2[w + 90] * partials2[v + 10];
sum += matrices2[w + 91] * partials2[v + 11];
sum += matrices2[w + 92] * partials2[v + 12];
sum += matrices2[w + 93] * partials2[v + 13];
sum += matrices2[w + 94] * partials2[v + 14];
sum += matrices2[w + 95] * partials2[v + 15];
sum += matrices2[w + 96] * partials2[v + 16];
sum += matrices2[w + 97] * partials2[v + 17];
sum += matrices2[w + 98] * partials2[v + 18];
sum += matrices2[w + 99] * partials2[v + 19];
partials3[u] =  sum;	u++;


sum = matrices2[w + 100] * partials2[v];
sum += matrices2[w + 101] * partials2[v + 1];
sum += matrices2[w + 102] * partials2[v + 2];
sum += matrices2[w + 103] * partials2[v + 3];
sum += matrices2[w + 104] * partials2[v + 4];
sum += matrices2[w + 105] * partials2[v + 5];
sum += matrices2[w + 106] * partials2[v + 6];
sum += matrices2[w + 107] * partials2[v + 7];
sum += matrices2[w + 108] * partials2[v + 8];
sum += matrices2[w + 109] * partials2[v + 9];
sum += matrices2[w + 110] * partials2[v + 10];
sum += matrices2[w + 111] * partials2[v + 11];
sum += matrices2[w + 112] * partials2[v + 12];
sum += matrices2[w + 113] * partials2[v + 13];
sum += matrices2[w + 114] * partials2[v + 14];
sum += matrices2[w + 115] * partials2[v + 15];
sum += matrices2[w + 116] * partials2[v + 16];
sum += matrices2[w + 117] * partials2[v + 17];
sum += matrices2[w + 118] * partials2[v + 18];
sum += matrices2[w + 119] * partials2[v + 19];
partials3[u] =  sum;	u++;


sum = matrices2[w + 120] * partials2[v];
sum += matrices2[w + 121] * partials2[v + 1];
sum += matrices2[w + 122] * partials2[v + 2];
sum += matrices2[w + 123] * partials2[v + 3];
sum += matrices2[w + 124] * partials2[v + 4];
sum += matrices2[w + 125] * partials2[v + 5];
sum += matrices2[w + 126] * partials2[v + 6];
sum += matrices2[w + 127] * partials2[v + 7];
sum += matrices2[w + 128] * partials2[v + 8];
sum += matrices2[w + 129] * partials2[v + 9];
sum += matrices2[w + 130] * partials2[v + 10];
sum += matrices2[w + 131] * partials2[v + 11];
sum += matrices2[w + 132] * partials2[v + 12];
sum += matrices2[w + 133] * partials2[v + 13];
sum += matrices2[w + 134] * partials2[v + 14];
sum += matrices2[w + 135] * partials2[v + 15];
sum += matrices2[w + 136] * partials2[v + 16];
sum += matrices2[w + 137] * partials2[v + 17];
sum += matrices2[w + 138] * partials2[v + 18];
sum += matrices2[w + 139] * partials2[v + 19];
partials3[u] =  sum;	u++;


sum = matrices2[w + 140] * partials2[v];
sum += matrices2[w + 141] * partials2[v + 1];
sum += matrices2[w + 142] * partials2[v + 2];
sum += matrices2[w + 143] * partials2[v + 3];
sum += matrices2[w + 144] * partials2[v + 4];
sum += matrices2[w + 145] * partials2[v + 5];
sum += matrices2[w + 146] * partials2[v + 6];
sum += matrices2[w + 147] * partials2[v + 7];
sum += matrices2[w + 148] * partials2[v + 8];
sum += matrices2[w + 149] * partials2[v + 9];
sum += matrices2[w + 150] * partials2[v + 10];
sum += matrices2[w + 151] * partials2[v + 11];
sum += matrices2[w + 152] * partials2[v + 12];
sum += matrices2[w + 153] * partials2[v + 13];
sum += matrices2[w + 154] * partials2[v + 14];
sum += matrices2[w + 155] * partials2[v + 15];
sum += matrices2[w + 156] * partials2[v + 16];
sum += matrices2[w + 157] * partials2[v + 17];
sum += matrices2[w + 158] * partials2[v + 18];
sum += matrices2[w + 159] * partials2[v + 19];
partials3[u] =  sum;	u++;


sum = matrices2[w + 160] * partials2[v];
sum += matrices2[w + 161] * partials2[v + 1];
sum += matrices2[w + 162] * partials2[v + 2];
sum += matrices2[w + 163] * partials2[v + 3];
sum += matrices2[w + 164] * partials2[v + 4];
sum += matrices2[w + 165] * partials2[v + 5];
sum += matrices2[w + 166] * partials2[v + 6];
sum += matrices2[w + 167] * partials2[v + 7];
sum += matrices2[w + 168] * partials2[v + 8];
sum += matrices2[w + 169] * partials2[v + 9];
sum += matrices2[w + 170] * partials2[v + 10];
sum += matrices2[w + 171] * partials2[v + 11];
sum += matrices2[w + 172] * partials2[v + 12];
sum += matrices2[w + 173] * partials2[v + 13];
sum += matrices2[w + 174] * partials2[v + 14];
sum += matrices2[w + 175] * partials2[v + 15];
sum += matrices2[w + 176] * partials2[v + 16];
sum += matrices2[w + 177] * partials2[v + 17];
sum += matrices2[w + 178] * partials2[v + 18];
sum += matrices2[w + 179] * partials2[v + 19];
partials3[u] =  sum;	u++;


sum = matrices2[w + 180] * partials2[v];
sum += matrices2[w + 181] * partials2[v + 1];
sum += matrices2[w + 182] * partials2[v + 2];
sum += matrices2[w + 183] * partials2[v + 3];
sum += matrices2[w + 184] * partials2[v + 4];
sum += matrices2[w + 185] * partials2[v + 5];
sum += matrices2[w + 186] * partials2[v + 6];
sum += matrices2[w + 187] * partials2[v + 7];
sum += matrices2[w + 188] * partials2[v + 8];
sum += matrices2[w + 189] * partials2[v + 9];
sum += matrices2[w + 190] * partials2[v + 10];
sum += matrices2[w + 191] * partials2[v + 11];
sum += matrices2[w + 192] * partials2[v + 12];
sum += matrices2[w + 193] * partials2[v + 13];
sum += matrices2[w + 194] * partials2[v + 14];
sum += matrices2[w + 195] * partials2[v + 15];
sum += matrices2[w + 196] * partials2[v + 16];
sum += matrices2[w + 197] * partials2[v + 17];
sum += matrices2[w + 198] * partials2[v + 18];
sum += matrices2[w + 199] * partials2[v + 19];
partials3[u] =  sum;	u++;


sum = matrices2[w + 200] * partials2[v];
sum += matrices2[w + 201] * partials2[v + 1];
sum += matrices2[w + 202] * partials2[v + 2];
sum += matrices2[w + 203] * partials2[v + 3];
sum += matrices2[w + 204] * partials2[v + 4];
sum += matrices2[w + 205] * partials2[v + 5];
sum += matrices2[w + 206] * partials2[v + 6];
sum += matrices2[w + 207] * partials2[v + 7];
sum += matrices2[w + 208] * partials2[v + 8];
sum += matrices2[w + 209] * partials2[v + 9];
sum += matrices2[w + 210] * partials2[v + 10];
sum += matrices2[w + 211] * partials2[v + 11];
sum += matrices2[w + 212] * partials2[v + 12];
sum += matrices2[w + 213] * partials2[v + 13];
sum += matrices2[w + 214] * partials2[v + 14];
sum += matrices2[w + 215] * partials2[v + 15];
sum += matrices2[w + 216] * partials2[v + 16];
sum += matrices2[w + 217] * partials2[v + 17];
sum += matrices2[w + 218] * partials2[v + 18];
sum += matrices2[w + 219] * partials2[v + 19];
partials3[u] =  sum;	u++;


sum = matrices2[w + 220] * partials2[v];
sum += matrices2[w + 221] * partials2[v + 1];
sum += matrices2[w + 222] * partials2[v + 2];
sum += matrices2[w + 223] * partials2[v + 3];
sum += matrices2[w + 224] * partials2[v + 4];
sum += matrices2[w + 225] * partials2[v + 5];
sum += matrices2[w + 226] * partials2[v + 6];
sum += matrices2[w + 227] * partials2[v + 7];
sum += matrices2[w + 228] * partials2[v + 8];
sum += matrices2[w + 229] * partials2[v + 9];
sum += matrices2[w + 230] * partials2[v + 10];
sum += matrices2[w + 231] * partials2[v + 11];
sum += matrices2[w + 232] * partials2[v + 12];
sum += matrices2[w + 233] * partials2[v + 13];
sum += matrices2[w + 234] * partials2[v + 14];
sum += matrices2[w + 235] * partials2[v + 15];
sum += matrices2[w + 236] * partials2[v + 16];
sum += matrices2[w + 237] * partials2[v + 17];
sum += matrices2[w + 238] * partials2[v + 18];
sum += matrices2[w + 239] * partials2[v + 19];
partials3[u] =  sum;	u++;


sum = matrices2[w + 240] * partials2[v];
sum += matrices2[w + 241] * partials2[v + 1];
sum += matrices2[w + 242] * partials2[v + 2];
sum += matrices2[w + 243] * partials2[v + 3];
sum += matrices2[w + 244] * partials2[v + 4];
sum += matrices2[w + 245] * partials2[v + 5];
sum += matrices2[w + 246] * partials2[v + 6];
sum += matrices2[w + 247] * partials2[v + 7];
sum += matrices2[w + 248] * partials2[v + 8];
sum += matrices2[w + 249] * partials2[v + 9];
sum += matrices2[w + 250] * partials2[v + 10];
sum += matrices2[w + 251] * partials2[v + 11];
sum += matrices2[w + 252] * partials2[v + 12];
sum += matrices2[w + 253] * partials2[v + 13];
sum += matrices2[w + 254] * partials2[v + 14];
sum += matrices2[w + 255] * partials2[v + 15];
sum += matrices2[w + 256] * partials2[v + 16];
sum += matrices2[w + 257] * partials2[v + 17];
sum += matrices2[w + 258] * partials2[v + 18];
sum += matrices2[w + 259] * partials2[v + 19];
partials3[u] =  sum;	u++;


sum = matrices2[w + 260] * partials2[v];
sum += matrices2[w + 261] * partials2[v + 1];
sum += matrices2[w + 262] * partials2[v + 2];
sum += matrices2[w + 263] * partials2[v + 3];
sum += matrices2[w + 264] * partials2[v + 4];
sum += matrices2[w + 265] * partials2[v + 5];
sum += matrices2[w + 266] * partials2[v + 6];
sum += matrices2[w + 267] * partials2[v + 7];
sum += matrices2[w + 268] * partials2[v + 8];
sum += matrices2[w + 269] * partials2[v + 9];
sum += matrices2[w + 270] * partials2[v + 10];
sum += matrices2[w + 271] * partials2[v + 11];
sum += matrices2[w + 272] * partials2[v + 12];
sum += matrices2[w + 273] * partials2[v + 13];
sum += matrices2[w + 274] * partials2[v + 14];
sum += matrices2[w + 275] * partials2[v + 15];
sum += matrices2[w + 276] * partials2[v + 16];
sum += matrices2[w + 277] * partials2[v + 17];
sum += matrices2[w + 278] * partials2[v + 18];
sum += matrices2[w + 279] * partials2[v + 19];
partials3[u] =  sum;	u++;


sum = matrices2[w + 280] * partials2[v];
sum += matrices2[w + 281] * partials2[v + 1];
sum += matrices2[w + 282] * partials2[v + 2];
sum += matrices2[w + 283] * partials2[v + 3];
sum += matrices2[w + 284] * partials2[v + 4];
sum += matrices2[w + 285] * partials2[v + 5];
sum += matrices2[w + 286] * partials2[v + 6];
sum += matrices2[w + 287] * partials2[v + 7];
sum += matrices2[w + 288] * partials2[v + 8];
sum += matrices2[w + 289] * partials2[v + 9];
sum += matrices2[w + 290] * partials2[v + 10];
sum += matrices2[w + 291] * partials2[v + 11];
sum += matrices2[w + 292] * partials2[v + 12];
sum += matrices2[w + 293] * partials2[v + 13];
sum += matrices2[w + 294] * partials2[v + 14];
sum += matrices2[w + 295] * partials2[v + 15];
sum += matrices2[w + 296] * partials2[v + 16];
sum += matrices2[w + 297] * partials2[v + 17];
sum += matrices2[w + 298] * partials2[v + 18];
sum += matrices2[w + 299] * partials2[v + 19];
partials3[u] =  sum;	u++;


sum = matrices2[w + 300] * partials2[v];
sum += matrices2[w + 301] * partials2[v + 1];
sum += matrices2[w + 302] * partials2[v + 2];
sum += matrices2[w + 303] * partials2[v + 3];
sum += matrices2[w + 304] * partials2[v + 4];
sum += matrices2[w + 305] * partials2[v + 5];
sum += matrices2[w + 306] * partials2[v + 6];
sum += matrices2[w + 307] * partials2[v + 7];
sum += matrices2[w + 308] * partials2[v + 8];
sum += matrices2[w + 309] * partials2[v + 9];
sum += matrices2[w + 310] * partials2[v + 10];
sum += matrices2[w + 311] * partials2[v + 11];
sum += matrices2[w + 312] * partials2[v + 12];
sum += matrices2[w + 313] * partials2[v + 13];
sum += matrices2[w + 314] * partials2[v + 14];
sum += matrices2[w + 315] * partials2[v + 15];
sum += matrices2[w + 316] * partials2[v + 16];
sum += matrices2[w + 317] * partials2[v + 17];
sum += matrices2[w + 318] * partials2[v + 18];
sum += matrices2[w + 319] * partials2[v + 19];
partials3[u] =  sum;	u++;


sum = matrices2[w + 320] * partials2[v];
sum += matrices2[w + 321] * partials2[v + 1];
sum += matrices2[w + 322] * partials2[v + 2];
sum += matrices2[w + 323] * partials2[v + 3];
sum += matrices2[w + 324] * partials2[v + 4];
sum += matrices2[w + 325] * partials2[v + 5];
sum += matrices2[w + 326] * partials2[v + 6];
sum += matrices2[w + 327] * partials2[v + 7];
sum += matrices2[w + 328] * partials2[v + 8];
sum += matrices2[w + 329] * partials2[v + 9];
sum += matrices2[w + 330] * partials2[v + 10];
sum += matrices2[w + 331] * partials2[v + 11];
sum += matrices2[w + 332] * partials2[v + 12];
sum += matrices2[w + 333] * partials2[v + 13];
sum += matrices2[w + 334] * partials2[v + 14];
sum += matrices2[w + 335] * partials2[v + 15];
sum += matrices2[w + 336] * partials2[v + 16];
sum += matrices2[w + 337] * partials2[v + 17];
sum += matrices2[w + 338] * partials2[v + 18];
sum += matrices2[w + 339] * partials2[v + 19];
partials3[u] =  sum;	u++;


sum = matrices2[w + 340] * partials2[v];
sum += matrices2[w + 341] * partials2[v + 1];
sum += matrices2[w + 342] * partials2[v + 2];
sum += matrices2[w + 343] * partials2[v + 3];
sum += matrices2[w + 344] * partials2[v + 4];
sum += matrices2[w + 345] * partials2[v + 5];
sum += matrices2[w + 346] * partials2[v + 6];
sum += matrices2[w + 347] * partials2[v + 7];
sum += matrices2[w + 348] * partials2[v + 8];
sum += matrices2[w + 349] * partials2[v + 9];
sum += matrices2[w + 350] * partials2[v + 10];
sum += matrices2[w + 351] * partials2[v + 11];
sum += matrices2[w + 352] * partials2[v + 12];
sum += matrices2[w + 353] * partials2[v + 13];
sum += matrices2[w + 354] * partials2[v + 14];
sum += matrices2[w + 355] * partials2[v + 15];
sum += matrices2[w + 356] * partials2[v + 16];
sum += matrices2[w + 357] * partials2[v + 17];
sum += matrices2[w + 358] * partials2[v + 18];
sum += matrices2[w + 359] * partials2[v + 19];
partials3[u] =  sum;	u++;


sum = matrices2[w + 360] * partials2[v];
sum += matrices2[w + 361] * partials2[v + 1];
sum += matrices2[w + 362] * partials2[v + 2];
sum += matrices2[w + 363] * partials2[v + 3];
sum += matrices2[w + 364] * partials2[v + 4];
sum += matrices2[w + 365] * partials2[v + 5];
sum += matrices2[w + 366] * partials2[v + 6];
sum += matrices2[w + 367] * partials2[v + 7];
sum += matrices2[w + 368] * partials2[v + 8];
sum += matrices2[w + 369] * partials2[v + 9];
sum += matrices2[w + 370] * partials2[v + 10];
sum += matrices2[w + 371] * partials2[v + 11];
sum += matrices2[w + 372] * partials2[v + 12];
sum += matrices2[w + 373] * partials2[v + 13];
sum += matrices2[w + 374] * partials2[v + 14];
sum += matrices2[w + 375] * partials2[v + 15];
sum += matrices2[w + 376] * partials2[v + 16];
sum += matrices2[w + 377] * partials2[v + 17];
sum += matrices2[w + 378] * partials2[v + 18];
sum += matrices2[w + 379] * partials2[v + 19];
partials3[u] =  sum;	u++;


sum = matrices2[w + 380] * partials2[v];
sum += matrices2[w + 381] * partials2[v + 1];
sum += matrices2[w + 382] * partials2[v + 2];
sum += matrices2[w + 383] * partials2[v + 3];
sum += matrices2[w + 384] * partials2[v + 4];
sum += matrices2[w + 385] * partials2[v + 5];
sum += matrices2[w + 386] * partials2[v + 6];
sum += matrices2[w + 387] * partials2[v + 7];
sum += matrices2[w + 388] * partials2[v + 8];
sum += matrices2[w + 389] * partials2[v + 9];
sum += matrices2[w + 390] * partials2[v + 10];
sum += matrices2[w + 391] * partials2[v + 11];
sum += matrices2[w + 392] * partials2[v + 12];
sum += matrices2[w + 393] * partials2[v + 13];
sum += matrices2[w + 394] * partials2[v + 14];
sum += matrices2[w + 395] * partials2[v + 15];
sum += matrices2[w + 396] * partials2[v + 16];
sum += matrices2[w + 397] * partials2[v + 17];
sum += matrices2[w + 398] * partials2[v + 18];
sum += matrices2[w + 399] * partials2[v + 19];
partials3[u] =  sum;	u++;

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
sum1 += matrices1[w + 8] * partials1[v + 8];
sum2 += matrices2[w + 8] * partials2[v + 8];
sum1 += matrices1[w + 9] * partials1[v + 9];
sum2 += matrices2[w + 9] * partials2[v + 9];
sum1 += matrices1[w + 10] * partials1[v + 10];
sum2 += matrices2[w + 10] * partials2[v + 10];
sum1 += matrices1[w + 11] * partials1[v + 11];
sum2 += matrices2[w + 11] * partials2[v + 11];
sum1 += matrices1[w + 12] * partials1[v + 12];
sum2 += matrices2[w + 12] * partials2[v + 12];
sum1 += matrices1[w + 13] * partials1[v + 13];
sum2 += matrices2[w + 13] * partials2[v + 13];
sum1 += matrices1[w + 14] * partials1[v + 14];
sum2 += matrices2[w + 14] * partials2[v + 14];
sum1 += matrices1[w + 15] * partials1[v + 15];
sum2 += matrices2[w + 15] * partials2[v + 15];
sum1 += matrices1[w + 16] * partials1[v + 16];
sum2 += matrices2[w + 16] * partials2[v + 16];
sum1 += matrices1[w + 17] * partials1[v + 17];
sum2 += matrices2[w + 17] * partials2[v + 17];
sum1 += matrices1[w + 18] * partials1[v + 18];
sum2 += matrices2[w + 18] * partials2[v + 18];
sum1 += matrices1[w + 19] * partials1[v + 19];
sum2 += matrices2[w + 19] * partials2[v + 19];
partials3[u] = sum1 * sum2; u++;

sum1 = matrices1[w + 20] * partials1[v];
sum2 = matrices2[w + 20] * partials2[v];
sum1 += matrices1[w + 21] * partials1[v + 1];
sum2 += matrices2[w + 21] * partials2[v + 1];
sum1 += matrices1[w + 22] * partials1[v + 2];
sum2 += matrices2[w + 22] * partials2[v + 2];
sum1 += matrices1[w + 23] * partials1[v + 3];
sum2 += matrices2[w + 23] * partials2[v + 3];
sum1 += matrices1[w + 24] * partials1[v + 4];
sum2 += matrices2[w + 24] * partials2[v + 4];
sum1 += matrices1[w + 25] * partials1[v + 5];
sum2 += matrices2[w + 25] * partials2[v + 5];
sum1 += matrices1[w + 26] * partials1[v + 6];
sum2 += matrices2[w + 26] * partials2[v + 6];
sum1 += matrices1[w + 27] * partials1[v + 7];
sum2 += matrices2[w + 27] * partials2[v + 7];
sum1 += matrices1[w + 28] * partials1[v + 8];
sum2 += matrices2[w + 28] * partials2[v + 8];
sum1 += matrices1[w + 29] * partials1[v + 9];
sum2 += matrices2[w + 29] * partials2[v + 9];
sum1 += matrices1[w + 30] * partials1[v + 10];
sum2 += matrices2[w + 30] * partials2[v + 10];
sum1 += matrices1[w + 31] * partials1[v + 11];
sum2 += matrices2[w + 31] * partials2[v + 11];
sum1 += matrices1[w + 32] * partials1[v + 12];
sum2 += matrices2[w + 32] * partials2[v + 12];
sum1 += matrices1[w + 33] * partials1[v + 13];
sum2 += matrices2[w + 33] * partials2[v + 13];
sum1 += matrices1[w + 34] * partials1[v + 14];
sum2 += matrices2[w + 34] * partials2[v + 14];
sum1 += matrices1[w + 35] * partials1[v + 15];
sum2 += matrices2[w + 35] * partials2[v + 15];
sum1 += matrices1[w + 36] * partials1[v + 16];
sum2 += matrices2[w + 36] * partials2[v + 16];
sum1 += matrices1[w + 37] * partials1[v + 17];
sum2 += matrices2[w + 37] * partials2[v + 17];
sum1 += matrices1[w + 38] * partials1[v + 18];
sum2 += matrices2[w + 38] * partials2[v + 18];
sum1 += matrices1[w + 39] * partials1[v + 19];
sum2 += matrices2[w + 39] * partials2[v + 19];
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
sum1 += matrices1[w + 48] * partials1[v + 8];
sum2 += matrices2[w + 48] * partials2[v + 8];
sum1 += matrices1[w + 49] * partials1[v + 9];
sum2 += matrices2[w + 49] * partials2[v + 9];
sum1 += matrices1[w + 50] * partials1[v + 10];
sum2 += matrices2[w + 50] * partials2[v + 10];
sum1 += matrices1[w + 51] * partials1[v + 11];
sum2 += matrices2[w + 51] * partials2[v + 11];
sum1 += matrices1[w + 52] * partials1[v + 12];
sum2 += matrices2[w + 52] * partials2[v + 12];
sum1 += matrices1[w + 53] * partials1[v + 13];
sum2 += matrices2[w + 53] * partials2[v + 13];
sum1 += matrices1[w + 54] * partials1[v + 14];
sum2 += matrices2[w + 54] * partials2[v + 14];
sum1 += matrices1[w + 55] * partials1[v + 15];
sum2 += matrices2[w + 55] * partials2[v + 15];
sum1 += matrices1[w + 56] * partials1[v + 16];
sum2 += matrices2[w + 56] * partials2[v + 16];
sum1 += matrices1[w + 57] * partials1[v + 17];
sum2 += matrices2[w + 57] * partials2[v + 17];
sum1 += matrices1[w + 58] * partials1[v + 18];
sum2 += matrices2[w + 58] * partials2[v + 18];
sum1 += matrices1[w + 59] * partials1[v + 19];
sum2 += matrices2[w + 59] * partials2[v + 19];
partials3[u] = sum1 * sum2; u++;

sum1 = matrices1[w + 60] * partials1[v];
sum2 = matrices2[w + 60] * partials2[v];
sum1 += matrices1[w + 61] * partials1[v + 1];
sum2 += matrices2[w + 61] * partials2[v + 1];
sum1 += matrices1[w + 62] * partials1[v + 2];
sum2 += matrices2[w + 62] * partials2[v + 2];
sum1 += matrices1[w + 63] * partials1[v + 3];
sum2 += matrices2[w + 63] * partials2[v + 3];
sum1 += matrices1[w + 64] * partials1[v + 4];
sum2 += matrices2[w + 64] * partials2[v + 4];
sum1 += matrices1[w + 65] * partials1[v + 5];
sum2 += matrices2[w + 65] * partials2[v + 5];
sum1 += matrices1[w + 66] * partials1[v + 6];
sum2 += matrices2[w + 66] * partials2[v + 6];
sum1 += matrices1[w + 67] * partials1[v + 7];
sum2 += matrices2[w + 67] * partials2[v + 7];
sum1 += matrices1[w + 68] * partials1[v + 8];
sum2 += matrices2[w + 68] * partials2[v + 8];
sum1 += matrices1[w + 69] * partials1[v + 9];
sum2 += matrices2[w + 69] * partials2[v + 9];
sum1 += matrices1[w + 70] * partials1[v + 10];
sum2 += matrices2[w + 70] * partials2[v + 10];
sum1 += matrices1[w + 71] * partials1[v + 11];
sum2 += matrices2[w + 71] * partials2[v + 11];
sum1 += matrices1[w + 72] * partials1[v + 12];
sum2 += matrices2[w + 72] * partials2[v + 12];
sum1 += matrices1[w + 73] * partials1[v + 13];
sum2 += matrices2[w + 73] * partials2[v + 13];
sum1 += matrices1[w + 74] * partials1[v + 14];
sum2 += matrices2[w + 74] * partials2[v + 14];
sum1 += matrices1[w + 75] * partials1[v + 15];
sum2 += matrices2[w + 75] * partials2[v + 15];
sum1 += matrices1[w + 76] * partials1[v + 16];
sum2 += matrices2[w + 76] * partials2[v + 16];
sum1 += matrices1[w + 77] * partials1[v + 17];
sum2 += matrices2[w + 77] * partials2[v + 17];
sum1 += matrices1[w + 78] * partials1[v + 18];
sum2 += matrices2[w + 78] * partials2[v + 18];
sum1 += matrices1[w + 79] * partials1[v + 19];
sum2 += matrices2[w + 79] * partials2[v + 19];
partials3[u] = sum1 * sum2; u++;

sum1 = matrices1[w + 80] * partials1[v];
sum2 = matrices2[w + 80] * partials2[v];
sum1 += matrices1[w + 81] * partials1[v + 1];
sum2 += matrices2[w + 81] * partials2[v + 1];
sum1 += matrices1[w + 82] * partials1[v + 2];
sum2 += matrices2[w + 82] * partials2[v + 2];
sum1 += matrices1[w + 83] * partials1[v + 3];
sum2 += matrices2[w + 83] * partials2[v + 3];
sum1 += matrices1[w + 84] * partials1[v + 4];
sum2 += matrices2[w + 84] * partials2[v + 4];
sum1 += matrices1[w + 85] * partials1[v + 5];
sum2 += matrices2[w + 85] * partials2[v + 5];
sum1 += matrices1[w + 86] * partials1[v + 6];
sum2 += matrices2[w + 86] * partials2[v + 6];
sum1 += matrices1[w + 87] * partials1[v + 7];
sum2 += matrices2[w + 87] * partials2[v + 7];
sum1 += matrices1[w + 88] * partials1[v + 8];
sum2 += matrices2[w + 88] * partials2[v + 8];
sum1 += matrices1[w + 89] * partials1[v + 9];
sum2 += matrices2[w + 89] * partials2[v + 9];
sum1 += matrices1[w + 90] * partials1[v + 10];
sum2 += matrices2[w + 90] * partials2[v + 10];
sum1 += matrices1[w + 91] * partials1[v + 11];
sum2 += matrices2[w + 91] * partials2[v + 11];
sum1 += matrices1[w + 92] * partials1[v + 12];
sum2 += matrices2[w + 92] * partials2[v + 12];
sum1 += matrices1[w + 93] * partials1[v + 13];
sum2 += matrices2[w + 93] * partials2[v + 13];
sum1 += matrices1[w + 94] * partials1[v + 14];
sum2 += matrices2[w + 94] * partials2[v + 14];
sum1 += matrices1[w + 95] * partials1[v + 15];
sum2 += matrices2[w + 95] * partials2[v + 15];
sum1 += matrices1[w + 96] * partials1[v + 16];
sum2 += matrices2[w + 96] * partials2[v + 16];
sum1 += matrices1[w + 97] * partials1[v + 17];
sum2 += matrices2[w + 97] * partials2[v + 17];
sum1 += matrices1[w + 98] * partials1[v + 18];
sum2 += matrices2[w + 98] * partials2[v + 18];
sum1 += matrices1[w + 99] * partials1[v + 19];
sum2 += matrices2[w + 99] * partials2[v + 19];
partials3[u] = sum1 * sum2; u++;

sum1 = matrices1[w + 100] * partials1[v];
sum2 = matrices2[w + 100] * partials2[v];
sum1 += matrices1[w + 101] * partials1[v + 1];
sum2 += matrices2[w + 101] * partials2[v + 1];
sum1 += matrices1[w + 102] * partials1[v + 2];
sum2 += matrices2[w + 102] * partials2[v + 2];
sum1 += matrices1[w + 103] * partials1[v + 3];
sum2 += matrices2[w + 103] * partials2[v + 3];
sum1 += matrices1[w + 104] * partials1[v + 4];
sum2 += matrices2[w + 104] * partials2[v + 4];
sum1 += matrices1[w + 105] * partials1[v + 5];
sum2 += matrices2[w + 105] * partials2[v + 5];
sum1 += matrices1[w + 106] * partials1[v + 6];
sum2 += matrices2[w + 106] * partials2[v + 6];
sum1 += matrices1[w + 107] * partials1[v + 7];
sum2 += matrices2[w + 107] * partials2[v + 7];
sum1 += matrices1[w + 108] * partials1[v + 8];
sum2 += matrices2[w + 108] * partials2[v + 8];
sum1 += matrices1[w + 109] * partials1[v + 9];
sum2 += matrices2[w + 109] * partials2[v + 9];
sum1 += matrices1[w + 110] * partials1[v + 10];
sum2 += matrices2[w + 110] * partials2[v + 10];
sum1 += matrices1[w + 111] * partials1[v + 11];
sum2 += matrices2[w + 111] * partials2[v + 11];
sum1 += matrices1[w + 112] * partials1[v + 12];
sum2 += matrices2[w + 112] * partials2[v + 12];
sum1 += matrices1[w + 113] * partials1[v + 13];
sum2 += matrices2[w + 113] * partials2[v + 13];
sum1 += matrices1[w + 114] * partials1[v + 14];
sum2 += matrices2[w + 114] * partials2[v + 14];
sum1 += matrices1[w + 115] * partials1[v + 15];
sum2 += matrices2[w + 115] * partials2[v + 15];
sum1 += matrices1[w + 116] * partials1[v + 16];
sum2 += matrices2[w + 116] * partials2[v + 16];
sum1 += matrices1[w + 117] * partials1[v + 17];
sum2 += matrices2[w + 117] * partials2[v + 17];
sum1 += matrices1[w + 118] * partials1[v + 18];
sum2 += matrices2[w + 118] * partials2[v + 18];
sum1 += matrices1[w + 119] * partials1[v + 19];
sum2 += matrices2[w + 119] * partials2[v + 19];
partials3[u] = sum1 * sum2; u++;

sum1 = matrices1[w + 120] * partials1[v];
sum2 = matrices2[w + 120] * partials2[v];
sum1 += matrices1[w + 121] * partials1[v + 1];
sum2 += matrices2[w + 121] * partials2[v + 1];
sum1 += matrices1[w + 122] * partials1[v + 2];
sum2 += matrices2[w + 122] * partials2[v + 2];
sum1 += matrices1[w + 123] * partials1[v + 3];
sum2 += matrices2[w + 123] * partials2[v + 3];
sum1 += matrices1[w + 124] * partials1[v + 4];
sum2 += matrices2[w + 124] * partials2[v + 4];
sum1 += matrices1[w + 125] * partials1[v + 5];
sum2 += matrices2[w + 125] * partials2[v + 5];
sum1 += matrices1[w + 126] * partials1[v + 6];
sum2 += matrices2[w + 126] * partials2[v + 6];
sum1 += matrices1[w + 127] * partials1[v + 7];
sum2 += matrices2[w + 127] * partials2[v + 7];
sum1 += matrices1[w + 128] * partials1[v + 8];
sum2 += matrices2[w + 128] * partials2[v + 8];
sum1 += matrices1[w + 129] * partials1[v + 9];
sum2 += matrices2[w + 129] * partials2[v + 9];
sum1 += matrices1[w + 130] * partials1[v + 10];
sum2 += matrices2[w + 130] * partials2[v + 10];
sum1 += matrices1[w + 131] * partials1[v + 11];
sum2 += matrices2[w + 131] * partials2[v + 11];
sum1 += matrices1[w + 132] * partials1[v + 12];
sum2 += matrices2[w + 132] * partials2[v + 12];
sum1 += matrices1[w + 133] * partials1[v + 13];
sum2 += matrices2[w + 133] * partials2[v + 13];
sum1 += matrices1[w + 134] * partials1[v + 14];
sum2 += matrices2[w + 134] * partials2[v + 14];
sum1 += matrices1[w + 135] * partials1[v + 15];
sum2 += matrices2[w + 135] * partials2[v + 15];
sum1 += matrices1[w + 136] * partials1[v + 16];
sum2 += matrices2[w + 136] * partials2[v + 16];
sum1 += matrices1[w + 137] * partials1[v + 17];
sum2 += matrices2[w + 137] * partials2[v + 17];
sum1 += matrices1[w + 138] * partials1[v + 18];
sum2 += matrices2[w + 138] * partials2[v + 18];
sum1 += matrices1[w + 139] * partials1[v + 19];
sum2 += matrices2[w + 139] * partials2[v + 19];
partials3[u] = sum1 * sum2; u++;

sum1 = matrices1[w + 140] * partials1[v];
sum2 = matrices2[w + 140] * partials2[v];
sum1 += matrices1[w + 141] * partials1[v + 1];
sum2 += matrices2[w + 141] * partials2[v + 1];
sum1 += matrices1[w + 142] * partials1[v + 2];
sum2 += matrices2[w + 142] * partials2[v + 2];
sum1 += matrices1[w + 143] * partials1[v + 3];
sum2 += matrices2[w + 143] * partials2[v + 3];
sum1 += matrices1[w + 144] * partials1[v + 4];
sum2 += matrices2[w + 144] * partials2[v + 4];
sum1 += matrices1[w + 145] * partials1[v + 5];
sum2 += matrices2[w + 145] * partials2[v + 5];
sum1 += matrices1[w + 146] * partials1[v + 6];
sum2 += matrices2[w + 146] * partials2[v + 6];
sum1 += matrices1[w + 147] * partials1[v + 7];
sum2 += matrices2[w + 147] * partials2[v + 7];
sum1 += matrices1[w + 148] * partials1[v + 8];
sum2 += matrices2[w + 148] * partials2[v + 8];
sum1 += matrices1[w + 149] * partials1[v + 9];
sum2 += matrices2[w + 149] * partials2[v + 9];
sum1 += matrices1[w + 150] * partials1[v + 10];
sum2 += matrices2[w + 150] * partials2[v + 10];
sum1 += matrices1[w + 151] * partials1[v + 11];
sum2 += matrices2[w + 151] * partials2[v + 11];
sum1 += matrices1[w + 152] * partials1[v + 12];
sum2 += matrices2[w + 152] * partials2[v + 12];
sum1 += matrices1[w + 153] * partials1[v + 13];
sum2 += matrices2[w + 153] * partials2[v + 13];
sum1 += matrices1[w + 154] * partials1[v + 14];
sum2 += matrices2[w + 154] * partials2[v + 14];
sum1 += matrices1[w + 155] * partials1[v + 15];
sum2 += matrices2[w + 155] * partials2[v + 15];
sum1 += matrices1[w + 156] * partials1[v + 16];
sum2 += matrices2[w + 156] * partials2[v + 16];
sum1 += matrices1[w + 157] * partials1[v + 17];
sum2 += matrices2[w + 157] * partials2[v + 17];
sum1 += matrices1[w + 158] * partials1[v + 18];
sum2 += matrices2[w + 158] * partials2[v + 18];
sum1 += matrices1[w + 159] * partials1[v + 19];
sum2 += matrices2[w + 159] * partials2[v + 19];
partials3[u] = sum1 * sum2; u++;

sum1 = matrices1[w + 160] * partials1[v];
sum2 = matrices2[w + 160] * partials2[v];
sum1 += matrices1[w + 161] * partials1[v + 1];
sum2 += matrices2[w + 161] * partials2[v + 1];
sum1 += matrices1[w + 162] * partials1[v + 2];
sum2 += matrices2[w + 162] * partials2[v + 2];
sum1 += matrices1[w + 163] * partials1[v + 3];
sum2 += matrices2[w + 163] * partials2[v + 3];
sum1 += matrices1[w + 164] * partials1[v + 4];
sum2 += matrices2[w + 164] * partials2[v + 4];
sum1 += matrices1[w + 165] * partials1[v + 5];
sum2 += matrices2[w + 165] * partials2[v + 5];
sum1 += matrices1[w + 166] * partials1[v + 6];
sum2 += matrices2[w + 166] * partials2[v + 6];
sum1 += matrices1[w + 167] * partials1[v + 7];
sum2 += matrices2[w + 167] * partials2[v + 7];
sum1 += matrices1[w + 168] * partials1[v + 8];
sum2 += matrices2[w + 168] * partials2[v + 8];
sum1 += matrices1[w + 169] * partials1[v + 9];
sum2 += matrices2[w + 169] * partials2[v + 9];
sum1 += matrices1[w + 170] * partials1[v + 10];
sum2 += matrices2[w + 170] * partials2[v + 10];
sum1 += matrices1[w + 171] * partials1[v + 11];
sum2 += matrices2[w + 171] * partials2[v + 11];
sum1 += matrices1[w + 172] * partials1[v + 12];
sum2 += matrices2[w + 172] * partials2[v + 12];
sum1 += matrices1[w + 173] * partials1[v + 13];
sum2 += matrices2[w + 173] * partials2[v + 13];
sum1 += matrices1[w + 174] * partials1[v + 14];
sum2 += matrices2[w + 174] * partials2[v + 14];
sum1 += matrices1[w + 175] * partials1[v + 15];
sum2 += matrices2[w + 175] * partials2[v + 15];
sum1 += matrices1[w + 176] * partials1[v + 16];
sum2 += matrices2[w + 176] * partials2[v + 16];
sum1 += matrices1[w + 177] * partials1[v + 17];
sum2 += matrices2[w + 177] * partials2[v + 17];
sum1 += matrices1[w + 178] * partials1[v + 18];
sum2 += matrices2[w + 178] * partials2[v + 18];
sum1 += matrices1[w + 179] * partials1[v + 19];
sum2 += matrices2[w + 179] * partials2[v + 19];
partials3[u] = sum1 * sum2; u++;

sum1 = matrices1[w + 180] * partials1[v];
sum2 = matrices2[w + 180] * partials2[v];
sum1 += matrices1[w + 181] * partials1[v + 1];
sum2 += matrices2[w + 181] * partials2[v + 1];
sum1 += matrices1[w + 182] * partials1[v + 2];
sum2 += matrices2[w + 182] * partials2[v + 2];
sum1 += matrices1[w + 183] * partials1[v + 3];
sum2 += matrices2[w + 183] * partials2[v + 3];
sum1 += matrices1[w + 184] * partials1[v + 4];
sum2 += matrices2[w + 184] * partials2[v + 4];
sum1 += matrices1[w + 185] * partials1[v + 5];
sum2 += matrices2[w + 185] * partials2[v + 5];
sum1 += matrices1[w + 186] * partials1[v + 6];
sum2 += matrices2[w + 186] * partials2[v + 6];
sum1 += matrices1[w + 187] * partials1[v + 7];
sum2 += matrices2[w + 187] * partials2[v + 7];
sum1 += matrices1[w + 188] * partials1[v + 8];
sum2 += matrices2[w + 188] * partials2[v + 8];
sum1 += matrices1[w + 189] * partials1[v + 9];
sum2 += matrices2[w + 189] * partials2[v + 9];
sum1 += matrices1[w + 190] * partials1[v + 10];
sum2 += matrices2[w + 190] * partials2[v + 10];
sum1 += matrices1[w + 191] * partials1[v + 11];
sum2 += matrices2[w + 191] * partials2[v + 11];
sum1 += matrices1[w + 192] * partials1[v + 12];
sum2 += matrices2[w + 192] * partials2[v + 12];
sum1 += matrices1[w + 193] * partials1[v + 13];
sum2 += matrices2[w + 193] * partials2[v + 13];
sum1 += matrices1[w + 194] * partials1[v + 14];
sum2 += matrices2[w + 194] * partials2[v + 14];
sum1 += matrices1[w + 195] * partials1[v + 15];
sum2 += matrices2[w + 195] * partials2[v + 15];
sum1 += matrices1[w + 196] * partials1[v + 16];
sum2 += matrices2[w + 196] * partials2[v + 16];
sum1 += matrices1[w + 197] * partials1[v + 17];
sum2 += matrices2[w + 197] * partials2[v + 17];
sum1 += matrices1[w + 198] * partials1[v + 18];
sum2 += matrices2[w + 198] * partials2[v + 18];
sum1 += matrices1[w + 199] * partials1[v + 19];
sum2 += matrices2[w + 199] * partials2[v + 19];
partials3[u] = sum1 * sum2; u++;

sum1 = matrices1[w + 200] * partials1[v];
sum2 = matrices2[w + 200] * partials2[v];
sum1 += matrices1[w + 201] * partials1[v + 1];
sum2 += matrices2[w + 201] * partials2[v + 1];
sum1 += matrices1[w + 202] * partials1[v + 2];
sum2 += matrices2[w + 202] * partials2[v + 2];
sum1 += matrices1[w + 203] * partials1[v + 3];
sum2 += matrices2[w + 203] * partials2[v + 3];
sum1 += matrices1[w + 204] * partials1[v + 4];
sum2 += matrices2[w + 204] * partials2[v + 4];
sum1 += matrices1[w + 205] * partials1[v + 5];
sum2 += matrices2[w + 205] * partials2[v + 5];
sum1 += matrices1[w + 206] * partials1[v + 6];
sum2 += matrices2[w + 206] * partials2[v + 6];
sum1 += matrices1[w + 207] * partials1[v + 7];
sum2 += matrices2[w + 207] * partials2[v + 7];
sum1 += matrices1[w + 208] * partials1[v + 8];
sum2 += matrices2[w + 208] * partials2[v + 8];
sum1 += matrices1[w + 209] * partials1[v + 9];
sum2 += matrices2[w + 209] * partials2[v + 9];
sum1 += matrices1[w + 210] * partials1[v + 10];
sum2 += matrices2[w + 210] * partials2[v + 10];
sum1 += matrices1[w + 211] * partials1[v + 11];
sum2 += matrices2[w + 211] * partials2[v + 11];
sum1 += matrices1[w + 212] * partials1[v + 12];
sum2 += matrices2[w + 212] * partials2[v + 12];
sum1 += matrices1[w + 213] * partials1[v + 13];
sum2 += matrices2[w + 213] * partials2[v + 13];
sum1 += matrices1[w + 214] * partials1[v + 14];
sum2 += matrices2[w + 214] * partials2[v + 14];
sum1 += matrices1[w + 215] * partials1[v + 15];
sum2 += matrices2[w + 215] * partials2[v + 15];
sum1 += matrices1[w + 216] * partials1[v + 16];
sum2 += matrices2[w + 216] * partials2[v + 16];
sum1 += matrices1[w + 217] * partials1[v + 17];
sum2 += matrices2[w + 217] * partials2[v + 17];
sum1 += matrices1[w + 218] * partials1[v + 18];
sum2 += matrices2[w + 218] * partials2[v + 18];
sum1 += matrices1[w + 219] * partials1[v + 19];
sum2 += matrices2[w + 219] * partials2[v + 19];
partials3[u] = sum1 * sum2; u++;

sum1 = matrices1[w + 220] * partials1[v];
sum2 = matrices2[w + 220] * partials2[v];
sum1 += matrices1[w + 221] * partials1[v + 1];
sum2 += matrices2[w + 221] * partials2[v + 1];
sum1 += matrices1[w + 222] * partials1[v + 2];
sum2 += matrices2[w + 222] * partials2[v + 2];
sum1 += matrices1[w + 223] * partials1[v + 3];
sum2 += matrices2[w + 223] * partials2[v + 3];
sum1 += matrices1[w + 224] * partials1[v + 4];
sum2 += matrices2[w + 224] * partials2[v + 4];
sum1 += matrices1[w + 225] * partials1[v + 5];
sum2 += matrices2[w + 225] * partials2[v + 5];
sum1 += matrices1[w + 226] * partials1[v + 6];
sum2 += matrices2[w + 226] * partials2[v + 6];
sum1 += matrices1[w + 227] * partials1[v + 7];
sum2 += matrices2[w + 227] * partials2[v + 7];
sum1 += matrices1[w + 228] * partials1[v + 8];
sum2 += matrices2[w + 228] * partials2[v + 8];
sum1 += matrices1[w + 229] * partials1[v + 9];
sum2 += matrices2[w + 229] * partials2[v + 9];
sum1 += matrices1[w + 230] * partials1[v + 10];
sum2 += matrices2[w + 230] * partials2[v + 10];
sum1 += matrices1[w + 231] * partials1[v + 11];
sum2 += matrices2[w + 231] * partials2[v + 11];
sum1 += matrices1[w + 232] * partials1[v + 12];
sum2 += matrices2[w + 232] * partials2[v + 12];
sum1 += matrices1[w + 233] * partials1[v + 13];
sum2 += matrices2[w + 233] * partials2[v + 13];
sum1 += matrices1[w + 234] * partials1[v + 14];
sum2 += matrices2[w + 234] * partials2[v + 14];
sum1 += matrices1[w + 235] * partials1[v + 15];
sum2 += matrices2[w + 235] * partials2[v + 15];
sum1 += matrices1[w + 236] * partials1[v + 16];
sum2 += matrices2[w + 236] * partials2[v + 16];
sum1 += matrices1[w + 237] * partials1[v + 17];
sum2 += matrices2[w + 237] * partials2[v + 17];
sum1 += matrices1[w + 238] * partials1[v + 18];
sum2 += matrices2[w + 238] * partials2[v + 18];
sum1 += matrices1[w + 239] * partials1[v + 19];
sum2 += matrices2[w + 239] * partials2[v + 19];
partials3[u] = sum1 * sum2; u++;

sum1 = matrices1[w + 240] * partials1[v];
sum2 = matrices2[w + 240] * partials2[v];
sum1 += matrices1[w + 241] * partials1[v + 1];
sum2 += matrices2[w + 241] * partials2[v + 1];
sum1 += matrices1[w + 242] * partials1[v + 2];
sum2 += matrices2[w + 242] * partials2[v + 2];
sum1 += matrices1[w + 243] * partials1[v + 3];
sum2 += matrices2[w + 243] * partials2[v + 3];
sum1 += matrices1[w + 244] * partials1[v + 4];
sum2 += matrices2[w + 244] * partials2[v + 4];
sum1 += matrices1[w + 245] * partials1[v + 5];
sum2 += matrices2[w + 245] * partials2[v + 5];
sum1 += matrices1[w + 246] * partials1[v + 6];
sum2 += matrices2[w + 246] * partials2[v + 6];
sum1 += matrices1[w + 247] * partials1[v + 7];
sum2 += matrices2[w + 247] * partials2[v + 7];
sum1 += matrices1[w + 248] * partials1[v + 8];
sum2 += matrices2[w + 248] * partials2[v + 8];
sum1 += matrices1[w + 249] * partials1[v + 9];
sum2 += matrices2[w + 249] * partials2[v + 9];
sum1 += matrices1[w + 250] * partials1[v + 10];
sum2 += matrices2[w + 250] * partials2[v + 10];
sum1 += matrices1[w + 251] * partials1[v + 11];
sum2 += matrices2[w + 251] * partials2[v + 11];
sum1 += matrices1[w + 252] * partials1[v + 12];
sum2 += matrices2[w + 252] * partials2[v + 12];
sum1 += matrices1[w + 253] * partials1[v + 13];
sum2 += matrices2[w + 253] * partials2[v + 13];
sum1 += matrices1[w + 254] * partials1[v + 14];
sum2 += matrices2[w + 254] * partials2[v + 14];
sum1 += matrices1[w + 255] * partials1[v + 15];
sum2 += matrices2[w + 255] * partials2[v + 15];
sum1 += matrices1[w + 256] * partials1[v + 16];
sum2 += matrices2[w + 256] * partials2[v + 16];
sum1 += matrices1[w + 257] * partials1[v + 17];
sum2 += matrices2[w + 257] * partials2[v + 17];
sum1 += matrices1[w + 258] * partials1[v + 18];
sum2 += matrices2[w + 258] * partials2[v + 18];
sum1 += matrices1[w + 259] * partials1[v + 19];
sum2 += matrices2[w + 259] * partials2[v + 19];
partials3[u] = sum1 * sum2; u++;

sum1 = matrices1[w + 260] * partials1[v];
sum2 = matrices2[w + 260] * partials2[v];
sum1 += matrices1[w + 261] * partials1[v + 1];
sum2 += matrices2[w + 261] * partials2[v + 1];
sum1 += matrices1[w + 262] * partials1[v + 2];
sum2 += matrices2[w + 262] * partials2[v + 2];
sum1 += matrices1[w + 263] * partials1[v + 3];
sum2 += matrices2[w + 263] * partials2[v + 3];
sum1 += matrices1[w + 264] * partials1[v + 4];
sum2 += matrices2[w + 264] * partials2[v + 4];
sum1 += matrices1[w + 265] * partials1[v + 5];
sum2 += matrices2[w + 265] * partials2[v + 5];
sum1 += matrices1[w + 266] * partials1[v + 6];
sum2 += matrices2[w + 266] * partials2[v + 6];
sum1 += matrices1[w + 267] * partials1[v + 7];
sum2 += matrices2[w + 267] * partials2[v + 7];
sum1 += matrices1[w + 268] * partials1[v + 8];
sum2 += matrices2[w + 268] * partials2[v + 8];
sum1 += matrices1[w + 269] * partials1[v + 9];
sum2 += matrices2[w + 269] * partials2[v + 9];
sum1 += matrices1[w + 270] * partials1[v + 10];
sum2 += matrices2[w + 270] * partials2[v + 10];
sum1 += matrices1[w + 271] * partials1[v + 11];
sum2 += matrices2[w + 271] * partials2[v + 11];
sum1 += matrices1[w + 272] * partials1[v + 12];
sum2 += matrices2[w + 272] * partials2[v + 12];
sum1 += matrices1[w + 273] * partials1[v + 13];
sum2 += matrices2[w + 273] * partials2[v + 13];
sum1 += matrices1[w + 274] * partials1[v + 14];
sum2 += matrices2[w + 274] * partials2[v + 14];
sum1 += matrices1[w + 275] * partials1[v + 15];
sum2 += matrices2[w + 275] * partials2[v + 15];
sum1 += matrices1[w + 276] * partials1[v + 16];
sum2 += matrices2[w + 276] * partials2[v + 16];
sum1 += matrices1[w + 277] * partials1[v + 17];
sum2 += matrices2[w + 277] * partials2[v + 17];
sum1 += matrices1[w + 278] * partials1[v + 18];
sum2 += matrices2[w + 278] * partials2[v + 18];
sum1 += matrices1[w + 279] * partials1[v + 19];
sum2 += matrices2[w + 279] * partials2[v + 19];
partials3[u] = sum1 * sum2; u++;

sum1 = matrices1[w + 280] * partials1[v];
sum2 = matrices2[w + 280] * partials2[v];
sum1 += matrices1[w + 281] * partials1[v + 1];
sum2 += matrices2[w + 281] * partials2[v + 1];
sum1 += matrices1[w + 282] * partials1[v + 2];
sum2 += matrices2[w + 282] * partials2[v + 2];
sum1 += matrices1[w + 283] * partials1[v + 3];
sum2 += matrices2[w + 283] * partials2[v + 3];
sum1 += matrices1[w + 284] * partials1[v + 4];
sum2 += matrices2[w + 284] * partials2[v + 4];
sum1 += matrices1[w + 285] * partials1[v + 5];
sum2 += matrices2[w + 285] * partials2[v + 5];
sum1 += matrices1[w + 286] * partials1[v + 6];
sum2 += matrices2[w + 286] * partials2[v + 6];
sum1 += matrices1[w + 287] * partials1[v + 7];
sum2 += matrices2[w + 287] * partials2[v + 7];
sum1 += matrices1[w + 288] * partials1[v + 8];
sum2 += matrices2[w + 288] * partials2[v + 8];
sum1 += matrices1[w + 289] * partials1[v + 9];
sum2 += matrices2[w + 289] * partials2[v + 9];
sum1 += matrices1[w + 290] * partials1[v + 10];
sum2 += matrices2[w + 290] * partials2[v + 10];
sum1 += matrices1[w + 291] * partials1[v + 11];
sum2 += matrices2[w + 291] * partials2[v + 11];
sum1 += matrices1[w + 292] * partials1[v + 12];
sum2 += matrices2[w + 292] * partials2[v + 12];
sum1 += matrices1[w + 293] * partials1[v + 13];
sum2 += matrices2[w + 293] * partials2[v + 13];
sum1 += matrices1[w + 294] * partials1[v + 14];
sum2 += matrices2[w + 294] * partials2[v + 14];
sum1 += matrices1[w + 295] * partials1[v + 15];
sum2 += matrices2[w + 295] * partials2[v + 15];
sum1 += matrices1[w + 296] * partials1[v + 16];
sum2 += matrices2[w + 296] * partials2[v + 16];
sum1 += matrices1[w + 297] * partials1[v + 17];
sum2 += matrices2[w + 297] * partials2[v + 17];
sum1 += matrices1[w + 298] * partials1[v + 18];
sum2 += matrices2[w + 298] * partials2[v + 18];
sum1 += matrices1[w + 299] * partials1[v + 19];
sum2 += matrices2[w + 299] * partials2[v + 19];
partials3[u] = sum1 * sum2; u++;

sum1 = matrices1[w + 300] * partials1[v];
sum2 = matrices2[w + 300] * partials2[v];
sum1 += matrices1[w + 301] * partials1[v + 1];
sum2 += matrices2[w + 301] * partials2[v + 1];
sum1 += matrices1[w + 302] * partials1[v + 2];
sum2 += matrices2[w + 302] * partials2[v + 2];
sum1 += matrices1[w + 303] * partials1[v + 3];
sum2 += matrices2[w + 303] * partials2[v + 3];
sum1 += matrices1[w + 304] * partials1[v + 4];
sum2 += matrices2[w + 304] * partials2[v + 4];
sum1 += matrices1[w + 305] * partials1[v + 5];
sum2 += matrices2[w + 305] * partials2[v + 5];
sum1 += matrices1[w + 306] * partials1[v + 6];
sum2 += matrices2[w + 306] * partials2[v + 6];
sum1 += matrices1[w + 307] * partials1[v + 7];
sum2 += matrices2[w + 307] * partials2[v + 7];
sum1 += matrices1[w + 308] * partials1[v + 8];
sum2 += matrices2[w + 308] * partials2[v + 8];
sum1 += matrices1[w + 309] * partials1[v + 9];
sum2 += matrices2[w + 309] * partials2[v + 9];
sum1 += matrices1[w + 310] * partials1[v + 10];
sum2 += matrices2[w + 310] * partials2[v + 10];
sum1 += matrices1[w + 311] * partials1[v + 11];
sum2 += matrices2[w + 311] * partials2[v + 11];
sum1 += matrices1[w + 312] * partials1[v + 12];
sum2 += matrices2[w + 312] * partials2[v + 12];
sum1 += matrices1[w + 313] * partials1[v + 13];
sum2 += matrices2[w + 313] * partials2[v + 13];
sum1 += matrices1[w + 314] * partials1[v + 14];
sum2 += matrices2[w + 314] * partials2[v + 14];
sum1 += matrices1[w + 315] * partials1[v + 15];
sum2 += matrices2[w + 315] * partials2[v + 15];
sum1 += matrices1[w + 316] * partials1[v + 16];
sum2 += matrices2[w + 316] * partials2[v + 16];
sum1 += matrices1[w + 317] * partials1[v + 17];
sum2 += matrices2[w + 317] * partials2[v + 17];
sum1 += matrices1[w + 318] * partials1[v + 18];
sum2 += matrices2[w + 318] * partials2[v + 18];
sum1 += matrices1[w + 319] * partials1[v + 19];
sum2 += matrices2[w + 319] * partials2[v + 19];
partials3[u] = sum1 * sum2; u++;

sum1 = matrices1[w + 320] * partials1[v];
sum2 = matrices2[w + 320] * partials2[v];
sum1 += matrices1[w + 321] * partials1[v + 1];
sum2 += matrices2[w + 321] * partials2[v + 1];
sum1 += matrices1[w + 322] * partials1[v + 2];
sum2 += matrices2[w + 322] * partials2[v + 2];
sum1 += matrices1[w + 323] * partials1[v + 3];
sum2 += matrices2[w + 323] * partials2[v + 3];
sum1 += matrices1[w + 324] * partials1[v + 4];
sum2 += matrices2[w + 324] * partials2[v + 4];
sum1 += matrices1[w + 325] * partials1[v + 5];
sum2 += matrices2[w + 325] * partials2[v + 5];
sum1 += matrices1[w + 326] * partials1[v + 6];
sum2 += matrices2[w + 326] * partials2[v + 6];
sum1 += matrices1[w + 327] * partials1[v + 7];
sum2 += matrices2[w + 327] * partials2[v + 7];
sum1 += matrices1[w + 328] * partials1[v + 8];
sum2 += matrices2[w + 328] * partials2[v + 8];
sum1 += matrices1[w + 329] * partials1[v + 9];
sum2 += matrices2[w + 329] * partials2[v + 9];
sum1 += matrices1[w + 330] * partials1[v + 10];
sum2 += matrices2[w + 330] * partials2[v + 10];
sum1 += matrices1[w + 331] * partials1[v + 11];
sum2 += matrices2[w + 331] * partials2[v + 11];
sum1 += matrices1[w + 332] * partials1[v + 12];
sum2 += matrices2[w + 332] * partials2[v + 12];
sum1 += matrices1[w + 333] * partials1[v + 13];
sum2 += matrices2[w + 333] * partials2[v + 13];
sum1 += matrices1[w + 334] * partials1[v + 14];
sum2 += matrices2[w + 334] * partials2[v + 14];
sum1 += matrices1[w + 335] * partials1[v + 15];
sum2 += matrices2[w + 335] * partials2[v + 15];
sum1 += matrices1[w + 336] * partials1[v + 16];
sum2 += matrices2[w + 336] * partials2[v + 16];
sum1 += matrices1[w + 337] * partials1[v + 17];
sum2 += matrices2[w + 337] * partials2[v + 17];
sum1 += matrices1[w + 338] * partials1[v + 18];
sum2 += matrices2[w + 338] * partials2[v + 18];
sum1 += matrices1[w + 339] * partials1[v + 19];
sum2 += matrices2[w + 339] * partials2[v + 19];
partials3[u] = sum1 * sum2; u++;

sum1 = matrices1[w + 340] * partials1[v];
sum2 = matrices2[w + 340] * partials2[v];
sum1 += matrices1[w + 341] * partials1[v + 1];
sum2 += matrices2[w + 341] * partials2[v + 1];
sum1 += matrices1[w + 342] * partials1[v + 2];
sum2 += matrices2[w + 342] * partials2[v + 2];
sum1 += matrices1[w + 343] * partials1[v + 3];
sum2 += matrices2[w + 343] * partials2[v + 3];
sum1 += matrices1[w + 344] * partials1[v + 4];
sum2 += matrices2[w + 344] * partials2[v + 4];
sum1 += matrices1[w + 345] * partials1[v + 5];
sum2 += matrices2[w + 345] * partials2[v + 5];
sum1 += matrices1[w + 346] * partials1[v + 6];
sum2 += matrices2[w + 346] * partials2[v + 6];
sum1 += matrices1[w + 347] * partials1[v + 7];
sum2 += matrices2[w + 347] * partials2[v + 7];
sum1 += matrices1[w + 348] * partials1[v + 8];
sum2 += matrices2[w + 348] * partials2[v + 8];
sum1 += matrices1[w + 349] * partials1[v + 9];
sum2 += matrices2[w + 349] * partials2[v + 9];
sum1 += matrices1[w + 350] * partials1[v + 10];
sum2 += matrices2[w + 350] * partials2[v + 10];
sum1 += matrices1[w + 351] * partials1[v + 11];
sum2 += matrices2[w + 351] * partials2[v + 11];
sum1 += matrices1[w + 352] * partials1[v + 12];
sum2 += matrices2[w + 352] * partials2[v + 12];
sum1 += matrices1[w + 353] * partials1[v + 13];
sum2 += matrices2[w + 353] * partials2[v + 13];
sum1 += matrices1[w + 354] * partials1[v + 14];
sum2 += matrices2[w + 354] * partials2[v + 14];
sum1 += matrices1[w + 355] * partials1[v + 15];
sum2 += matrices2[w + 355] * partials2[v + 15];
sum1 += matrices1[w + 356] * partials1[v + 16];
sum2 += matrices2[w + 356] * partials2[v + 16];
sum1 += matrices1[w + 357] * partials1[v + 17];
sum2 += matrices2[w + 357] * partials2[v + 17];
sum1 += matrices1[w + 358] * partials1[v + 18];
sum2 += matrices2[w + 358] * partials2[v + 18];
sum1 += matrices1[w + 359] * partials1[v + 19];
sum2 += matrices2[w + 359] * partials2[v + 19];
partials3[u] = sum1 * sum2; u++;

sum1 = matrices1[w + 360] * partials1[v];
sum2 = matrices2[w + 360] * partials2[v];
sum1 += matrices1[w + 361] * partials1[v + 1];
sum2 += matrices2[w + 361] * partials2[v + 1];
sum1 += matrices1[w + 362] * partials1[v + 2];
sum2 += matrices2[w + 362] * partials2[v + 2];
sum1 += matrices1[w + 363] * partials1[v + 3];
sum2 += matrices2[w + 363] * partials2[v + 3];
sum1 += matrices1[w + 364] * partials1[v + 4];
sum2 += matrices2[w + 364] * partials2[v + 4];
sum1 += matrices1[w + 365] * partials1[v + 5];
sum2 += matrices2[w + 365] * partials2[v + 5];
sum1 += matrices1[w + 366] * partials1[v + 6];
sum2 += matrices2[w + 366] * partials2[v + 6];
sum1 += matrices1[w + 367] * partials1[v + 7];
sum2 += matrices2[w + 367] * partials2[v + 7];
sum1 += matrices1[w + 368] * partials1[v + 8];
sum2 += matrices2[w + 368] * partials2[v + 8];
sum1 += matrices1[w + 369] * partials1[v + 9];
sum2 += matrices2[w + 369] * partials2[v + 9];
sum1 += matrices1[w + 370] * partials1[v + 10];
sum2 += matrices2[w + 370] * partials2[v + 10];
sum1 += matrices1[w + 371] * partials1[v + 11];
sum2 += matrices2[w + 371] * partials2[v + 11];
sum1 += matrices1[w + 372] * partials1[v + 12];
sum2 += matrices2[w + 372] * partials2[v + 12];
sum1 += matrices1[w + 373] * partials1[v + 13];
sum2 += matrices2[w + 373] * partials2[v + 13];
sum1 += matrices1[w + 374] * partials1[v + 14];
sum2 += matrices2[w + 374] * partials2[v + 14];
sum1 += matrices1[w + 375] * partials1[v + 15];
sum2 += matrices2[w + 375] * partials2[v + 15];
sum1 += matrices1[w + 376] * partials1[v + 16];
sum2 += matrices2[w + 376] * partials2[v + 16];
sum1 += matrices1[w + 377] * partials1[v + 17];
sum2 += matrices2[w + 377] * partials2[v + 17];
sum1 += matrices1[w + 378] * partials1[v + 18];
sum2 += matrices2[w + 378] * partials2[v + 18];
sum1 += matrices1[w + 379] * partials1[v + 19];
sum2 += matrices2[w + 379] * partials2[v + 19];
partials3[u] = sum1 * sum2; u++;

sum1 = matrices1[w + 380] * partials1[v];
sum2 = matrices2[w + 380] * partials2[v];
sum1 += matrices1[w + 381] * partials1[v + 1];
sum2 += matrices2[w + 381] * partials2[v + 1];
sum1 += matrices1[w + 382] * partials1[v + 2];
sum2 += matrices2[w + 382] * partials2[v + 2];
sum1 += matrices1[w + 383] * partials1[v + 3];
sum2 += matrices2[w + 383] * partials2[v + 3];
sum1 += matrices1[w + 384] * partials1[v + 4];
sum2 += matrices2[w + 384] * partials2[v + 4];
sum1 += matrices1[w + 385] * partials1[v + 5];
sum2 += matrices2[w + 385] * partials2[v + 5];
sum1 += matrices1[w + 386] * partials1[v + 6];
sum2 += matrices2[w + 386] * partials2[v + 6];
sum1 += matrices1[w + 387] * partials1[v + 7];
sum2 += matrices2[w + 387] * partials2[v + 7];
sum1 += matrices1[w + 388] * partials1[v + 8];
sum2 += matrices2[w + 388] * partials2[v + 8];
sum1 += matrices1[w + 389] * partials1[v + 9];
sum2 += matrices2[w + 389] * partials2[v + 9];
sum1 += matrices1[w + 390] * partials1[v + 10];
sum2 += matrices2[w + 390] * partials2[v + 10];
sum1 += matrices1[w + 391] * partials1[v + 11];
sum2 += matrices2[w + 391] * partials2[v + 11];
sum1 += matrices1[w + 392] * partials1[v + 12];
sum2 += matrices2[w + 392] * partials2[v + 12];
sum1 += matrices1[w + 393] * partials1[v + 13];
sum2 += matrices2[w + 393] * partials2[v + 13];
sum1 += matrices1[w + 394] * partials1[v + 14];
sum2 += matrices2[w + 394] * partials2[v + 14];
sum1 += matrices1[w + 395] * partials1[v + 15];
sum2 += matrices2[w + 395] * partials2[v + 15];
sum1 += matrices1[w + 396] * partials1[v + 16];
sum2 += matrices2[w + 396] * partials2[v + 16];
sum1 += matrices1[w + 397] * partials1[v + 17];
sum2 += matrices2[w + 397] * partials2[v + 17];
sum1 += matrices1[w + 398] * partials1[v + 18];
sum2 += matrices2[w + 398] * partials2[v + 18];
sum1 += matrices1[w + 399] * partials1[v + 19];
sum2 += matrices2[w + 399] * partials2[v + 19];
partials3[u] = sum1 * sum2; u++;

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

/*
 * Class:     dr_evomodel_treelikelihood_NativeAminoAcidLikelihoodCore
 * Method:    nativeStatesPartialsPruning
 * Signature: ([I[D[D[DII[D)V
 */
/*
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
}*/
							

/*
 * Class:     dr_evomodel_treelikelihood_NativeAminoAcidLikelihoodCore
 * Method:    nativePartialsPartialsPruning
 * Signature: ([D[D[D[DII[D)V
 */
/*
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
}*/
