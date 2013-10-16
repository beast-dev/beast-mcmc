#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "NativeSubstitutionModel.h"

//#define STATECOUNT	4

int stateCount = 0;
int stateCountSquared = 0;

#define REAL		double
#define SIZE_REAL	sizeof(REAL)

#define MINPROB		1E-10f;

#define MEMCPY(to,from,length,toType)	   { int m; \
										for(m=0; m<length; m++) { \
											to[m] = (toType) from[m]; \
										}}


JNIEXPORT jlong JNICALL Java_dr_evomodel_substmodel_NativeSubstitutionModel_allocateNativeMemoryArray
  (JNIEnv *env, jobject obj, jint length) {
	REAL *ptr = (REAL *)malloc(length * SIZE_REAL);
	return (jlong) ptr;
  }
  
  JNIEXPORT jint JNICALL Java_dr_evomodel_substmodel_NativeSubstitutionModel_getNativeRealSize
  (JNIEnv *env, jobject obj) {
	return (jint) SIZE_REAL;
}
  
JNIEXPORT void JNICALL Java_dr_evomodel_substmodel_NativeSubstitutionModel_copyNativeMemoryArray
  (JNIEnv *env, jobject obj, jlong inFrom, jlong inTo, jint length) {
  
  REAL *from = (REAL *) inFrom;
  REAL *to = (REAL *) inTo;
  memcpy(to,from,length*SIZE_REAL);
} 
 

JNIEXPORT void JNICALL Java_dr_evomodel_substmodel_NativeSubstitutionModel_nativeStoreState
  (JNIEnv *env, jobject obj) {}

JNIEXPORT void JNICALL Java_dr_evomodel_substmodel_NativeSubstitutionModel_nativeRestoreState
  (JNIEnv *env, jobject obj) {}




REAL *deviceCache(int stateCount) {
	return (REAL *)malloc(stateCount*stateCount*SIZE_REAL);
}

  
JNIEXPORT void JNICALL Java_dr_evomodel_substmodel_NativeSubstitutionModel_nativeGetTransitionProbabilities
   (JNIEnv *env, jobject obj, jlong inCache, jdouble distance, jlong inPtr, jint states) {
	
	REAL *matrix, *iexp, *Ievc, *Evec, *Eval, *EvalImag;
	int stateCount = states;
	
	Ievc  = (REAL *) inCache;	
	Evec = Ievc + states*stateCount;
	Eval = Evec + stateCount*stateCount;
	EvalImag = Eval + stateCount;
	
	matrix = (REAL *) inPtr;
	
	iexp = deviceCache(stateCount);
		
	int i,j,k;
	
	
		for (i = 0; i < stateCount; i++) {

			if (EvalImag[i] == 0) {
				// 1x1 block
				REAL temp = exp(distance * Eval[i]);
				for (j = 0; j < stateCount; j++) {
					iexp[i*stateCount + j] = Ievc[i*stateCount + j] * temp;
				}
			} else {
				// 2x2 conjugate block
				// If A is 2x2 with complex conjugate pair eigenvalues a +/- bi, then
				// exp(At) = exp(at)*( cos(bt)I + \frac{sin(bt)}{b}(A - aI)).
				int i2 = i + 1;
				REAL b = EvalImag[i];
				REAL expat = exp(distance * Eval[i]);
				REAL expatcosbt = expat * cos(distance * b);
				REAL expatsinbt = expat * sin(distance * b);

				for (j = 0; j < stateCount; j++) {
					iexp[i*stateCount + j] = expatcosbt * Ievc[i*stateCount + j] + expatsinbt * Ievc[i2*stateCount + j];
					iexp[i2*stateCount + j] = expatcosbt * Ievc[i2*stateCount + j] - expatsinbt * Ievc[i*stateCount + j];
				}
				i++; // processed two conjugate rows
			}
		}

		int u = 0;
		for (i = 0; i < stateCount; i++) {
			for (j = 0; j < stateCount; j++) {
				REAL temp = 0.0;
				for (k = 0; k < stateCount; k++) {
					temp += Evec[i*stateCount + k] * iexp[k*stateCount + j];
				}
				if (temp < 0.0)
					matrix[u] = 1E-10;
				else
					matrix[u] = temp;
				u++;
			}
		}
	
	free(iexp);
	
   }
  
  JNIEXPORT void JNICALL Java_dr_evomodel_substmodel_NativeSubstitutionModel_getNativeMemoryArray
  (JNIEnv *env, jobject obj, jlong inNativePtr, jint fromOffset, jdoubleArray inToJavaArray, jint toOffset, jint length) 
  {

	jdouble *toJavaArray = (jdouble*)(*env)->GetPrimitiveArrayCritical(env, inToJavaArray, 0);
	jdouble *toOffsetJavaArray = toJavaArray + toOffset;
  
	REAL *fromPtr = (REAL *)inNativePtr + fromOffset;
	MEMCPY(toOffsetJavaArray,fromPtr,length,jdouble);
  
	(*env)->ReleasePrimitiveArrayCritical(env, inToJavaArray, toJavaArray, 0); // need to save, not abort
	
  }


JNIEXPORT void JNICALL Java_dr_evomodel_substmodel_NativeSubstitutionModel_nativeSetup
  (JNIEnv *env, jobject obj, jlong inPtr, jobjectArray inIevc, jobjectArray inEvec, jdoubleArray inEval, 
  jdoubleArray inEvalImag, jint stateCount) {
  
	if( inPtr == 0 ) {
		printf("Null pointer sent.\n");
		exit(-1);
	}
	
	REAL *Ievc, *tmpIevc, *Evec, *tmpEvec, *Eval, *EvalImag;
	
	tmpIevc = Ievc = (REAL *) inPtr;
	tmpEvec = Evec = Ievc + stateCount*stateCount;
	Eval = Evec + stateCount*stateCount;
	EvalImag = Eval + stateCount;
	
	jdouble *javaEval = (jdouble*)(*env)->GetPrimitiveArrayCritical(env,inEval,0);
	jdouble *javaEvalImag = (jdouble*)(*env)->GetPrimitiveArrayCritical(env,inEvalImag,0);
	
	int i,j;
	for(i=0; i<stateCount; i++) {
		jdoubleArray oneDim = (*env)->GetObjectArrayElement(env,inIevc,i);
		jdouble *element = (jdouble*)(*env)->GetPrimitiveArrayCritical(env,oneDim,0);
		MEMCPY(tmpIevc,element,stateCount,REAL);
		tmpIevc += stateCount;			
		(*env)->ReleasePrimitiveArrayCritical(env,oneDim,element,JNI_ABORT);
		
		oneDim = (*env)->GetObjectArrayElement(env,inEvec,i);
		element = (jdouble*)(*env)->GetPrimitiveArrayCritical(env,oneDim,0);
		MEMCPY(tmpEvec,element,stateCount,REAL);
		tmpEvec += stateCount;
		(*env)->ReleasePrimitiveArrayCritical(env,oneDim,element,JNI_ABORT);
		
	}
	
	MEMCPY(Eval,javaEval,stateCount,REAL);
	MEMCPY(EvalImag, javaEvalImag,stateCount,REAL);
	
	(*env)->ReleasePrimitiveArrayCritical(env,inEval, javaEval, JNI_ABORT);
	(*env)->ReleasePrimitiveArrayCritical(env,inEvalImag, javaEvalImag, JNI_ABORT);
	
}

	 
