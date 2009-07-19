#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "NativeMemoryLikelihoodCore.h"

/**** COMMENT OUT TO TURN OFF DEBUGGING ****/
//#define DEBUG

#define DLIM	100

/* Definition of REAL can be switched between 'double' and 'float' */

#define REAL		double
#define SIZE_REAL	sizeof(REAL)
#define INT		int
#define SIZE_INT	sizeof(int)


JNIEXPORT jint JNICALL Java_dr_evomodel_treelikelihood_NativeMemoryLikelihoodCore_getNativeRealSize
  (JNIEnv *env, jobject obj) {
	return (jint) SIZE_REAL;
}

JNIEXPORT jint JNICALL Java_dr_evomodel_treelikelihood_NativeMemoryLikelihoodCore_getNativeIntSize
  (JNIEnv *env, jobject obj) {
	return (jint) SIZE_INT;
}

JNIEXPORT jlong JNICALL Java_dr_evomodel_treelikelihood_NativeMemoryLikelihoodCore_allocateNativeMemoryArray
  (JNIEnv *env, jobject obj, jint length) 
{  



#ifdef DEBUG
	fprintf(stderr,"Entering ANMA\n");
#endif

  REAL *data = (REAL *)malloc(SIZE_REAL * length);
  if (data == NULL) {
	fprintf(stderr,"Failed to allocate memory!");
	exit(-1);
  }
  
#ifdef DEBUG
	fprintf(stderr,"Leaving ANMA\n");
#endif  
  
  return (jlong) data;
} 

JNIEXPORT jlong JNICALL Java_dr_evomodel_treelikelihood_NativeMemoryLikelihoodCore_allocateNativeIntMemoryArray
  (JNIEnv *env, jobject obj, jint length) 
{  

#ifdef DEBUG
	fprintf(stderr,"Entering ANMA\n");
#endif

  INT *data = (INT *)malloc(SIZE_INT * length);
  if (data == NULL) {
	fprintf(stderr,"Failed to allocate memory!");
	exit(-1);
  }
  
#ifdef DEBUG
	fprintf(stderr,"Leaving ANMA\n");
#endif  
  
  return (jlong) data;
} 


JNIEXPORT void JNICALL Java_dr_evomodel_treelikelihood_NativeMemoryLikelihoodCore_freeNativeMemoryArray
  (JNIEnv *env, jobject obj, jlong inNativePtr)
  {

#ifdef DEBUG
	fprintf(stderr,"Entering FNMA\n");
#endif
  
	if (inNativePtr != NULL) {	
		free((void*)inNativePtr);
	}  

#ifdef DEBUG
	fprintf(stderr,"Leaving ANMA\n");
#endif	
	
  }

	
JNIEXPORT void JNICALL Java_dr_evomodel_treelikelihood_NativeMemoryLikelihoodCore_setNativeMemoryArray___3DIJII
  (JNIEnv *env, jobject obj, jdoubleArray inFromJavaArray, jint fromOffset, jlong inNativePtr, jint toOffset, jint length) 
  {
  
#ifdef DEBUG
	fprintf(stderr,"Entering SNMA\n");
#endif  

	jdouble *fromJavaArray = (jdouble*)(*env)->GetPrimitiveArrayCritical(env, inFromJavaArray, 0);
	
	//	if (inNativePtr == 0) {
	//		fprintf(stderr,"Throwing null pointer native exception\n");
	//		jclass nullPtr = (*env)->FindClass(env,"java/lang/NullPointerException");
	//		(*env)->ThrowNew(env,nullPtr,"Null native pointer in setNativeMemoryArray(double[])");
	//	}
	
	// copy entry-by-entry to ensure conversion to REAL
	REAL *toPtr = (REAL *)(inNativePtr + SIZE_REAL*toOffset);
	int i;
	for(i=0; i<length; i++) {
		toPtr[i] = (REAL) fromJavaArray[i+fromOffset]; // Can speed-up by moving starting pointers by offset, TODO
	}
	
  	(*env)->ReleasePrimitiveArrayCritical(env, inFromJavaArray, fromJavaArray, JNI_ABORT);  

#ifdef DEBUG
	fprintf(stderr,"Leaving SNMA\n");
#endif	
	
  }
  
JNIEXPORT void JNICALL Java_dr_evomodel_treelikelihood_NativeMemoryLikelihoodCore_setNativeMemoryArray___3IIJII
  (JNIEnv *env, jobject obj, jintArray inFromJavaArray, jint fromOffset, jlong inNativePtr, jint toOffset, jint length) 
  {  
  
#ifdef DEBUG
	fprintf(stderr,"Entering SNMA\n");
#endif  
  
	jint *fromJavaArray = (jint*)(*env)->GetPrimitiveArrayCritical(env, inFromJavaArray, 0);
	
	if (inNativePtr == 0) {
		fprintf(stderr,"Throwing null pointer native exception\n");
		jclass nullPtr = (*env)->FindClass(env,"java/lang/NullPointerException");
		(*env)->ThrowNew(env,nullPtr,"Null native pointer in setNativeMemoryArray(double[])");
	}
	// copy entry-by-entry to ensure conversion to REAL
	INT *toPtr = (INT *)inNativePtr;
	int i;
	for(i=0; i<length; i++) {
		toPtr[i+toOffset] = (INT) fromJavaArray[i+fromOffset];  // Can speed-up by moving starting pointers by offset, TODO
	}
	
  	(*env)->ReleasePrimitiveArrayCritical(env, inFromJavaArray, fromJavaArray, JNI_ABORT);  

#ifdef DEBUG
	fprintf(stderr,"Leaving SNMA\n");
#endif	
	
  }

JNIEXPORT void JNICALL Java_dr_evomodel_treelikelihood_NativeMemoryLikelihoodCore_getNativeMemoryArray__JI_3DII
  (JNIEnv *env, jobject obj, jlong inNativePtr, jint fromOffset, jdoubleArray inToJavaArray, jint toOffset, jint length) 
  {

#ifdef DEBUG  
	fprintf(stderr,"Entering GA\n");
#endif 

	jdouble *toJavaArray = (jdouble*)(*env)->GetPrimitiveArrayCritical(env, inToJavaArray, 0);
  
	// copy entry-by-entry to ensure conversion to Java type
	REAL *fromPtr = (REAL *)inNativePtr;
	int i;
	for(i=0; i<length; i++) {
		toJavaArray[i+toOffset] = (jdouble) fromPtr[i+fromOffset];
	}
  
	(*env)->ReleasePrimitiveArrayCritical(env, inToJavaArray, toJavaArray, 0); // need to save, not abort
	
#ifdef DEBUG  
	fprintf(stderr,"Leaving GA\n");
#endif	

  }

JNIEXPORT void JNICALL Java_dr_evomodel_treelikelihood_NativeMemoryLikelihoodCore_getNativeMemoryArray__JI_3III
  (JNIEnv *env, jobject obj, jlong inNativePtr, jint fromOffset, jintArray inToJavaArray, jint toOffset, jint length)
  {
  	jint *toJavaArray = (jint*)(*env)->GetPrimitiveArrayCritical(env, inToJavaArray, 0);
	
#ifdef DEBUG  
	fprintf(stderr,"Entering GA\n");
//	if (inNativePtr < 100) fprintf(stderr, "inNativePtr is null in getNativeMemoryArray\n");
//	if (inNativePtr == 0) {
//		fprintf(stderr,"Throwing null pointer native exception\n");
//		jclass nullPtr = (*env)->FindClass(env,"java/lang/NullPointerException");
//		(*env)->ThrowNew(env,nullPtr,"Null native pointer in setNativeMemoryArray(double[])");
//	}
#endif	
  	  
	// copy entry-by-entry to ensure conversion to Java type
	INT *fromPtr = (INT *)inNativePtr;
	int i;
	for(i=0; i<length; i++) {
		toJavaArray[i+toOffset] = (jint) fromPtr[i+fromOffset]; // TODO precompute offsets
	}
  
	(*env)->ReleasePrimitiveArrayCritical(env, inToJavaArray, toJavaArray, 0); // need to save, not abort
	
#ifdef DEBUG  
	fprintf(stderr,"Leaving GA\n");
#endif
  
  }



JNIEXPORT void JNICALL Java_dr_evomodel_treelikelihood_NativeMemoryLikelihoodCore_nativePartialsPartialsPruning
  (JNIEnv *env, jobject obj, jlong inPartials1, jlong inMatrices1, jlong inPartials2, jlong inMatrices2, 
  jint patternCount, jint matrixCount, jlong outPartials, jint stateCount)
{
	REAL *partials1 = (REAL *)inPartials1;
	REAL *matrices1 = (REAL *)inMatrices1;
	REAL *partials2 = (REAL *)inPartials2;
	REAL *matrices2 = (REAL *)inMatrices2;
	REAL *partials3 = (REAL *)outPartials;
	
#ifdef DEBUG
	fprintf(stderr,"Entering PP\n");
//		fprintf(stderr,"Checking PartialsPartials\n");
//		if (partials1 < 100)
//			fprintf(stderr,"partials1 is null in PartialsPartials\n");
//		if (partials2 < 100)
//			fprintf(stderr,"partials2 is null in PartialsPartials\n");
//		if (partials3 < 100)
//			fprintf(stderr,"partials3 is null in PartialsPartials\n");
//		if (matrices1 < 100)
//			fprintf(stderr,"matrices1 is null in PartialsPartials\n");
//		if (matrices2 < 100)
//			fprintf(stderr,"matrices2 is null in PartialsPartials\n");
#endif
						
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
	
#ifdef DEBUG
	fprintf(stderr,"Leaving SP\n");
#endif	
}


JNIEXPORT void JNICALL Java_dr_evomodel_treelikelihood_NativeMemoryLikelihoodCore_nativeStatesPartialsPruning
  (JNIEnv *env, jobject obj, jlong inStates1, jlong inMatrices1, jlong inPartials2, jlong inMatrices2, 
  jint patternCount, jint matrixCount, jlong outPartials, jint stateCount) 
{
	INT  *states1   = (INT *) inStates1;
	REAL *matrices1 = (REAL *)inMatrices1;
	REAL *partials2 = (REAL *)inPartials2;
	REAL *matrices2 = (REAL *)inMatrices2;
	REAL *partials3 = (REAL *)outPartials;
	
#ifdef DEBUG
	fprintf(stderr,"Entering SP\n");
	if (states1 < DLIM)
		fprintf(stderr,"states1\n");
	if (matrices1 < DLIM)
		fprintf(stderr,"matrices1\n");	
	if (partials2 < DLIM)
		fprintf(stderr,"partials2\n");
	if (matrices2 < DLIM)
		fprintf(stderr,"matrices2\n");
	if (partials3 < DLIM)
		fprintf(stderr,"partials3\n");
		
	fprintf(stderr,"states1 = %d\nmatrices1 = %d\npartials2 = %d\nmatrices2 = %d\npartials3 = %d\n",
		states1,matrices1,partials2,matrices2,partials3);
	fprintf(stderr,"<done>\n");
#endif	
	
	int i, l, k, state1, m;
	REAL sum;
	int u = 0;
	int v = 0;
	int w = 0;
	int x, y;

	for (l = 0; l < matrixCount; l++) {
		for (k = 0; k < patternCount; k++) {
		
			int state1 = states1[k];

#ifdef DEBUG			
			fprintf(stderr,"Read datum %d for m = %d, p = %d\n",state1,l,k);
#endif			
				if (state1 < stateCount) {
			
					x = w;
					for (i = 0; i < stateCount; i++) {
					
						y = v;
						REAL value = matrices1[x + state1];	
#ifdef DEBUG
						fprintf(stderr,"Read m1 %1.2e at index %d\n",value,(x+state1));
#endif						
						
						sum = 0;
						for (m = 0; m < stateCount; m++) {
							sum += matrices2[x] * partials2[y]; 
#ifdef DEBUG
							fprintf(stderr,"Read m2 %1.2e at index %d and p2 %1.2e at index %d\n",matrices2[x],x,partials2[y],y);
#endif							
							x++; y++;
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
		
#ifdef DEBUG
	fprintf(stderr,"Leaving SP\n");
#endif			
		
}



JNIEXPORT void JNICALL Java_dr_evomodel_treelikelihood_NativeMemoryLikelihoodCore_nativeIntegratePartials
  (JNIEnv *env, jobject obj, jlong inPartials, jdoubleArray inProportions, jint patternCount, jint matrixCount,
   jdoubleArray outPartials, jint stateCount)
{
	REAL *partials1 = (REAL *) inPartials;
	
#ifdef DEBUG
	fprintf(stderr,"Entering IP\n");
//		fprintf(stderr,"Checking IntegratePartials\n");
//		if (partials1 < 100) fprintf(stderr,"partials1 in IntegratePartials is null\n");
//		fprintf(stderr,"Done IntegratePartials\n");
#endif
	
//	REAL *tmp = reinterpret_cast<REAL *>(inPartials);  // This statement should work without warning in C++
	
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
	
#ifdef DEBUG
	fprintf(stderr,"Leaving IP\n");
#endif
}


JNIEXPORT void JNICALL Java_dr_evomodel_treelikelihood_NativeMemoryLikelihoodCore_nativeStatesStatesPruning
  (JNIEnv *env, jobject obj, jlong inStates1, jlong inMatrices1, jlong inStates2, jlong inMatrices2, jint patternCount,
   jint matrixCount, jlong inPartials3, jint stateCount) 
{
	int j, k, u, v, w, state1, state2, m;
	
	INT  *states1   = (INT *)  inStates1;
	INT  *states2   = (INT *)  inStates2;
	
	REAL *matrices1 = (REAL *) inMatrices1;
	REAL *matrices2 = (REAL *) inMatrices2;
	
	REAL *partials3 = (REAL *) inPartials3;
	
#ifdef DEBUG
	fprintf(stderr,"Entering SS\n");
#endif	
	
//#ifdef DEBUG
//		fprintf(stderr,"Checking StatesStates\n");
//		if (states1 < 100)
//			fprintf(stderr,"states1 is null in StatesStates\n");
//		if (states2 < 100)
//			fprintf(stderr,"statess2 is null in StatesStates\n");
//		if (partials3 < 100)
//			fprintf(stderr,"partials3 is null in StatesStates\n");
//		if (matrices1 < 100)
//			fprintf(stderr,"matrices1 is null in StatesStates\n");
//		if (matrices2 < 100)
//			fprintf(stderr,"matrices2 is null in StatesStates\n");
//#endif	

	v = 0;
	u = 0;
	for (j = 0; j < matrixCount; j++) {

		for (k = 0; k < patternCount; k++) {
		
			w = u;
			
			state1 = states1[k];
			state2 = states2[k];
			
	//		fprintf(stderr,"p = %d:%d\n",state1,state2);
							
			if (state1 < stateCount && state2 < stateCount) {
			
				for (m = 0; m < stateCount; m++) {
					partials3[v] = matrices1[w + state1] * matrices2[w + state2];
	//				fprintf(stderr,"m: %d, p:%d, m:%d, %1.2e = %1.2e * %1.2e\n",j,k,m,
	//					partials3[v], 
	//					matrices1[w+state1],
	//					matrices2[w+state2]);
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
//	fprintf(stderr,"NATIVE: v_last = %d\n",v);

#ifdef DEBUG
	fprintf(stderr,"Leaving SS\n");
#endif

}

