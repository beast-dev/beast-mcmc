#include "SparseMatrixExponential.h"
#include <jni.h>


extern "C" void dgexpvms_(long *n, long *nz, long *m, double *t, double *v, double *w, double *tol, double *anorm, long *ia, long *ja, double *a, double *wsp, long *lwsp, long *iwsp, long *liwsp, long *itrace, long *iflag);

JNIEXPORT void JNICALL Java_dr_math_SparseMatrixExponential_executeDGEXPV
(JNIEnv *env, jobject, jint order, jint nz, jint maxBasis, jdouble time, jdoubleArray operandVector, jdoubleArray resultVector, jdouble tolerance, jdouble matrixNorm, jdoubleArray rate, jintArray indexX, jintArray indexY, jdoubleArray workspace, jint lengthWorkspace, jintArray intWorkspace, jint lengthIntWorkspace, jint flag) {




//  jdouble* operandVectorPtr = env->GetDoubleArrayElements(operandVector,0);
//  jdouble* resultVectorPtr = env->GetDoubleArrayElements(resultVector,0);
//  jdouble* ratePtr = env->GetDoubleArrayElements(rate,0);
//  jdouble* workspacePtr = env->GetDoubleArrayElements(workspace,0);
//  jint* indexXPtr = env->GetIntArrayElements(indexX,0);
//  jint* indexYPtr = env->GetIntArrayElements(indexY,0);
//  jint* intWorkspacePtr = env->GetIntArrayElements(intWorkspace,0);

  jdouble* operandVectorPtr = (jdouble*)env->GetPrimitiveArrayCritical(operandVector,0);
  jdouble* resultVectorPtr = (jdouble*)env->GetPrimitiveArrayCritical(resultVector,0);
  jdouble* ratePtr = (jdouble*)env->GetPrimitiveArrayCritical(rate,0);
  jdouble* workspacePtr = (jdouble*)env->GetPrimitiveArrayCritical(workspace,0);
  jint* indexXPtr = (jint*)env->GetPrimitiveArrayCritical(indexX,0);
  jint* indexYPtr = (jint*)env->GetPrimitiveArrayCritical(indexY,0);
  jint* intWorkspacePtr = (jint*)env->GetPrimitiveArrayCritical(intWorkspace,0);

  long itrace = 0;

    dgexpvms_(&order, &nz, &maxBasis, &time, operandVectorPtr, resultVectorPtr,
  	   &tolerance, &matrixNorm, indexXPtr, indexYPtr, ratePtr,
  	   workspacePtr, &lengthWorkspace, intWorkspacePtr, 
  	   &lengthIntWorkspace, &itrace, &flag);

  env->ReleasePrimitiveArrayCritical(operandVector,operandVectorPtr,0);
  env->ReleasePrimitiveArrayCritical(resultVector,resultVectorPtr,0);
  env->ReleasePrimitiveArrayCritical(rate,ratePtr,0);
  env->ReleasePrimitiveArrayCritical(workspace,workspacePtr,0);
  env->ReleasePrimitiveArrayCritical(indexX,indexXPtr,0);
  env->ReleasePrimitiveArrayCritical(indexY,indexYPtr,0);
  env->ReleasePrimitiveArrayCritical(intWorkspace,intWorkspacePtr,0);

//  env->ReleaseDoubleArrayElements(operandVector,operandVectorPtr,0);
//  env->ReleaseDoubleArrayElements(resultVector,resultVectorPtr,0);
//  env->ReleaseDoubleArrayElements(rate,ratePtr,0);
//  env->ReleaseDoubleArrayElements(workspace,workspacePtr,0);
//  env->ReleaseIntArrayElements(indexX,indexXPtr,0);
//  env->ReleaseIntArrayElements(indexY,indexYPtr,0);
//  env->ReleaseIntArrayElements(intWorkspace,intWorkspacePtr,0);


}
