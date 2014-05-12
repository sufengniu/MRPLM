#include <jni.h>
#include <stdio.h>
#include <time.h>

#include "jniWrapper.h"


JNIEXPORT jdoubleArray JNICALL Java_jniWrapper_wlsAcc
(JNIEnv *env, jobject thisObj, jdoubleArray inJNIArray, jint y_rows, jint y_cols) {
	// Step 1: Convert the incoming JNI jintarray to C's jint[]
	jint *inCArray = (*env)->GetIntArrayElements(env, inJNIArray, NULL);
	if (NULL == inCArray) return NULL;
	jsize length = (*env)->GetArrayLength(env, inJNIArray);

	
	(*env)->ReleaseIntArrayElements(env, inJNIArray, inCArray, 0); // release resources

	jdouble outCArray[];
	
	
	// Step 3: Convert the C's Native jdouble[] to JNI jdoublearray, and return
	jdoubleArray outJNIArray = (*env)->NewDoubleArray(env, y_rows+y_cols-1);  // allocate
	if (NULL == outJNIArray) return NULL;
	(*env)->SetDoubleArrayRegion(env, outJNIArray, 0 , y_rows+y_cols-1, outCArray);  // copy
	return outJNIArray;

}

