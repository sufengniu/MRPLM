#ifndef JNIWRAPPER_H_
#define JNIWRAPPER_H_

JNIEXPORT void JNICALL Java_jniWrapper_wlsAcc(JNIEnv *env, jobject thisObj, jdoubleArray inJNIWeights, jdoubleArray outJNIY, jdoubleArray outJNIBeta, jint y_rows, jint y_cols);
JNIEXPORT jdoubleArray JNICALL Java_jniWrapper_seAcc(JNIEnv *env, jobject thisObj, jdoubleArray inJNIArray, jint y_rows, jint y_cols);

#endif
