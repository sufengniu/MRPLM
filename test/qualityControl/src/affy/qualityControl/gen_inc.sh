#!/bin/bash

# first make sure you already set up the $JAVA_HOME for java
javac -d ../../ ./jniWrapper.java

javah -jni -classpath ../../ affy.qualityControl.jniWrapper

mv affy_qualityControl_jniWrapper.h ../../plmAcc/
rm *.class
