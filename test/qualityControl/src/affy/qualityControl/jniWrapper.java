package affy.qualityControl;

public class jniWrapper{
	
	static {
		System.loadLibrary("jniWrapper");
	}
	
	// application callback variables
	
	/* list accelerated functions by using JNI */
	public native void wlsAcc(double[] weights, double[] y, double[] out_beta, int y_rows, int y_cols);
	public native double[] seAcc(double[] weights, int y_rows, int y_cols);

}
