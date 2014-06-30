package affy.qualityControl;

public class jniWrapper{
	
	static {
		System.loadLibrary("jniWrapper");
	}
	
	// application callback variables
	
	/* list accelerated functions by using JNI */
	public native void rlm_fit_anova(double[] y, int y_rows, int y_cols, double[] out_beta, double[] resids, double[] weights, double psi_k, int max_iter);
	public native void rlm_compute_se_anova(int y_rows, int y_cols, double[] resids, double[] weights, double[] se_estimates, double[] varcov, double psi_k);
	
	
	//public native void wlsAcc(double[] weights, double[] y, double[] out_beta, int y_rows, int y_cols);
	//public native double[] seAcc(double[] weights, double[] xtx, int y_rows, int y_cols);

}
