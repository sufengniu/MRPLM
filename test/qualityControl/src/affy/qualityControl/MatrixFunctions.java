package affy.qualityControl;

public class MatrixFunctions {
	private final static boolean use_lapack = false;

	private static boolean Choleski_decompose(double []X, double []L, int n, boolean lapack){
		boolean error_code = false;
		// char upper = 'U';

		for (int i=0; i < n; i++){
			for (int j=0; j < n; j++){
				if (i > j)
					L[j*n+i] = 0.0;
				else {
					L[j*n+i] = X[j*n + i];
				}
			}
		}
		//if (!lapack){
		return Linpack.dpofa(L,n,n)==1 ? true : false;
		/*		  } else {
		    //dpotrf_(&upper,&n,L,&n,&error_code);
		  }*/

		//return error_code;
	}

	private static boolean Choleski_2_inverse(double [] X, double [] Xinv, int n,boolean upperonly, boolean lapack){

		boolean error_code=false;
		double [] d = new double [2];
		//char upper = 'U';

		for (int i=0; i < n; i++){ 
			/* check for a zero or close to zero diagonal element */ 
			if(Math.abs(X[i*n+ i]) < 1e-06){
				error_code = true;
				return error_code;
			}

			for (int j=i; j < n; j++){
				Xinv[j*n + i] = X[j*n + i];
			}
		}

		int inverseonly = 1;
		if (!lapack){
			Linpack.dpodi(Xinv,n,n,d,inverseonly);
		} else {
			//dpotri_(&upper,&n,Xinv,&n,&error_code);
		}
		if (!upperonly){
			for (int i=0; i < n; i++){
				for (int j=0; j <= i-1; j++){
					Xinv[j*n+i] = Xinv[i*n+j];
				}
			}
		}
		return error_code;
	}

	public static boolean Choleski_inverse(double [] X, double [] Xinv, double [] work, int n, boolean upperonly){

		boolean error_code = Choleski_decompose(X, work, n,use_lapack);
		if (!error_code){
			error_code = Choleski_2_inverse(work, Xinv, n,upperonly,use_lapack);
		}
		return error_code;
	}

	private boolean SVD_compute(double [] X, int n, double [] s, double [] u, double [] v,boolean lapack){

		int lwork = 7*n*n + 4*n;
		int job = 21;
		//boolean error_code = false;
		//char jobz = 'A';
		double []Xcopy= new double [n*n];              /* Calloc(n*n,double); */
		double []e =    new double [n];                /* Calloc(n,double); */
		double []work =  new double [n];               /* Calloc(n,double); */
		//double []work2 =  new double [lwork];
		//int []iwork = new int [8*n];


		for (int i=0; i < n; i++){
			for (int j=0; j < n; j++){
				Xcopy[j*n + i] = X[j*n+i];
			}
		}
		// if (!lapack){
		return Linpack.dsvdc(Xcopy,n,n,n,s,e,u,n,v,n,work,job)==1 ? true: false;
		/*		  } else {
		    //dgesdd_(&jobz,&n,&n,Xcopy,&n,s,u,&n,v,&n,work2,&lwork,iwork,&error_code);
		  }*/


	}

	private int SVD_2_inverse(double [] Xinv, int n, double [] s, double [] u, double []v, boolean lapack){
		double tolerance = 1e-7; /* 1.490116e-08; */
		int nonzero = n;


		for (int i = 0; i < n; i++){
			if (s[i] < tolerance*s[0]){
				nonzero = i;
				/* printf("nonzero %d",i); */
				break;
			}
		}


		/* for all columns where $d is not to small do */
		/*  svdu$v %*% (t(svdu$u)* 1/svdu$d); */

		for (int i = 0; i < n; i++){
			for (int j = 0; j < nonzero; j++){
				u[j*n + i] = u[j*n+i] * 1.0/s[j];
			}
		}
		if (!lapack){
			for (int i = 0; i < n; i++){
				for (int j = 0; j < n; j++){
					Xinv[j*n+i] = 0.0;
					for (int k=0; k < nonzero; k++){
						Xinv[j*n+i]+= v[k*n+i] * u[k*n+j];
					}
				}
			}
		} else {
			/* lapack so v is transposed */
			for (int i = 0; i < n; i++){
				for (int j = 0; j < n; j++){
					Xinv[j*n+i] = 0.0;
					for (int k=0; k < nonzero; k++){
						Xinv[j*n+i]+= v[i*n+k] * u[k*n+j];
					}
				}
			}
		}
		return 0;
	}

	private boolean SVD_inverse(double [] X, double [] Xinv, int n){

		boolean error_code = false;
		double [] s = new double [n+1];
		double [] v = new double [n*n];
		double [] u = new double [n*n];

		error_code = SVD_compute(X, n, s, u, v,use_lapack);
		SVD_2_inverse(Xinv,n, s, u, v,use_lapack);

		return error_code;
	}
}
