package affy.qualityControl;

public class PsiFunction {
	public static double huber(double u, double k,int deriv){

		if (deriv == 0){
			if ( 1 < k/Math.abs(u)){
				return 1.0;
			} else {
				return  k/Math.abs(u);
			}
		} else if (deriv == 1){
			if (Math.abs(u) <= k){
				return 1.0;
			} else {
				return 0.0;
			}
		} else {
			if (Math.abs(u) <= k){
				return u;
			} else {
				if (u < 0){
					return -k;
				} else {
					return k;
				}
			}
		}
	}

	public double fair(double u, double k,int deriv){

		if (deriv == 0){
			return 1.0/(1.0+Math.abs(u)/k);
		} else if (deriv == 1){
			if (u >=0){
				return 1.0/(1.0+Math.abs(u)/k) - u/(k*(1.0+Math.abs(u)/k)*(1.0+Math.abs(u)/k));
			} else {
				return 1.0/(1.0+Math.abs(u)/k) + u/(k*(1.0+Math.abs(u)/k)*(1.0+Math.abs(u)/k));
			}
		} else {    
			return u/(1.0+Math.abs(u)/k);
		}
	}

	public double cauchy(double u, double k,int deriv){

		if (deriv == 0){
			return 1.0/(1.0+(u/k)*(u/k));
		} else if (deriv == 1){
			return k*k*(k*k - u*u)/((k*k+u*u)*(k*k+u*u));
		} else {    
			return u/(1.0+(u/k)*(u/k));
		}
	}

	public double GemanMcClure(double u, double k,int deriv){

		if (deriv == 0){
			return 1.0/((1.0 + u*u)*(1.0 + u*u));
		} else if (deriv == 1){
			return (1.0 - 3.0*u*u)/((1.0+u*u)*(1.0+u*u)*(1.0+u*u));
		} else {    
			return u/((1.0 + u*u)*(1.0 + u*u));
		}
	}

	public double Welsch(double u, double k,int deriv){

		if (deriv == 0){
			return Math.exp(-(u/k)*(u/k));
		} else if (deriv == 1){
			return Math.exp(-(u/k)*(u/k))*(1 - 2.0*(u*u)/(k*k));
		} else {    
			return u*Math.exp(-(u/k)*(u/k));
		}
	}

	public double Tukey(double u, double k,int deriv){

		if (deriv == 0){
			if (Math.abs(u) <= k){
				return Math.pow((1.0 - (u/k)*(u/k)),2.0);
			} else {
				return 0;
			}
		} else if (deriv == 1){
			if (Math.abs(u) <= k){
				return (1.0 - (u/k)*(u/k))*(1.0-5.0*(u/k)*(u/k));
			} else {
				return 0;
			}
		} else {    
			if (Math.abs(u) <= k){
				return u*(1.0 - (u/k)*(u/k))* (1.0 - (u/k)*(u/k));
			} else {
				return 0;
			}
		}
	}

	public double Andrews(double u, double k,int deriv){

		if (deriv == 0){
			if (Math.abs(u) <= k*Math.PI){
				return Math.sin(u/k)/(u/k);
			} else {
				return 0;
			}
		} else if (deriv == 1){
			if (Math.abs(u) <= k*Math.PI){
				return Math.cos(u/k);
			} else {
				return 0;
			}
		} else {    
			if (Math.abs(u) <= k*Math.PI){
				return k*Math.sin(u/k);
			} else {
				return 0;
			}
		}
	}
}
