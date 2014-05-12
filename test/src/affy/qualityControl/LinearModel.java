package affy.qualityControl;

import org.apache.commons.math3.distribution.TDistribution;

import java.util.List;


public class LinearModel {
	private double intercept;
	private double slope;
	private double slopeStdErr;
	private double interceptStdErr;
	private double interceptProT;
	private double slopeProT;
	
	public double getIntercept(){
		return intercept;
	}
	
	public double getSlope(){
		return slope;
	}
	
	public double getSlopeStdErr(){
		return slopeStdErr;
	}
	
	public double getInterceptStdErr(){
		return interceptStdErr;
	}
	
	public double getInterceptProT(){
		return interceptProT;
	}
	
	public double getSlopeProT(){
		return slopeProT;
	}
	
	public void linearRegression(List<Double> x, List<Double> y){
		double sumx = 0.0;
		double sumy = 0.0;
		//double sumx2 = 0.0;
		int n = x.size();

		// first pass: compute xbar and ybar
		for (int i=0; i<n; i++){
			sumx += x.get(i);
			//sumx2 += x.get(i)*x.get(i);
			sumy += y.get(i);
		}
		double xbar = sumx / n;
		double ybar = sumy / n;

		// second pass: compute summary statistics
		double xxbar = 0.0;
		//double yybar = 0.0;
		double xybar = 0.0;
		
		for (int i = 0; i < n; i++) {
			xxbar += (x.get(i) - xbar) * (x.get(i) - xbar);
			//yybar += (y.get(i) - ybar) * (y.get(i) - ybar);
			xybar += (x.get(i) - xbar) * (y.get(i) - ybar);
		}
		//double beta1 = xybar / xxbar;
		//double beta0 = ybar - beta1 * xbar;
		
		slope = xybar / xxbar;
		intercept = ybar - slope * xbar;

		// print results
		//System.out.println("y   = " + beta1 + " * x + " + beta0);
		
        // analyze results
        int df = n - 2;        // degree of freedom
        double rss = 0.0;      // residual sum of squares
        //double ssr = 0.0;      // regression sum of squares
        for (int i = 0; i < n; i++) {
            double fit = slope*x.get(i) + intercept;
            rss += (fit - y.get(i)) * (fit - y.get(i));
           // ssr += (fit - ybar) * (fit - ybar);
        }
        //double R2    = ssr / yybar;
        double svar  = rss / df;
        double svar1 = svar / xxbar;
        double svar0 = svar/n + xbar*xbar*svar1;
        
        slopeStdErr = Math.sqrt(svar1);
        interceptStdErr = Math.sqrt(svar0);
        
/*      System.out.println("R^2                 = " + R2);
        System.out.println("std error of beta_1 = " + Math.sqrt(svar1));
        System.out.println("std error of beta_0 = " + Math.sqrt(svar0));
        svar0 = svar * sumx2 / (n * xxbar);
        System.out.println("std error of beta_0 = " + Math.sqrt(svar0));

        System.out.println("SSTO = " + yybar);
        System.out.println("SSE  = " + rss);
        System.out.println("SSR  = " + ssr);*/
       
        
        TDistribution tDistribution = new TDistribution(df);
        slopeProT = 2*(1-tDistribution.cumulativeProbability(Math.abs(slope/slopeStdErr)));
        interceptProT = 2*(1-tDistribution.cumulativeProbability(Math.abs(intercept/interceptStdErr)));
	}
	
	
}
