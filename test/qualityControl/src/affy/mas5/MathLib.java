package affy.mas5;

import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.List;

public class MathLib {

	final double PI = 3.1415926535897932384626433832795;
	final double TINY = 1e-20;
	final static double Eps = Double.MIN_VALUE;
	final static double RelativeBound = 5 * Eps;
	final static int MaxIterations = 5000;
	final int MeaningLess = -9999;

	double median(final List<Double> v) {
		List<Double> u = new ArrayList<Double> (v);
		int len = u.size();
		if (len < 1) {
			return 0.0;
		}
		Collections.sort(u);
		int half = len / 2;
		if (len % 2 == 1) {
			return u.get(half);
		}
		else {
			return (u.get(half - 1) + u.get(half)) / 2.0;
		}
	}

	DI oneSidedSignRank3(final List<Double>  x, final double alpha) throws IOException 
	{
		//		int len = x.length; // pk unused
		List<Double>  newdiff = new ArrayList<Double> (x);
		int n = 0;
		int i;
		for (i = 0; i < x.size(); ++i) {
			if (x.get(i) != 0.0) {
				newdiff.add(x.get(i));
				n++;
			}
		}
		//newdiff.setSize(n);
		List<Double> ranks = new ArrayList<Double> (n);
		for (i=0; i<n; ++i) {
			ranks.add((double)(i + 1));
		}

		if (n == 0)
			return new DI(0.5, 0);
		else {
			List<DI> absdiff = new ArrayList<DI> (n);
			for (i = 0; i < n; ++i) {
				absdiff.get(i).d = Math.abs(newdiff.get(i));
				absdiff.get(i).i = i;
			}
			Collections.sort(absdiff, new DIComparator());
			int nTies = 0;
			List<Integer> ties = new ArrayList<Integer> (n-1);
			for (i = 0; i < n - 1; ++i) {
				if (absdiff.get(i).d == absdiff.get(i+1).d) {
					ties.add(i);
					nTies ++;
				}
			}
			//ties.setSize(nTies);
			//int tieGroup = 0;
			double doubleVarMod = 0;
			if (nTies!=0) {
				i = 0;
				while ( i < n - 1) {
					double initElement = absdiff.get(i).d;
					int tieGroupSize = 1;
					for (int j = i + 1; j < n; ++j) {
						if (absdiff.get(j).d == initElement) {
							tieGroupSize ++;
							if (j == n - 1) {
								i = j;
								//tieGroup ++;
								for (int m = j - tieGroupSize + 1; m <= j; ++m) {
									ranks.set(m, (2*j - tieGroupSize + 3) / 2.0);
								}
								doubleVarMod += tieGroupSize *
										((double)tieGroupSize * tieGroupSize - 1);
								break;
							}
						}
						else {
							i = j;
							if (tieGroupSize > 1) {
								//tieGroup ++;
								for (int m = j - tieGroupSize; m <= j - 1; ++m) {
									ranks.set(m, (2 * j - tieGroupSize + 1) / 2.0);
								}
								doubleVarMod += tieGroupSize *
										((double)tieGroupSize * tieGroupSize - 1);
							}
							break;
						}
					}
				}
			}
			List<Double> invr = new ArrayList<Double> (n);
			for(int k=0; k<n; k++)
				invr.add(0.0);
			for (i = 0; i < n; ++i) {
				invr.set(absdiff.get(i).i, ranks.get(i));
			}

			double w = 0;
			for (i = 0; i < n; ++i) {
				if (newdiff.get(i) > 0)
					w += invr.get(i);
			}
			DI ans = new DI();
			if (n > 11) {
				double dw = w - ((double) n) * (n + 1) / 4.0;
				double denom2 = (((double)n)*(n+1)*(2*n+1) - 0.5*doubleVarMod)/24.0;
				if (denom2 <=0) {
					System.err.println("denom2="+denom2+" dw="+dw);
					System.exit(1);
				}
				double z = dw / Math.sqrt(denom2);
				ans.d = 1 - normalCDF(z);
			}
			else {
				if (nTies == 0) {
					int myCode = 0;
					for (i = 0; i < n; ++i) {
						if (newdiff.get(i) > 0) {
							myCode += 1 << (invr.get(i).intValue() - 1);
						}
					}
					PTable pTable = new PTable();
					ans.d = pTable.fGetPValue(n-1, myCode);
				}
				else {
					int twoToN = 1 << n;
					List<Integer> mask = new ArrayList<Integer> (n);
					for (i = 0; i < n; ++i) {
						mask.add(1 << i);
					}
					List<Double> posRanks = new ArrayList<Double> (twoToN);
					for (i = 0; i < twoToN; ++i) {
						double sum = 0;
						for (int j = 0; j < n; ++j) {
							if ((i & mask.get(j))!=0)
								sum += ranks.get(j);
						}
						posRanks.add(sum);
					}
					double tail = 0;
					for (i = 0; i < twoToN; ++i) {
						if (posRanks.get(i) > w) {
							tail ++;
						}
						else if (posRanks.get(i) == w) {
							tail += 0.5;
						}
					}
					ans.d = tail / (double) twoToN;
				}
			}
			ans.i = (ans.d < alpha) ? 1 : 0;
			return ans;
		}
	}

	static double normalCDF(final double x) {
		final double sqrt2 = 1.414213562373095048801689;
		double unAdjusted = 0.5 - 0.5 * erf( - x / sqrt2);
		if (unAdjusted > 1) {
			return 1;
		}
		else if (unAdjusted < 0) {
			return 0;
		}
		else {
			return unAdjusted;
		}
	}

	static double erf(final double x) 
	{
		final int erfCoeffSize[] = {3, 7, 4};
		final double erfA1[] = {-3.5609843701815385e-2,
				6.9963834886191355,    2.1979261618294152e1,
				2.4266795523053175e2};
		final double erfB1[] = {1.5082797630407787e1,
				9.1164905404514901e1,  2.1505887586986120e2};
		final double erfA2[] = {-1.368648573827167067e-7,
				5.641955174789739711e-1, 7.211758250883093659,
				4.316222722205673530e1,	1.529892850469404039e2,
				3.393208167343436870e2,	4.519189537118729422e2,
				3.004592610201616005e2};
		final double erfB2[] = {1.278272731962942351e1,
				7.700015293522947295e1,  2.775854447439876434e2,
				6.389802644656311665e2,  9.313540948506096211e2,
				7.909509253278980272e2,  3.004592609569832933e2};
		final double erfA3[] ={2.23192459734184686e-2,
				2.78661308609647788e-1, 2.26956593539686930e-1,
				4.94730910623250734e-2, 2.99610707703542174e-3};
		final double erfB3[] = {1.98733201817135256,
				1.05167510706793207,  1.91308926107829841e-1,
				1.06209230528467918e-2};
		final double erfXBounds[] = {0.46875, 4.0};
		final double invSqrtPi = 0.56418958354775627928; // 1 / sqrt(pi)

		double absX = Math.abs(x);
		double xSquared = absX * absX;
		double temp = 0;


		if (absX <= erfXBounds[0]) {
			double num = erfA1[0] * xSquared;
			double den = xSquared;
			int last = erfCoeffSize[0];
			for (int i = 1; i < last; i++) {
				num = (num + erfA1[i]) * xSquared;
				den = (den + erfB1[i - 1]) * xSquared;
			}
			return x * (num + erfA1[last]) / (den + erfB1[last - 1]);
		}
		else {
			if (absX <= erfXBounds[1]) {
				double num = erfA2[0] * absX;
				double den = absX;
				int last = erfCoeffSize[1];
				for (int i = 1; i < last; i++) {
					num = (num + erfA2[i]) * absX;
					den = (den + erfB2[i - 1]) * absX;
				}
				temp = (num + erfA2[last]) / (den + erfB2[last - 1]);
			}
			else {
				double xInvSquared = 1.0 / xSquared;
				double num = erfA3[0] * xInvSquared;
				double den = xInvSquared;
				int last = erfCoeffSize[2];
				for (int i = 1; i < last; i++) {
					num = (num + erfA3[i]) * xInvSquared;
					den = (den + erfB3[i - 1]) * xInvSquared;
				}
				temp = xInvSquared * (num + erfA3[last]) / (den + erfB3[last - 1]);
				temp = (invSqrtPi - temp) / absX;
			}
			temp = 1 - Math.exp(-xSquared) * temp;
			if (x > 0)      // in fact, we may use if (x > erfXBounds[0])
				return temp;
			else
				return -temp;
		}
	}


	static double cFraction(double x, double a, double b) {
		double ai = 1, bi = 1, y = 1;
		double aPlus1 = a + 1, aPlusB = a + b;
		//		double = aMinus1 = a - 1; // pk unused
		double z0 = 1 - x * aPlusB / aPlus1;
		double error = 100;
		for (int i = 1; i < MaxIterations; ++i) {
			double aPlus2I = a + 2 * i;
			double xModified = x / aPlus2I;
			double c = xModified * i * (b - i) / (aPlus2I - 1);
			double d = - xModified * (a + i) * (aPlusB + i) / (aPlus2I + 1);
			double y1 = y + ai * c;
			double y2 = y1 + y * d;
			double z1 = 0, z2 = 0;
			if (i == 1) {
				z1 = bi * c + z0;
				z2 = z1 + d * z0;
			}
			else {
				z1 = bi * c + 1;
				z2 = z1 + d;
			}
			ai = y1 / z2;
			bi = z1 / z2;
			double yold = y;
			y = y2 / z2;
			error = Math.abs(y - yold);
			if (error < RelativeBound * Math.abs(y)) {
				return y;
			}
		}
		return y;
	}

	static double incompleteBeta(final double x, final double a, final double b)  {
		if (x == 0) {
			return 0;
		}
		else if (x == 1) {
			return 1;
		}
		else {
			double coeff = Math.pow(x, a) * Math.pow(1 - x, b) /
					Math.exp(logGamma(a) + logGamma(b) - logGamma(a + b));
			if (x < (1.0 + a) / (2.0 + a + b)) { 
				return coeff * cFraction(x, a, b) / a;
			}
			else {
				return 1.0 - coeff * cFraction(1.0 - x, b, a) / b;
			}
		}
	}
	static double tCDF(double t, double df) {
		double p;
		if (df == 1) {
			p = 0.5 + Math.atan(t) / Math.PI;
		}
		else {
			if (t == 0) {
				p = 0.5;
			}
			else {
				double x = df / (t*t + df);
				double y = 1 - incompleteBeta(x, df/2.0, 0.5);
				p = (1 - y) / 2.0;
				if (t > 0) {
					p = 1 - p;
				}
			}
		}
		if (p > 1)
			return 1;
		else if (p < 0)
			return 0;
		else
			return p;
	}

	static double tCDFinversed(double p, double df) {
		double t = 0;
		if (p == 0) {
			t = -9999.0;
		}
		else if (p == 1) {
			t = 9999.0;
		}
		else if (p == 0.5) {
			t = 0;
		}
		else if (df == 1) {
			t = Math.tan(Math.PI * (p - 0.5));
		}	
		else if (p < 0.5) {
			t = - tCDFinversed(1 - p, df);
		}
		else {
			int maxIter = 120;
			double rB = 0.00001;
			double aB = 0.00005;
			//			double coeff = exp(logGamma((df + 1) / 2.0) - logGamma(df / 2.0))
			//					/ sqrt(df * PI); // pk unused
			//			bool success = false; // pk unused
			double t0 = 1;
			double t1 = 2;
			double p0 = tCDF(t0, df);
			if (p0 == p) {
				return t0;
			}
			double p1 = tCDF(t1, df);
			if (p1 == p) {
				return t1;
			}
			for (int i = 0; i < maxIter; ++i) {
				if (p1 == p) {
					return t1;
				}
				else if (p0 == p) {
					return t0;
				}
				else if (p0 > p) {
					t1 = t0;
					p1 = p0;
					t0 = t0 / 2.0;
					p0 = tCDF(t0, df);
				}
				else if (p1 < p) {
					t0 = t1;
					p0 = p1;
					t1 = 2 * t1;
					p1 = tCDF(t1, df);
				}
				else {
					double midT = (t0 + t1) / 2.0;
					double midP = tCDF(midT, df);
					double dT = Math.abs(t1 - midT);
					if (midP == p || dT < aB || dT < rB * t1) {
						return midT;
					}
					else if (midP < p) {
						t0 = midT;
						p0 = tCDF(t0, df);
					}
					else {
						t1 = midT;
						p1 = tCDF(t1, df);
					}
				}
			}
		}
		return t;
	}

	static double logGamma(final double x) {
		double a[] = {5.7083835261e-03, -1.910444077728e-03,
				8.4171387781295e-04, -5.952379913043012e-04,
				7.93650793500350248e-04, -2.777777777777681622553e-03,
				8.333333333333333331554247e-02, 0.9189385332046727417803297};

		double c1[] = {4.945235359296727046734888e0,
				2.018112620856775083915565e2, 2.290838373831346393026739e3,
				1.131967205903380828685045e4, 2.855724635671635335736389e4,
				3.848496228443793359990269e4, 2.637748787624195437963534e4,
				7.225813979700288197698961e3, -5.772156649015328605195174e-1};

		double d1[] = {6.748212550303777196073036e1,
				1.113332393857199323513008e3, 7.738757056935398733233834e3,
				2.763987074403340708898585e4, 5.499310206226157329794414e4,
				6.161122180066002127833352e4, 3.635127591501940507276287e4,
				8.785536302431013170870835e3};

		double c2[] = {4.974607845568932035012064e0,
				5.424138599891070494101986e2, 1.550693864978364947665077e4,
				1.847932904445632425417223e5, 1.088204769468828767498470e6,
				3.338152967987029735917223e6, 5.106661678927352456275255e6,
				3.074109054850539556250927e6, 4.227843350984671393993777e-1};

		double d2[] = {1.830328399370592604055942e2,
				7.765049321445005871323047e3, 1.331903827966074194402448e5,
				1.136705821321969608938755e6, 5.267964117437946917577538e6,
				1.346701454311101692290052e7, 1.782736530353274213975932e7,
				9.533095591844353613395747e6};

		double c4[] = {-1.474502166059939948905062e4,
				-2.426813369486704502836312e6, -1.214755574045093227939592e8,
				-2.663432449630976949898078e9, -2.940378956634553899906876e10,
				-1.702665737765398868392998e11, -4.926125793377430887588120e11,
				-5.606251856223951465078242e11, 1.791759469228055000094023e0};

		double d4[] = {-2.690530175870899333379843e3,
				-6.393885654300092398984238e5, -4.135599930241388052042842e7,
				-1.120872109616147941376570e9, -1.488613728678813811542398e10,
				-1.016803586272438228077304e11, -3.417476345507377132798597e11,
				-4.463158187419713286462081e11};

		double numerator = 0;
		double denominator = 1;
		int n = 8;

		if (x <= Eps) {
			return - Math.log(x);
		}
		else if (x <= 0.5) {
			for (int i = 0; i < n; ++i) {
				numerator = x * numerator + c1[i];
				denominator = x * denominator + d1[i];
			}
			return - Math.log(x) + x * (x * numerator / denominator + c1[n]);
		}
		else if (x <= 0.6796875) {
			double z = x - 1;
			for (int i = 0; i < n; ++i) {
				numerator = z * numerator + c2[i];
				denominator = z * denominator + d2[i];
			}
			return -Math.log(x) + z * (z * numerator / denominator + c2[n]);
		}
		else if (x <= 1.5) {
			double z = x - 1;
			for (int i = 0; i < n; ++i) {
				numerator = z * numerator + c1[i];
				denominator = z * denominator + d1[i];
			}
			return z * (z * numerator / denominator + c1[n]);
		}
		else if (x <= 4) {
			double z = x - 2;
			for (int i = 0; i < n; ++i) {
				numerator = z * numerator + c2[i];
				denominator = z * denominator + d2[i];
			}
			return z * (z * numerator / denominator + c2[n]);
		}
		else if (x <= 12) {
			double z = x - 4;
			for (int i = 0; i < n; ++i) {
				numerator = z * numerator + c4[i];
				denominator = z * denominator + d4[i];
			}
			return z * numerator / denominator + c4[n];
		}
		else {
			double z = x * x;
			numerator = a[0];
			for (int i = 1; i <= 6; ++i) {
				numerator = numerator / z + a[i];
			}
			numerator /= x;
			return numerator + Math.log(x) * (x - 0.5) - x + a[7];
		}	
	}

	double t_distribution_lookup(final int df) {

		final double dfFrom1To30[] = {
				6.314, 2.920, 2.353, 2.132, 2.015,
				1.943, 1.895, 1.860, 1.833, 1.812,
				1.796, 1.782, 1.771, 1.761, 1.753,
				1.746, 1.740, 1.734, 1.729, 1.725,
				1.721, 1.717, 1.714, 1.711,	1.708,
				1.706, 1.703, 1.701, 1.699, 1.697
		};

		final double df40 = 1.684;
		final double df50 = 1.676;
		final double df75 = 1.665;
		final double df100 = 1.660;
		final double df200 = 1.653;
		final double df1000 = 1.646;

		if (df >=1 && df <= 30)
			return dfFrom1To30[df - 1];
		else if (df > 30 && df <= 40)
			return df40;
		else if (df > 40 && df <= 50)
			return df50;
		else if (df > 50 && df <= 75)
			return df75;
		else if (df > 75 && df <= 100)
			return df100;
		else if (df > 100 && df <= 200)
			return df200;
		else
			return df1000;
	}

}

class DI {
	public double d;
	public int i;
	DI(double dv, int iv)
	{
		d = dv;
		i = iv;
	}
	DI()
	{
		d = 0.0;
		i = 0;
	}
};

class DIComparator implements Comparator<DI>{

	public int compare(DI o1, DI o2) {
		return o1.d > o2.d? 1 : -1;
	}
}
