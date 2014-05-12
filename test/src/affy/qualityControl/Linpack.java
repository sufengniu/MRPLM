package affy.qualityControl;

public class Linpack {

	//****************************************************************************80

	int dchdc ( double a[], int lda, int p, double work[], int ipvt[], int job )
	{
		int info;
		int j;
		int jp;
		int jt;
		int k;
		int l;
		double maxdia;
		int maxl;
		boolean negk;
		int pl;
		int pu;
		boolean swapk;
		double temp;

		pl = 1;
		pu = 0;
		info = p;
		//
		//  Pivoting has been requested.
		//  Rearrange the the elements according to IPVT.
		//
		if ( job != 0 )
		{
			for ( k = 1; k <= p; k++ )
			{
				swapk = ( 0 < ipvt[k-1] );
				negk = ( ipvt[k-1] < 0 );

				if ( negk )
				{
					ipvt[k-1] = -k;
				}
				else
				{
					ipvt[k-1] = k;
				}

				if ( swapk )
				{
					if ( k != pl )
					{
						Blas1.dswap  ( pl-1, a, 0+(k-1)*lda, 1, a, 0+(pl-1)*lda, 1 );

						temp = a[k-1+(k-1)*lda];
						a[k-1+(k-1)*lda] = a[pl-1+(pl-1)*lda];
						a[pl-1+(pl-1)*lda] = temp;

						for ( j = pl+1; j <= p; j++ )
						{
							if ( j < k )
							{
								temp = a[pl-1+(j-1)*lda];
								a[pl-1+(j-1)*lda] = a[j-1+(k-1)*lda];
								a[j-1+(k-1)*lda] = temp;
							}
							else if ( k < j )
							{
								temp = a[k-1+(j-1)*lda];
								a[k-1+(j-1)*lda] = a[pl-1+(j-1)*lda];
								a[pl-1+(j-1)*lda] = temp;
							}
						}
						ipvt[k-1] = ipvt[pl-1];
						ipvt[pl-1] = k;
					}
					pl = pl + 1;
				}
			}
			pu = p;

			for ( k = p; pl <= k; k-- )
			{
				if ( ipvt[k-1] < 0 )
				{
					ipvt[k-1] = -ipvt[k-1];

					if ( pu != k )
					{
						Blas1.dswap  ( k-1, a, 0+(k-1)*lda, 1, a, 0+(pu-1)*lda, 1 );

						temp = a[k-1+(k-1)*lda];
						a[k-1+(k-1)*lda] = a[pu-1+(pu-1)*lda];
						a[pu-1+(pu-1)*lda] = temp;

						for ( j = k+1; j <= p; j++ )
						{
							if ( j < pu )
							{
								temp = a[k-1+(j-1)*lda];
								a[k-1+(j-1)*lda] = a[j-1+(pu-1)*lda];
								a[j-1+(pu-1)*lda] = temp;
							}
							else if ( pu < j )
							{
								temp = a[k-1+(j-1)*lda];
								a[k-1+(j-1)*lda] = a[pu-1+(j-1)*lda];
								a[pu-1+(j-1)*lda] = temp;
							}
						}
						jt = ipvt[k-1];
						ipvt[k-1] = ipvt[pu-1];
						ipvt[pu-1] = jt;
					}
					pu = pu - 1;
				}
			}
		}

		for ( k = 1; k <= p; k++ )
		{
			//
			//  Reduction loop.
			//
			maxdia = a[k-1+(k-1)*lda];
			maxl = k;
			//
			//  Determine the pivot element.
			//
			if ( pl <= k && k < pu )
			{
				for ( l = k+1; l <= pu; l++ )
				{
					if ( maxdia < a[l-1+(l-1)*lda] )
					{
						maxdia = a[l-1+(l-1)*lda];
						maxl = l;
					}
				}
			}
			//
			//  Quit if the pivot element is not positive.
			//
			if ( maxdia <= 0.0 )
			{
				info = k - 1;
				return info;
			}
			//
			//  Start the pivoting and update IPVT.
			//
			if ( k != maxl )
			{
				Blas1.dswap  ( k-1, a, 0+(k-1)*lda, 1, a, 0+(maxl-1)*lda, 1 );
				a[maxl-1+(maxl-1)*lda] = a[k-1+(k-1)*lda];
				a[k-1+(k-1)*lda] = maxdia;
				jp = ipvt[maxl-1];
				ipvt[maxl-1] = ipvt[k-1];
				ipvt[k-1] = jp;
			}
			//
			//  Reduction step.
			//  Pivoting is contained across the rows.
			//
			work[k-1] = Math.sqrt ( a[k-1+(k-1)*lda] );
			a[k-1+(k-1)*lda] = work[k-1];

			for ( j = k+1; j <= p; j++ )
			{
				if ( k != maxl )
				{
					if ( j < maxl )
					{
						temp = a[k-1+(j-1)*lda];
						a[k-1+(j-1)*lda] = a[j-1+(maxl-1)*lda];
						a[j-1+(maxl-1)*lda] = temp;
					}
					else if ( maxl < j )
					{
						temp = a[k-1+(j-1)*lda];
						a[k-1+(j-1)*lda] = a[maxl-1+(j-1)*lda];
						a[maxl-1+(j-1)*lda] = temp;
					}
				}
				a[k-1+(j-1)*lda] = a[k-1+(j-1)*lda] / work[k-1];
				work[j-1] = a[k-1+(j-1)*lda];
				temp = -a[k-1+(j-1)*lda];
				Blas1.daxpy ( j-k, temp, work, k, 1, a, k+(j-1)*lda, 1 );
			}
		}

		return info;
	}
	//****************************************************************************80

	int dchdd ( double r[], int ldr, int p, double x[], double z[], int ldz, 
			int nz, double y[], double rho[], double c[], double s[] )
	{
		double a;
		double alpha;
		double azeta;
		double b;
		int i;
		int ii;
		int info;
		int j;
		double norm;
		double scale;
		double t;
		double xx;
		double zeta;
		//
		//  Solve R' * A = X, placing the result in the array S.
		//
		info = 0;
		s[0] = x[0] / r[0+0*ldr];

		for ( j = 2; j <= p; j++ )
		{
			s[j-1] = x[j-1] - Blas1.ddot ( j-1, r, 0+(j-1)*ldr, 1, s, 0, 1 );
			s[j-1] = s[j-1] / r[j-1+(j-1)*ldr];
		}

		norm = Blas1.dnrm2 ( p, s, 0, 1 );

		if ( 1.0 <= norm )
		{
			info = -1;
			return info;
		}

		alpha = Math.sqrt ( 1.0 - norm * norm );
		//
		//  Determine the transformations.
		//
		for ( ii = 1; ii <= p; ii++ )
		{
			i = p - ii + 1;
			scale = alpha + Blas1.r8_abs ( s[i-1] );
			a = alpha / scale;
			b = s[i-1] / scale;
			norm = Math.sqrt ( a * a + b * b );
			c[i-1] = a / norm;
			s[i-1] = b / norm;
			alpha = scale * norm;
		}
		//
		//  Apply the transformations to R.
		//
		for ( j = 1; j <= p; j++ )
		{
			xx = 0.0;
			for ( ii = 1; ii <= j; ii++ )
			{
				i = j - ii + 1;
				t = c[i-1] * xx + s[i-1] * r[i-1+(j-1)*ldr];
				r[i-1+(j-1)*ldr] = c[i-1] * r[i-1+(j-1)*ldr] - s[i-1] * xx;
				xx = t;
			}
		}
		//
		//  If required, downdate Z and RHO.
		//
		for ( j = 1; j <= nz; j++ )
		{
			zeta = y[j-1];
			for ( i = 1; i <= p; i++ )
			{
				z[i-1+(j-1)*ldz] = ( z[i-1+(j-1)*ldz] - s[i-1] * zeta ) / c[i-1];
				zeta = c[i-1] * zeta - s[i-1] * z[i-1+(j-1)*ldz];
			}

			azeta = Blas1.r8_abs ( zeta );

			if ( rho[j-1] < azeta )
			{
				info = 1;
				rho[j-1] = -1.0;
			}
			else
			{
				rho[j-1] = rho[j-1] * Math.sqrt ( 1.0 - Math.pow ( azeta / rho[j-1], 2 ) );
			}
		}

		return info;
	}
	//****************************************************************************80

	void dchex ( double r[], int ldr, int p, int k, int l, double z[], int ldz, 
			int nz, double c[], double s[], int job )
	{
		int i;
		int ii;
		int il;
		int iu;
		int j;
		int jj;
		int lm1;
		int lmk;
		double t;
		//
		//  Initialize
		//
		lmk = l - k;
		lm1 = l - 1;
		//
		//  Right circular shift.
		//
		if ( job == 1 )
		{
			//
			//  Reorder the columns.
			//
			for ( i = 1; i <= l; i++ )
			{
				ii = l - i + 1;
				s[i-1] = r[ii-1+(l-1)*ldr];
			}

			for ( jj = k; jj <= lm1; jj++ )
			{
				j = lm1 - jj + k;
				for ( i = 1; i <= j; i++ )
				{
					r[i-1+(j)*ldr] = r[i-1+(j-1)*ldr];
				}
				r[j+(j)*ldr] = 0.0;
			}

			for ( i = 1; i <= k-1; i++ )
			{
				ii = l - i + 1;
				r[i-1+(k-1)*ldr] = s[ii-1];
			}
			//
			//  Calculate the rotations.
			//
			t = s[0];
			for ( i = 1; i <= lmk; i++ )
			{
				//////////////////////////////////
				double [] drotg = new double [4];
				drotg[0] = s[i];
				drotg[1] = t;
				drotg[2] = c[i-1];
				drotg[3] = s[i-1];
				///////////////////////////////////
				Blas1.drotg (drotg);
				///////////////////////////////////
				s[i] = drotg[0];
				t = drotg[1];
				c[i-1] = drotg[2];
				s[i-1] = drotg[3];
				///////////////////////////////////
				
				t = s[i];
			}

			r[k-1+(k-1)*ldr] = t;

			for ( j = k+1; j <= p; j++ )
			{
				il = Blas1.i4_max ( 1, l-j+1 );
				for ( ii = il; ii <= lmk; ii++ )
				{
					i = l - ii;
					t = c[ii-1] * r[i-1+(j-1)*ldr] + s[ii-1] * r[i+(j-1)*ldr];
					r[i+(j-1)*ldr] = c[ii-1] * r[i+(j-1)*ldr] - s[ii-1] * r[i-1+(j-1)*ldr];
					r[i-1+(j-1)*ldr] = t;
				}
			}
			//
			//  If required, apply the transformations to Z.
			//
			for ( j = 1; j <= nz; j++ )
			{
				for ( ii = 1; ii <= lmk; ii++ )
				{
					i = l - ii;
					t = c[ii-1] * z[i-1+(j-1)*ldr] + s[ii-1] * z[i+(j-1)*ldr];
					z[i+(j-1)*ldr] = c[ii-1] * z[i+(j-1)*ldr] - s[ii-1] * z[i-1+(j-1)*ldr];
					z[i-1+(j-1)*ldr] = t;
				}
			}
		}
		//
		//  Left circular shift.
		//
		else
		{
			//
			//  Reorder the columns.
			//
			for ( i = 1; i <= k; i++ )
			{
				ii = lmk + i;
				s[ii-1] = r[i-1+(k-1)*ldr];
			}

			for ( j = k; j <= lm1; j++ )
			{
				for ( i = 1; i <= j; i++ )
				{
					r[i-1+(j-1)*ldr] = r[i-1+(j)*ldr];
				}
				jj = j - k + 1;
				s[jj-1] = r[j+(j)*ldr];
			}

			for ( i = 1; i <= k; i++ )
			{
				ii = lmk + i;
				r[i-1+(l-1)*ldr] = s[ii-1];
			}

			for ( i = k+1; i <= l; i++ )
			{
				r[i-1+(l-1)*ldr] = 0.0;
			}
			//
			//  Reduction loop.
			//
			for ( j = k; j <= p; j++ )
			{
				//
				//  Apply the rotations.
				//
				if ( j != k )
				{
					iu = Blas1.i4_min ( j-1, l-1 );

					for ( i = k; i <= iu; i++ )
					{
						ii = i - k + 1;
						t = c[ii-1] * r[i-1+(j-1)*ldr] + s[ii-1] * r[i+(j-1)*ldr];
						r[i+(j-1)*ldr] = c[ii-1] * r[i+(j-1)*ldr] 
								- s[ii-1] * r[i-1+(j-1)*ldr];
						r[i-1+(j-1)*ldr] = t;
					}
				}

				if ( j < l )
				{
					jj = j - k + 1;
					t = s[jj-1];
					//Blas1.drotg ( r, j-1+(j-1)*ldr, &t, c+jj-1, s+jj-1 );
					//////////////////////////////////
					double [] drotg = new double [4];
					drotg[0] = r[j-1+(j-1)*ldr];
					drotg[1] = t;
					drotg[2] = c[jj-1];
					drotg[3] = s[jj-1];
					///////////////////////////////////
					Blas1.drotg (drotg);
					///////////////////////////////////
					r[j-1+(j-1)*ldr] = drotg[0];
					t = drotg[1];
					c[jj-1] = drotg[2];
					s[jj-1] = drotg[3];
					///////////////////////////////////
				}
			}
			//
			//  Apply the rotations to Z.
			//
			for ( j = 1; j <= nz; j++ )
			{
				for ( i = k; i <= lm1; i++ )
				{
					ii = i - k + 1;
					t = c[ii-1] * z[i-1+(j-1)*ldr] + s[ii-1] * z[i+(j-1)*ldr];
					z[i+(j-1)*ldr] = c[ii-1] * z[i+(j-1)*ldr] - s[ii-1] * z[i-1+(j-1)*ldr];
					z[i-1+(j-1)*ldr] = t;
				}
			}
		}

		return;
	}
	//****************************************************************************80

	void dchud ( double r[], int ldr, int p, double x[], double z[], int ldz, 
			int nz, double y[], double rho[], double c[], double s[] )
	{
		double azeta;
		int i;
		int j;
		double scale;
		double t;
		double xj;
		double zeta;
		//
		//  Update R.
		//
		for ( j = 1; j <= p; j++ )
		{
			xj = x[j-1];
			//
			//  Apply the previous rotations.
			//
			for ( i = 1; i <= j-1; i++ )
			{
				t = c[i-1] * r[i-1+(j-1)*ldz] + s[i-1] * xj;
				xj = c[i-1] * xj - s[i-1] * r[i-1+(j-1)*ldz];
				r[i-1+(j-1)*ldz] = t;
			}
			//
			//  Compute the next rotation.
			//
			//Blas1.drotg ( r, j-1+(j-1)*ldr, &xj, c+j-1, s+j-1 );
			//////////////////////////////////
			double [] drotg = new double [4];
			drotg[0] = r[j-1+(j-1)*ldr];
			drotg[1] = xj;
			drotg[2] = c[j-1];
			drotg[3] = s[j-1];
			///////////////////////////////////
			Blas1.drotg (drotg);
			///////////////////////////////////
			r[j-1+(j-1)*ldr] = drotg[0];
			xj = drotg[1];
			c[j-1] = drotg[2];
			s[j-1] = drotg[3];
			///////////////////////////////////
		}
		//
		//  If required, update Z and RHO.
		//
		for ( j = 1; j <= nz; j++ )
		{
			zeta = y[j-1];
			for ( i = 1; i <= p; i++ )
			{
				t =    c[i-1] * z[i-1+(j-1)*ldz] + s[i-1] * zeta;
				zeta = c[i-1] * zeta   - s[i-1] * z[i-1+(j-1)*ldz];
				z[i-1+(j-1)*ldz] = t;
			}

			azeta = Blas1.r8_abs ( zeta );

			if ( azeta != 0.0 && 0.0 <= rho[j-1] )
			{
				scale = azeta + rho[j-1];
				rho[j-1] = scale * Math.sqrt ( 
						Math.pow ( azeta / scale, 2 ) + Math.pow ( rho[j-1] / scale, 2 ) );
			}

		}

		return;
	}
	//****************************************************************************80

	double dgbco ( double abd[], int lda, int n, int ml, int mu, int ipvt[], 
			double z[] )
	{
		double anorm;
		double ek;
		int i;
		int info;
		int is;
		int j;
		int ju;
		int k;
		int l;
		int la;
		int lm;
		int lz;
		int m;
		int mm;
		double rcond;
		double s;
		double sm;
		double t;
		double wk;
		double wkm;
		double ynorm;
		//
		//  Compute the 1-norm of A.
		//
		anorm = 0.0;
		l = ml + 1;
		is = l + mu;

		for ( j = 1; j <= n; j++ )
		{
			anorm = Blas1.r8_max  ( anorm, Blas1.dasum ( l, abd, is-1+(j-1)*lda, 1 ) );
			if ( ml + 1 < is )
			{
				is = is - 1;
			}
			if ( j <= mu )
			{
				l = l + 1;
			}
			if ( n - ml <= j )
			{
				l = l - 1;
			}
		}
		//
		//  Factor.
		//
		info = dgbfa ( abd, lda, n, ml, mu, ipvt );
		//
		//  RCOND = 1/(norm(A)*(estimate of norm(inverse(A)))).
		//
		//  Estimate = norm(Z)/norm(Y) where  a*z = y  and A'*Y = E.
		//
		//  A' is the transpose of A.  The components of E are
		//  chosen to cause maximum local growth in the elements of W where
		//  U'*W = E.  The vectors are frequently rescaled to avoid
		//  overflow.
		//
		//  Solve U' * W = E.
		//
		ek = 1.0;
		for ( i = 1; i <= n; i++ )
		{
			z[i-1] = 0.0;
		}
		m = ml + mu + 1;
		ju = 0;

		for ( k = 1; k <= n; k++ )
		{
			if ( z[k-1] != 0.0 )
			{
				ek = ek * Blas1.r8_sign ( -z[k-1] );
			}

			if ( Blas1.r8_abs ( abd[m-1+(k-1)*lda] ) < Blas1.r8_abs ( ek - z[k-1] ) )
			{
				s = Blas1.r8_abs ( abd[m-1+(k-1)*lda] ) / Blas1.r8_abs ( ek - z[k-1] );
				for ( i = 1; i <= n; i++ )
				{
					z[i-1] = s * z[i-1];
				}
				ek = s * ek;
			}

			wk = ek - z[k-1];
			wkm = -ek - z[k-1];
			s = Blas1.r8_abs ( wk );
			sm = Blas1.r8_abs ( wkm );

			if ( abd[m-1+(k-1)*lda] != 0.0 )
			{
				wk = wk / abd[m-1+(k-1)*lda];
				wkm = wkm / abd[m-1+(k-1)*lda];
			}
			else
			{
				wk = 1.0;
				wkm = 1.0;
			}

			ju = Blas1.i4_min ( Blas1.i4_max ( ju, mu+ipvt[k-1] ), n );
			mm = m;

			if ( k+1 <= ju )
			{
				for ( j = k+1; j <= ju; j++ )
				{
					mm = mm - 1;
					sm = sm + Blas1.r8_abs ( z[j-1] + wkm * abd[mm-1+(j-1)*lda] );
					z[j-1] = z[j-1] + wk * abd[mm-1+(j-1)*lda];
					s = s + Blas1.r8_abs ( z[j-1] );
				}

				if ( s < sm )
				{
					t = wkm - wk;
					wk = wkm;
					mm = m;
					for ( j = k+1; j <= ju; ju++ )
					{
						mm = mm - 1;
						z[j-1] = z[j-1] + t * abd[mm-1+(j-1)*lda];
					}
				}
			}
			z[k-1] = wk;
		}

		s = Blas1.dasum ( n, z, 0, 1 );

		for ( i = 1; i <= n; i++ )
		{
			z[i-1] = z[i-1] / s;
		}
		//
		//  Solve L' * Y = W.
		//
		for ( k = n; 1 <= k; k-- )
		{
			lm = Blas1.i4_min ( ml, n-k );

			if ( k < m )
			{
				z[k-1] = z[k-1] + Blas1.ddot ( lm, abd, m+(k-1)*lda, 1, z, k, 1 );
			}

			if ( 1.0 < Blas1.r8_abs ( z[k-1] ) )
			{
				s = 1.0 / Blas1.r8_abs ( z[k-1] );
				for ( i = 1; i <= n; i++ )
				{
					z[i-1] = s * z[i-1];
				}
			}
			l = ipvt[k-1];
			t = z[l-1];
			z[l-1] = z[k-1];
			z[k-1] = t;
		}

		s = Blas1.dasum ( n, z, 0, 1 );
		for ( i = 1; i <= n; i++ )
		{
			z[i-1] = z[i-1] / s;
		}
		ynorm = 1.0;
		//
		//  Solve L * V = Y.
		//
		for ( k = 1; k <= n; k++ )
		{
			l = ipvt[k-1];
			t = z[l-1];
			z[l-1] = z[k-1];
			z[k-1] = t;
			lm = Blas1.i4_min ( ml, n-k );

			if ( k < n )
			{
				Blas1.daxpy ( lm, t, abd, m+(k-1)*lda, 1, z, k, 1 );
			}

			if ( 1.0 < Blas1.r8_abs ( z[k-1] ) )
			{
				s = 1.0 / Blas1.r8_abs ( z[k-1] );
				for ( i = 1; i <= n; i++ )
				{
					z[i-1] = s * z[i-1];
				}
				ynorm = s * ynorm;
			}
		}
		s = 1.0 / Blas1.dasum ( n, z, 0, 1 );
		for ( i = 1; i <= n; i++ )
		{
			z[i-1] = s * z[i-1];
		}
		ynorm = s * ynorm;
		//
		//  Solve U * Z = W.
		//
		for ( k = n; 1 <= k; k-- )
		{
			if ( Blas1.r8_abs ( abd[m-1+(k-1)*lda] ) < Blas1.r8_abs ( z[k-1] ) )
			{
				s = Blas1.r8_abs ( abd[m-1+(k-1)*lda] ) / Blas1.r8_abs ( z[k-1] );
				for ( i = 1; i <= n; i++ )
				{
					z[i-1] = s * z[i-1];
				}
				ynorm = s * ynorm;
			}

			if ( abd[m-1+(k-1)*lda] != 0.0 )
			{
				z[k-1] = z[k-1] / abd[m-1+(k-1)*lda];
			}
			else
			{
				z[k-1] = 1.0;
			}

			lm = Blas1.i4_min ( k, m ) - 1;
			la = m - lm;
			lz = k - lm;
			t = -z[k-1];
			Blas1.daxpy ( lm, t, abd, la-1+(k-1)*lda, 1, z, lz-1, 1 );
		}
		//
		//  Make ZNORM = 1.0.
		//
		s = 1.0 / Blas1.dasum ( n, z, 0, 1 );
		for ( i = 1; i <= n; i++ )
		{
			z[i-1] = s * z[i-1];
		}
		ynorm = s * ynorm;

		if ( anorm != 0.0 )
		{
			rcond = ynorm / anorm;
		}
		else
		{
			rcond = 0.0;
		}

		return rcond;
	}
	//****************************************************************************80

/*	void dgbdi ( double abd[], int lda, int n, int ml, int mu, int ipvt[], 
			double det[2] )
	{
		int i;
		int m;

		m = ml + mu + 1;
		det[0] = 1.0;
		det[1] = 0.0;

		for ( i = 1; i <= n; i++ )
		{
			if ( ipvt[i-1] != i )
			{
				det[0] = -det[0];
			}

			det[0] = det[0] * abd[m-1+(i-1)*lda];

			if ( det[0] == 0.0 )
			{
				return;
			}

			while ( Blas1.r8_abs ( det[0] ) < 1.0 )
			{
				det[0] = det[0] * 10.0;
				det[1] = det[1] - 1.0;
			}

			while ( 10.0 <= Blas1.r8_abs ( det[0] ) )
			{
				det[0] = det[0] / 10.0;
				det[1] = det[1] + 1.0;
			}
		}
		return;
	}*/
	//****************************************************************************80

	int dgbfa ( double abd[], int lda, int n, int ml, int mu, int ipvt[] )
	{
		int i;
		int i0;
		int info;
		int j;
		int j0;
		int j1;
		int ju;
		int jz;
		int k;
		int l;
		int lm;
		int m;
		int mm;
		double t;

		m = ml + mu + 1;
		info = 0;
		//
		//  Zero initial fill-in columns.
		//
		j0 = mu + 2;
		j1 = Blas1.i4_min ( n, m ) - 1;

		for ( jz = j0; jz <= j1; jz++ )
		{
			i0 = m + 1 - jz;
			for ( i = i0; i <= ml; i++ )
			{
				abd[i-1+(jz-1)*lda] = 0.0;
			}
		}

		jz = j1;
		ju = 0;
		//
		//  Gaussian elimination with partial pivoting.
		//
		for ( k = 1; k <= n-1; k++ )
		{
			//
			//  Zero out the next fill-in column.
			//
			jz = jz + 1;
			if ( jz <= n )
			{
				for ( i = 1; i <= ml; i++ )
				{
					abd[i-1+(jz-1)*lda] = 0.0;
				}
			}
			//
			//  Find L = pivot index.
			//
			lm = Blas1.i4_min ( ml, n-k );
			l = Blas1.idamax  ( lm+1, abd, m-1+(k-1)*lda, 1 ) + m - 1;
			ipvt[k-1] = l + k - m;
			//
			//  Zero pivot implies this column already triangularized.
			//
			if ( abd[l-1+(k-1)*lda] == 0.0 )
			{
				info = k;
			}
			//
			//  Interchange if necessary.
			//
			else
			{
				if ( l != m )
				{
					t = abd[l-1+(k-1)*lda];
					abd[l-1+(k-1)*lda] = abd[m-1+(k-1)*lda];
					abd[m-1+(k-1)*lda] = t;
				}
				//
				//  Compute multipliers.
				//
				t = -1.0 / abd[m-1+(k-1)*lda];
				Blas1.dscal  ( lm, t, abd, m+(k-1)*lda, 1 );
				//
				//  Row elimination with column indexing.
				//
				ju = Blas1.i4_min ( Blas1.i4_max ( ju, mu+ipvt[k-1] ), n );
				mm = m;

				for ( j = k+1; j <= ju; j++ )
				{
					l = l - 1;
					mm = mm - 1;
					t = abd[l-1+(j-1)*lda];
					if ( l != mm )
					{
						abd[l-1+(j-1)*lda] = abd[mm-1+(j-1)*lda];
						abd[mm-1+(j-1)*lda] = t;
					}
					Blas1.daxpy ( lm, t, abd, m+(k-1)*lda, 1, abd, mm+(j-1)*lda, 1 );
				}

			}

		}

		ipvt[n-1] = n;

		if ( abd[m-1+(n-1)*lda] == 0.0 )
		{
			info = n;
		}

		return info;
	}
	//****************************************************************************80

	void dgbsl ( double abd[], int lda, int n, int ml, int mu, int ipvt[], 
			double b[], int job )
	{
		int k;
		int l;
		int la;
		int lb;
		int lm;
		int m;
		double t;

		m = mu + ml + 1;
		//
		//  JOB = 0, Solve A * x = b.
		//
		//  First solve L * y = b.
		//
		if ( job == 0 )
		{
			if ( 0 < ml )
			{
				for ( k = 1; k <= n-1; k++ )
				{
					lm = Blas1.i4_min ( ml, n-k );
					l = ipvt[k-1];
					t = b[l-1];
					if ( l != k )
					{
						b[l-1] = b[k-1];
						b[k-1] = t;
					}
					Blas1.daxpy ( lm, t, abd, m+(k-1)*lda, 1, b, k, 1 );
				}
			}
			//
			//  Now solve U * x = y.
			//
			for ( k = n; 1 <= k; k-- )
			{
				b[k-1] = b[k-1] / abd[m-1+(k-1)*lda];
				lm = Blas1.i4_min ( k, m ) - 1;
				la = m - lm;
				lb = k - lm;
				t = -b[k-1];
				Blas1.daxpy ( lm, t, abd, la-1+(k-1)*lda, 1, b, lb-1, 1 );
			}
		}
		//
		//  JOB nonzero, solve A' * x = b.
		//
		//  First solve U' * y = b.
		//
		else
		{
			for ( k = 1; k <= n; k++ )
			{
				lm = Blas1.i4_min ( k, m ) - 1;
				la = m - lm;
				lb = k - lm;
				t = Blas1.ddot ( lm, abd, la-1+(k-1)*lda, 1, b, lb-1, 1 );
				b[k-1] = ( b[k-1] - t ) / abd[m-1+(k-1)*lda];
			}
			//
			//  Now solve L' * x = y.
			//
			if ( 0 < ml )
			{
				for ( k = n-1; 1 <= k; k-- )
				{
					lm = Blas1.i4_min ( ml, n-k );
					b[k-1] = b[k-1] + Blas1.ddot ( lm, abd, m+(k-1)*lda, 1, b, k, 1 );
					l = ipvt[k-1];
					if ( l != k )
					{
						t = b[l-1];
						b[l-1] = b[k-1];
						b[k-1] = t;
					}
				}
			}
		}

		return;
	}
	//****************************************************************************80

	double dgeco ( double a[], int lda, int n, int ipvt[], double z[] )
	{
		double anorm;
		double ek;
		int i;
		int info;
		int j;
		int k;
		int l;
		double rcond;
		double s;
		double sm;
		double t;
		double wk;
		double wkm;
		double ynorm;
		//
		//  Compute the L1 norm of A.
		//
		anorm = 0.0;
		for ( j = 1; j <= n; j++ )
		{
			anorm = Blas1.r8_max  ( anorm, Blas1.dasum ( n, a, 0+(j-1)*lda, 1 ) );
		}
		//
		//  Compute the LU factorization.
		//
		info = dgefa ( a, lda, n, ipvt );
		//
		//  RCOND = 1 / ( norm(A) * (estimate of norm(inverse(A))) )
		//
		//  estimate of norm(inverse(A)) = norm(Z) / norm(Y)
		//
		//  where
		//	    A * Z = Y
		//  and
		//	    A' * Y = E
		//
		//  The components of E are chosen to cause maximum local growth in the
		//  elements of W, where U'*W = E.  The vectors are frequently rescaled
		//  to avoid overflow.
		//
		//  Solve U' * W = E.
		//
		ek = 1.0;
		for ( i = 1; i <= n; i++ )
		{
			z[i-1] = 0.0;
		}

		for ( k = 1; k <= n; k++ )
		{
			if ( z[k-1] != 0.0 )
			{
				ek = ek * Blas1.r8_sign ( -z[k-1] );
			}

			if ( Blas1.r8_abs ( a[k-1+(k-1)*lda] ) < Blas1.r8_abs ( ek - z[k-1] ) )
			{
				s = Blas1.r8_abs ( a[k-1+(k-1)*lda] ) / Blas1.r8_abs ( ek - z[k-1] );
				for ( i = 1; i <= n; i++ )
				{
					z[i-1] = s * z[i-1];
				}
				ek = s * ek;
			}

			wk = ek - z[k-1];
			wkm = -ek - z[k-1];
			s = Blas1.r8_abs ( wk );
			sm = Blas1.r8_abs ( wkm );

			if ( a[k-1+(k-1)*lda] != 0.0 )
			{
				wk = wk / a[k-1+(k-1)*lda];
				wkm = wkm / a[k-1+(k-1)*lda];
			}
			else
			{
				wk = 1.0;
				wkm = 1.0;
			}

			if ( k+1 <= n )
			{
				for ( j = k+1; j <= n; j++ )
				{
					sm = sm + Blas1.r8_abs ( z[j-1] + wkm * a[k-1+(j-1)*lda] );
					z[j-1] = z[j-1] + wk * a[k-1+(j-1)*lda];
					s = s + Blas1.r8_abs ( z[j-1] );
				}

				if ( s < sm )
				{
					t = wkm - wk;
					wk = wkm;
					for ( i = k+1; i <= n; i++ )
					{
						z[i-1] = z[i-1] + t * a[k-1+(i-1)*lda];
					}
				}
			}
			z[k-1] = wk;
		}

		t = Blas1.dasum ( n, z, 0, 1 );
		for ( i = 1; i <= n; i++ )
		{
			z[i-1] = z[i-1] / t;
		}
		//
		//  Solve L' * Y = W
		//
		for ( k = n; 1 <= k; k-- )
		{
			z[k-1] = z[k-1] + Blas1.ddot ( n - k, a, k+(k-1)*lda, 1, z, k, 1 );

			if ( 1.0 < Blas1.r8_abs ( z[k-1] ) )
			{
				t = Blas1.r8_abs ( z[k-1] );
				for ( i = 1; i <= n; i++ )
				{
					z[i-1] = z[i-1] / t;
				}
			}

			l = ipvt[k-1];

			t = z[l-1];
			z[l-1] = z[k-1];
			z[k-1] = t;
		}
		t = Blas1.dasum ( n, z, 0, 1 );
		for ( i = 1; i <= n; i++ )
		{
			z[i-1] = z[i-1] / t;
		}

		ynorm = 1.0;
		//
		//  Solve L * V = Y.
		//
		for ( k = 1; k <= n; k++ )
		{
			l = ipvt[k-1];

			t = z[l-1];
			z[l-1] = z[k-1];
			z[k-1] = t;

			for ( i = k+1; i <= n; i++ )
			{
				z[i-1] = z[i-1] + t * a[i-1+(k-1)*lda];
			}

			if ( 1.0 < Blas1.r8_abs ( z[k-1] ) )
			{
				ynorm = ynorm / Blas1.r8_abs ( z[k-1] );
				t = Blas1.r8_abs ( z[k-1] );
				for ( i = 1; i <= n; i++ )
				{
					z[i-1] = z[i-1] / t;
				}
			}
		}
		s = Blas1.dasum ( n, z, 0, 1 );
		for ( i = 1; i <= n; i++ )
		{
			z[i-1] = z[i-1] / s;
		}
		ynorm = ynorm / s;
		//
		//  Solve U * Z = V.
		//
		for ( k = n; 1 <= k; k-- )
		{
			if ( Blas1.r8_abs ( a[k-1+(k-1)*lda] ) < Blas1.r8_abs ( z[k-1] ) )
			{
				s = Blas1.r8_abs ( a[k-1+(k-1)*lda] ) / Blas1.r8_abs ( z[k-1] );
				for ( i = 1; i <= n; i++ )
				{
					z[i-1] = s * z[i-1];
				}
				ynorm = s * ynorm;
			}

			if ( a[k-1+(k-1)*lda] != 0.0 )
			{
				z[k-1] = z[k-1] / a[k-1+(k-1)*lda];
			}
			else
			{
				z[k-1] = 1.0;
			}
			for ( i = 1; i <= k-1; i++ )
			{
				z[i-1] = z[i-1] - z[k-1] * a[i-1+(k-1)*lda];
			}
		}
		//
		//  Normalize Z in the L1 norm.
		//
		s = 1.0 / Blas1.dasum ( n, z, 0, 1 );
		for ( i = 1; i <= n; i++ )
		{
			z[i-1] = s * z[i-1];
		}
		ynorm = s * ynorm;

		if ( anorm != 0.0 )
		{
			rcond = ynorm / anorm;
		}
		else
		{
			rcond = 0.0;
		}

		return rcond;
	}
	//****************************************************************************80

	void dgedi ( double a[], int lda, int n, int ipvt[], double det[], 
			double work[], int job )
	{
		int i;
		int j;
		int k;
		int l;
		double t;
		//
		//  Compute the determinant.
		//
		if ( job / 10 != 0 )
		{
			det[0] = 1.0;
			det[1] = 0.0;

			for ( i = 1; i <= n; i++ )
			{
				if ( ipvt[i-1] != i )
				{
					det[0] = -det[0];
				}
				det[0] = det[0] * a[i-1+(i-1)*lda];

				if ( det[0] == 0.0 )
				{
					break;
				}

				while ( Blas1.r8_abs ( det[0] ) < 1.0 )
				{
					det[0] = det[0] * 10.0;
					det[1] = det[1] - 1.0;
				}
				while ( 10.0 <= Blas1.r8_abs ( det[0] ) )
				{
					det[0] = det[0] / 10.0;
					det[1] = det[1] + 1.0;
				}
			}
		}
		//
		//  Compute inverse(U).
		//
		if ( ( job % 10 ) != 0 )
		{
			for ( k = 1; k <= n; k++ )
			{
				a[k-1+(k-1)*lda] = 1.0 / a[k-1+(k-1)*lda];
				t = -a[k-1+(k-1)*lda];
				Blas1.dscal  ( k-1, t, a, 0+(k-1)*lda, 1 );

				for ( j = k+1; j <= n; j++ )
				{
					t = a[k-1+(j-1)*lda];
					a[k-1+(j-1)*lda] = 0.0;
					Blas1.daxpy ( k, t, a, 0+(k-1)*lda, 1, a, 0+(j-1)*lda, 1 );
				}
			}
			//
			//  Form inverse(U) * inverse(L).
			//
			for ( k = n-1; 1 <= k; k-- )
			{
				for ( i = k+1; i <= n; i++ )
				{
					work[i-1] = a[i-1+(k-1)*lda];
					a[i-1+(k-1)*lda] = 0.0;
				}

				for ( j = k+1; j <= n; j++ )
				{
					t = work[j-1];
					Blas1.daxpy ( n, t, a, 0+(j-1)*lda, 1, a, 0+(k-1)*lda, 1 );
				}

				l = ipvt[k-1];
				if ( l != k )
				{
					Blas1.dswap  ( n, a, 0+(k-1)*lda, 1, a, 0+(l-1)*lda, 1 );
				}
			}
		}

		return;
	}
	//****************************************************************************80

	int dgefa ( double a[], int lda, int n, int ipvt[] )
	{
		int info;
		int j;
		int k;
		int l;
		double t;
		//
		//  Gaussian elimination with partial pivoting.
		//
		info = 0;

		for ( k = 1; k <= n-1; k++ )
		{
			//
			//  Find L = pivot index.
			//
			l = Blas1.idamax  ( n-k+1, a, (k-1)+(k-1)*lda, 1 ) + k - 1;
			ipvt[k-1] = l;
			//
			//  Zero pivot implies this column already triangularized.
			//
			if ( a[l-1+(k-1)*lda] == 0.0 )
			{
				info = k;
				continue;
			}
			//
			//  Interchange if necessary.
			//
			if ( l != k )
			{
				t = a[l-1+(k-1)*lda];
				a[l-1+(k-1)*lda] = a[k-1+(k-1)*lda];
				a[k-1+(k-1)*lda] = t;
			}
			//
			//  Compute multipliers.
			//
			t = -1.0 / a[k-1+(k-1)*lda];

			Blas1.dscal  ( n-k, t, a, k+(k-1)*lda, 1 );
			//
			//  Row elimination with column indexing.
			//
			for ( j = k+1; j <= n; j++ )
			{
				t = a[l-1+(j-1)*lda];
				if ( l != k )
				{
					a[l-1+(j-1)*lda] = a[k-1+(j-1)*lda];
					a[k-1+(j-1)*lda] = t;
				}
				Blas1.daxpy ( n-k, t, a, k+(k-1)*lda, 1, a, k+(j-1)*lda, 1 );
			}

		}

		ipvt[n-1] = n;

		if ( a[n-1+(n-1)*lda] == 0.0 )
		{
			info = n;
		}

		return info;
	}
	//****************************************************************************80

	void dgesl ( double a[], int lda, int n, int ipvt[], double b[], int job )
	{
		int k;
		int l;
		double t;
		//
		//  Solve A * X = B.
		//
		if ( job == 0 )
		{
			for ( k = 1; k <= n-1; k++ )
			{
				l = ipvt[k-1];
				t = b[l-1];

				if ( l != k )
				{
					b[l-1] = b[k-1];
					b[k-1] = t;
				}

				Blas1.daxpy ( n-k, t, a, k+(k-1)*lda, 1, b, k, 1 );

			}

			for ( k = n; 1 <= k; k-- )
			{
				b[k-1] = b[k-1] / a[k-1+(k-1)*lda];
				t = -b[k-1];
				Blas1.daxpy ( k-1, t, a, 0+(k-1)*lda, 1, b, 0, 1 );
			}
		}
		//
		//  Solve A' * X = B.
		//
		else
		{
			for ( k = 1; k <= n; k++ )
			{
				t = Blas1.ddot ( k-1, a, 0+(k-1)*lda, 1, b, 0, 1 );
				b[k-1] = ( b[k-1] - t ) / a[k-1+(k-1)*lda];
			}

			for ( k = n-1; 1 <= k; k-- )
			{
				b[k-1] = b[k-1] + Blas1.ddot ( n-k, a, k+(k-1)*lda, 1, b, k, 1 );
				l = ipvt[k-1];

				if ( l != k )
				{
					t = b[l-1];
					b[l-1] = b[k-1];
					b[k-1] = t;
				}
			}
		}
		return;
	}
	//****************************************************************************80

	int dgtsl ( int n, double c[], double d[], double e[], double b[] )
	{
		int info;
		int k;
		double t;

		info = 0;
		c[0] = d[0];

		if ( 2 <= n )
		{
			d[0] = e[0];
			e[0] = 0.0;
			e[n-1] = 0.0;

			for ( k = 1; k <= n - 1; k++ )
			{
				//
				//  Find the larger of the two rows.
				//
				if ( Blas1.r8_abs ( c[k-1] ) <= Blas1.r8_abs ( c[k] ) )
				{
					//
					//  Interchange rows.
					//
					t = c[k];
					c[k] = c[k-1];
					c[k-1] = t;

					t = d[k];
					d[k] = d[k-1];
					d[k-1] = t;

					t = e[k];
					e[k] = e[k-1];
					e[k-1] = t;

					t = b[k];
					b[k] = b[k-1];
					b[k-1] = t;
				}
				//
				//  Zero elements.
				//
				if ( c[k-1] == 0.0 )
				{
					info = k;
					return info;
				}

				t = -c[k] / c[k-1];
				c[k] = d[k] + t * d[k-1];
				d[k] = e[k] + t * e[k-1];
				e[k] = 0.0;
				b[k] = b[k] + t * b[k-1];
			}
		}

		if ( c[n-1] == 0.0 )
		{
			info = n;
			return info;
		}
		//
		//  Back solve.
		//
		b[n-1] = b[n-1] / c[n-1];

		if ( 1 < n )
		{
			b[n-2] = ( b[n-2] - d[n-2] * b[n-1] ) / c[n-2];

			for ( k = n-2; 1 <= k; k-- )
			{
				b[k-1] = ( b[k-1] - d[k-1] * b[k] - e[k-1] * b[k+1] ) / c[k-1];
			}

		}

		return info;
	}
	//****************************************************************************80

	double dpbco ( double abd[], int lda, int n, int m, double z[] )
	{
		double anorm;
		double ek;
		int i;
		int info;
		int j;
		int j2;
		int k;
		int l;
		int la;
		int lb;
		int lm;
		int mu;
		double rcond;
		double s;
		double sm;
		double t;
		double wk;
		double wkm;
		double ynorm;
		//
		//  Find the norm of A.
		//
		for ( j = 1; j <= n; j++ )
		{
			l = Blas1.i4_min ( j, m+1 );
			mu = Blas1.i4_max ( m+2-j, 1 );
			z[j-1] = Blas1.dasum ( l, abd, mu-1+(j-1)*lda, 1 );
			k = j - l;
			for ( i = mu; i <= m; i++ )
			{
				k = k + 1;
				z[k-1] = z[k-1] + Blas1.r8_abs ( abd[i-1+(j-1)*lda] );
			}
		}
		anorm = 0.0;
		for ( i = 1; i <= n; i++ )
		{
			anorm = Blas1.r8_max  ( anorm, z[i-1] );
		}
		//
		//  Factor.
		//
		info = dpbfa ( abd, lda, n, m );

		if ( info != 0 )
		{
			rcond = 0.0;
			return rcond;
		}
		//
		//  RCOND = 1/(norm(A)*(estimate of norm(inverse(A)))).
		//
		//  Estimate = norm(Z)/norm(Y) where A*Z = Y and A*Y = E.
		//
		//  The components of E are chosen to cause maximum local
		//  growth in the elements of W where R'*W = E.
		//
		//  The vectors are frequently rescaled to avoid overflow.
		//
		//  Solve R' * W = E.
		//
		ek = 1.0;
		for ( i = 1; i <= n; i++ )
		{
			z[i-1] = 0.0;
		}

		for ( k = 1; k <= n; k++ )
		{
			if ( z[k-1] != 0.0 )
			{
				ek = ek * Blas1.r8_sign ( -z[k-1] );
			}

			if ( abd[m+(k-1)*lda] < Blas1.r8_abs ( ek - z[k-1] ) )
			{
				s = abd[m+(k-1)*lda] / Blas1.r8_abs ( ek - z[k-1] );
				for ( i = 1; i <= n; i++ )
				{
					z[i-1] = s * z[i-1];
				}
				ek = s * ek;
			}
			wk = ek - z[k-1];
			wkm = -ek - z[k-1];
			s = Blas1.r8_abs ( wk );
			sm = Blas1.r8_abs ( wkm );
			wk = wk / abd[m+(k-1)*lda];
			wkm = wkm / abd[m+(k-1)*lda];
			j2 = Blas1.i4_min ( k+m, n );
			i = m + 1;

			if ( k+1 <= j2 )
			{
				for ( j = k+1; j <= j2; j++ )
				{
					i = i - 1;
					sm = sm + Blas1.r8_abs ( z[j-1] + wkm * abd[i-1+(j-1)*lda] );
					z[j-1] = z[j-1] + wk * abd[i-1+(j-1)*lda];
					s = s + Blas1.r8_abs ( z[j-1] );
				}

				if ( s < sm )
				{
					t = wkm - wk;
					wk = wkm;
					i = m + 1;

					for ( j = k+1; j <= j2; j++ )
					{
						i = i - 1;
						z[j-1] = z[j-1] + t * abd[i-1+(j-1)*lda];
					}
				}
			}
			z[k-1] = wk;
		}
		s = Blas1.dasum ( n, z, 0, 1 );
		for ( i = 1; i <= n; i++ )
		{
			z[i-1] = z[i-1] / s;
		}
		//
		//  Solve R * Y = W.
		//
		for ( k = n; 1 <= k; k-- )
		{
			if ( abd[m+(k-1)*lda] < Blas1.r8_abs ( z[k-1] ) )
			{
				s = abd[m+(k-1)*lda] / Blas1.r8_abs ( z[k-1] );
				for ( i = 1; i <= n; i++ )
				{
					z[i-1] = s * z[i-1];
				}
			}
			z[k-1] = z[k-1] / abd[m+(k-1)*lda];
			lm = Blas1.i4_min ( k-1, m );
			la = m + 1 - lm;
			lb = k - lm;
			t = -z[k-1];
			Blas1.daxpy ( lm, t, abd, la-1+(k-1)*lda, 1, z, lb-1, 1 );
		}
		s = Blas1.dasum ( n, z, 0, 1 );
		for ( i = 1; i <= n; i++ )
		{
			z[i-1] = z[i-1] / s;
		}
		ynorm = 1.0;
		//
		//  Solve R' * V = Y.
		//
		for ( k = 1; k <= n; k++ )
		{
			lm = Blas1.i4_min ( k-1, m );
			la = m + 1 - lm;
			lb = k - lm;

			z[k-1] = z[k-1] - Blas1.ddot ( lm, abd, la-1+(k-1)*lda, 1, z, lb-1, 1 );

			if ( abd[m+(k-1)*lda] < Blas1.r8_abs ( z[k-1] ) )
			{
				s = abd[m+(k-1)*lda] / Blas1.r8_abs ( z[k-1] );
				for ( i = 1; i <= n; i++ )
				{
					z[i-1] = s * z[i-1];
				}
				ynorm = s * ynorm;
			}
			z[k-1] = z[k-1] / abd[m+(k-1)*lda];
		}
		s = 1.0 / Blas1.dasum ( n, z, 0, 1 );
		for ( i = 1; i <= n; i++ )
		{
			z[i-1] = s * z[i-1];
		}
		ynorm = s * ynorm;
		//
		//  Solve R * Z = W.
		//
		for ( k = n; 1 <= k; k-- )
		{
			if ( abd[m+(k-1)*lda] < Blas1.r8_abs ( z[k-1] ) )
			{
				s = abd[m+(k-1)*lda] / Blas1.r8_abs ( z[k-1] );
				for ( i = 1; i <= n; i++ )
				{
					z[i-1] = s * z[i-1];
				}
				ynorm = s * ynorm;
			}
			z[k-1] = z[k-1] / abd[m+(k-1)*lda];
			lm = Blas1.i4_min ( k-1, m );
			la = m + 1 - lm;
			lb = k - lm;
			t = -z[k-1];
			Blas1.daxpy ( lm, t, abd, la-1+(k-1)*lda, 1, z, lb-1, 1 );
		}
		//
		//  Make ZNORM = 1.0.
		//
		s = 1.0 / Blas1.dasum ( n, z, 0, 1 );
		for ( i = 1; i <= n; i++ )
		{
			z[i-1] = s * z[i-1];
		}
		ynorm = s * ynorm;

		if ( anorm != 0.0 )
		{
			rcond = ynorm / anorm;
		}
		else
		{
			rcond = 0.0;
		}

		return rcond;
	}
	//****************************************************************************80

	void dpbdi ( double abd[], int lda, int n, int m, double det[] )
	{
		int i;
		double s;
		//
		//  Compute the determinant.
		//
		det[0] = 1.0;
		det[1] = 0.0;
		s = 10.0;

		for ( i = 1; i <= n; i++ )
		{
			det[0] = det[0] * abd[m+(i-1)*lda] * abd[m+(i-1)*lda];

			if ( det[0] == 0.0 )
			{
				return;
			}

			while ( det[0] < 1.0 )
			{
				det[0] = det[0] * s;
				det[1] = det[1] - 1.0;
			}

			while ( s <= det[0] )
			{
				det[0] = det[0] / s;
				det[1] = det[1] + 1.0;
			}

		}

		return;
	}
	//****************************************************************************80

	int dpbfa ( double abd[], int lda, int n, int m )
	{
		int ik;
		int info;
		int j;
		int jk;
		int k;
		int mu;
		double s;
		double t;

		for ( j = 1; j <= n; j++ )
		{
			s = 0.0;
			ik = m + 1;
			jk = Blas1.i4_max ( j - m, 1 );
			mu = Blas1.i4_max ( m + 2 - j, 1 );

			for ( k = mu; k <= m; k++ )
			{
				t = abd[k-1+(j-1)*lda] 
						- Blas1.ddot ( k-mu, abd, ik-1+(jk-1)*lda, 1, abd, mu-1+(j-1)*lda, 1 );
				t = t / abd[m+(jk-1)*lda];
				abd[k-1+(j-1)*lda] = t;
				s = s + t * t;
				ik = ik - 1;
				jk = jk + 1;
			}
			s = abd[m+(j-1)*lda] - s;

			if ( s <= 0.0 )
			{
				info = j;
				return info;
			}
			abd[m+(j-1)*lda] = Math.sqrt ( s );
		}
		info = 0;

		return info;
	}
	//****************************************************************************80

	void dpbsl ( double abd[], int lda, int n, int m, double b[] )
	{
		int k;
		int la;
		int lb;
		int lm;
		double t;
		//
		//  Solve R'*Y = B.
		//
		for ( k = 1; k <= n; k++ )
		{
			lm = Blas1.i4_min ( k-1, m );
			la = m + 1 - lm;
			lb = k - lm;
			t = Blas1.ddot ( lm, abd, la-1+(k-1)*lda, 1, b, lb-1, 1 );
			b[k-1] = ( b[k-1] - t ) / abd[m+(k-1)*lda];
		}
		//
		//  Solve R*X = Y.
		//
		for ( k = n; 1 <= k; k-- )
		{
			lm = Blas1.i4_min ( k-1, m );
			la = m + 1 - lm;
			lb = k - lm;
			b[k-1] = b[k-1] / abd[m+(k-1)*lda];
			t = -b[k-1];
			Blas1.daxpy ( lm, t, abd, la-1+(k-1)*lda, 1, b, lb-1, 1 );
		}

		return;
	}
	//****************************************************************************80

	double dpoco ( double a[], int lda, int n, double z[] )
	{
		double anorm;
		double ek;
		int i;
		int info;
		int j;
		int k;
		double rcond;
		double s;
		double sm;
		double t;
		double wk;
		double wkm;
		double ynorm;
		//
		//  Find norm of A using only upper half.
		//
		for ( j = 1; j <= n; j++ )
		{
			z[j-1] = Blas1.dasum ( j, a, 0+(j-1)*lda, 1 );
			for ( i = 1; i <= j-1; i++ )
			{
				z[i-1] = z[i-1] + Blas1.r8_abs ( a[i-1+(j-1)*lda] );
			}
		}

		anorm = 0.0;
		for ( i = 1; i <= n; i++ )
		{
			anorm = Blas1.r8_max  ( anorm, z[i-1] );
		}
		//
		//  Factor.
		//
		info = dpofa ( a, lda, n );

		if ( info != 0 )
		{
			rcond = 0.0;
			return rcond;
		}
		//
		//  RCOND = 1/(norm(A)*(estimate of norm(inverse(A)))).
		//
		//  Estimate = norm(Z)/norm(Y) where A*Z = Y and A*Y = E.
		//
		//  The components of E are chosen to cause maximum local
		//  growth in the elements of W where R'*W = E.
		//
		//  The vectors are frequently rescaled to avoid overflow.
		//
		//  Solve R' * W = E.
		//
		ek = 1.0;
		for ( i = 1; i <= n; i++ )
		{
			z[i-1] = 0.0;
		}

		for ( k = 1; k <= n; k++ )
		{
			if ( z[k-1] != 0.0 )
			{
				ek = ek * Blas1.r8_sign ( -z[k-1] );
			}

			if ( a[k-1+(k-1)*lda] < Blas1.r8_abs ( ek - z[k-1] ) )
			{
				s = a[k-1+(k-1)*lda] / Blas1.r8_abs ( ek - z[k-1] );
				for ( i = 1; i <= n; i++ )
				{
					z[i-1] = s * z[i-1];
				}
				ek = s * ek;
			}

			wk = ek - z[k-1];
			wkm = -ek - z[k-1];
			s = Blas1.r8_abs ( wk );
			sm = Blas1.r8_abs ( wkm );
			wk = wk / a[k-1+(k-1)*lda];
			wkm = wkm / a[k-1+(k-1)*lda];

			if ( k + 1 <= n )
			{
				for ( j = k+1; j <= n; j++ )
				{
					sm = sm + Blas1.r8_abs ( z[j-1] + wkm * a[k-1+(j-1)*lda] );
					z[j-1] = z[j-1] + wk * a[k-1+(j-1)*lda];
					s = s + Blas1.r8_abs ( z[j-1] );
				}

				if ( s < sm )
				{
					t = wkm - wk;
					wk = wkm;
					for ( j = k+1; j <= n; j++ )
					{
						z[j-1] = z[j-1] + t * a[k-1+(j-1)*lda];
					}
				}
			}
			z[k-1] = wk;
		}
		s = Blas1.dasum ( n, z, 0, 1 );
		for ( i = 1; i <= n; i++ )
		{
			z[i-1] = z[i-1] / s;
		}
		//
		//  Solve R * Y = W.
		//
		for ( k = n; 1 <= k; k-- )
		{
			if ( a[k-1+(k-1)*lda] < Blas1.r8_abs ( z[k-1] ) )
			{
				s = a[k-1+(k-1)*lda] / Blas1.r8_abs ( z[k-1] );
				for ( i = 1; i <= n; i++ )
				{
					z[i-1] = s * z[i-1];
				}
			}
			z[k-1] = z[k-1] / a[k-1+(k-1)*lda];
			t = -z[k-1];
			Blas1.daxpy ( k-1, t, a, 0+(k-1)*lda, 1, z, 0, 1 );
		}
		s = Blas1.dasum ( n, z, 0, 1 );
		for ( i = 1; i <= n; i++ )
		{
			z[i-1] = z[i-1] / s;
		}
		ynorm = 1.0;
		//
		//  Solve R' * V = Y.
		//
		for ( k = 1; k <= n; k++ )
		{
			z[k-1] = z[k-1] - Blas1.ddot ( k-1, a, 0+(k-1)*lda, 1, z, 0, 1 );

			if ( a[k-1+(k-1)*lda] < Blas1.r8_abs ( z[k-1] ) )
			{
				s = a[k-1+(k-1)*lda] / Blas1.r8_abs ( z[k-1] );
				for ( i = 1; i <= n; i++ )
				{
					z[i-1] = s * z[i-1];
				}
				ynorm = s * ynorm;
			}
			z[k-1] = z[k-1] / a[k-1+(k-1)*lda];
		}
		s = 1.0 / Blas1.dasum ( n, z, 0, 1 );
		for ( i = 1; i <= n; i++ )
		{
			z[i-1] = s * z[i-1];
		}
		ynorm = s * ynorm;
		//
		//  Solve R * Z = V.
		//
		for ( k = n; 1 <= k; k-- )
		{
			if ( a[k-1+(k-1)*lda] < Blas1.r8_abs ( z[k-1] ) )
			{
				s = a[k-1+(k-1)*lda] / Blas1.r8_abs ( z[k-1] );
				for ( i = 1; i <= n; i++ )
				{
					z[i-1] = s * z[i-1];
				}
				ynorm = s * ynorm;
			}
			z[k-1] = z[k-1] / a[k-1+(k-1)*lda];
			t = -z[k-1];
			Blas1.daxpy ( k-1, t, a, 0+(k-1)*lda, 1, z, 0, 1 );
		}
		//
		//  Make ZNORM = 1.0.
		//
		s = 1.0 / Blas1.dasum ( n, z, 0, 1 );
		for ( i = 1; i <= n; i++ )
		{
			z[i-1] = s * z[i-1];
		}
		ynorm = s * ynorm;

		if ( anorm != 0.0 )
		{
			rcond = ynorm / anorm;
		}
		else
		{
			rcond = 0.0;
		}

		return rcond;
	}
	//****************************************************************************80

	static void dpodi ( double a[], int lda, int n, double det[], int job )
	{
		int i;
		int j;
		int k;
		double s;
		double t;
		//
		//  Compute the determinant.
		//
		if ( job / 10 != 0 )
		{
			det[0] = 1.0;
			det[1] = 0.0;
			s = 10.0;

			for ( i = 1; i <= n; i++ )
			{
				det[0] = det[0] * a[i-1+(i-1)*lda] * a[i-1+(i-1)*lda];

				if ( det[0] == 0.0 )
				{
					break;
				}

				while ( det[0] < 1.0 )
				{
					det[0] = det[0] * s;
					det[1] = det[1] - 1.0;
				}

				while ( s <= det[0] )
				{
					det[0] = det[0] / s;
					det[1] = det[1] + 1.0;
				}
			}
		}
		//
		//  Compute inverse(R).
		//
		if ( ( job % 10 ) != 0 )
		{
			for ( k = 1; k <= n; k++ )
			{
				a[k-1+(k-1)*lda] = 1.0 / a[k-1+(k-1)*lda];
				t = -a[k-1+(k-1)*lda];
				Blas1.dscal  ( k-1, t, a, 0+(k-1)*lda, 1 );

				for ( j = k+1; j <= n; j++ )
				{
					t = a[k-1+(j-1)*lda];
					a[k-1+(j-1)*lda] = 0.0;
					Blas1.daxpy ( k, t, a, 0+(k-1)*lda, 1, a, 0+(j-1)*lda, 1 );
				}
			}
			//
			//  Form inverse(R) * (inverse(R))'.
			//
			for ( j = 1; j <= n; j++ )
			{
				for ( k = 1; k <= j-1; k++ )
				{
					t = a[k-1+(j-1)*lda];
					Blas1.daxpy ( k, t, a, 0+(j-1)*lda, 1, a, 0+(k-1)*lda, 1 );
				}
				t = a[j-1+(j-1)*lda];
				Blas1.dscal  ( j, t, a, 0+(j-1)*lda, 1 );
			}
		}

		return;
	}
	//****************************************************************************80

	static int dpofa ( double a[], int lda, int n )
	{
		int info;
		int j;
		int k;
		double s;
		double t;

		for ( j = 1; j <= n; j++ )
		{
			s = 0.0;

			for ( k = 1; k <= j-1; k++ )
			{
				t = a[k-1+(j-1)*lda] - Blas1.ddot ( k-1, a, 0+(k-1)*lda, 1, a, 0+(j-1)*lda, 1 );
				t = t / a[k-1+(k-1)*lda];
				a[k-1+(j-1)*lda] = t;
				s = s + t * t;
			}

			s = a[j-1+(j-1)*lda] - s;

			if ( s <= 0.0 )
			{
				info = j;
				return info;
			}

			a[j-1+(j-1)*lda] = Math.sqrt ( s );
		}

		info = 0;

		return info;
	}
	//****************************************************************************80

	void dposl ( double a[], int lda, int n, double b[] )
	{
		int k;
		double t;
		//
		//  Solve R' * Y = B.
		//
		for ( k = 1; k <= n; k++ )
		{
			t = Blas1.ddot ( k-1, a, 0+(k-1)*lda, 1, b, 0, 1 );
			b[k-1] = ( b[k-1] - t ) / a[k-1+(k-1)*lda];
		}
		//
		//  Solve R * X = Y.
		//
		for ( k = n; 1 <= k; k-- )
		{
			b[k-1] = b[k-1] / a[k-1+(k-1)*lda];
			t = -b[k-1];
			Blas1.daxpy ( k-1, t, a, 0+(k-1)*lda, 1, b, 0, 1 );
		}

		return;
	}
	//****************************************************************************80

	double dppco ( double ap[], int n, double z[] )
	{
		double anorm;
		double ek;
		int i;
		int ij;
		int info;
		int j;
		int j1;
		int k;
		int kj;
		int kk;
		double rcond;
		double s;
		double sm;
		double t;
		double wk;
		double wkm;
		double ynorm;
		//
		//  Find the norm of A.
		//
		j1 = 1;
		for ( j = 1; j <= n; j++ )
		{
			z[j-1] = Blas1.dasum ( j, ap, j1-1, 1 );
			ij = j1;
			j1 = j1 + j;
			for ( i = 1; i <= j-1; i++ )
			{
				z[i-1] = z[i-1] + Blas1.r8_abs ( ap[ij-1] );
				ij = ij + 1;
			}
		}
		anorm = 0.0;
		for ( i = 1; i <= n; i++ )
		{
			anorm = Blas1.r8_max  ( anorm, z[i-1] );
		}
		//
		//  Factor.
		//
		info = dppfa ( ap, n );

		if ( info != 0 )
		{
			rcond = 0.0;
			return rcond;
		}
		//
		//  RCOND = 1/(norm(A)*(estimate of norm(inverse(A)))).
		//
		//  Estimate = norm(Z)/norm(Y) where A * Z = Y and A * Y = E.
		//
		//  The components of E are chosen to cause maximum local
		//  growth in the elements of W where R'*W = E.
		//
		//  The vectors are frequently rescaled to avoid overflow.
		//
		//  Solve R' * W = E.
		//
		ek = 1.0;
		for ( i = 1; i <= n; i++ )
		{
			z[i-1] = 0.0;
		}

		kk = 0;

		for ( k = 1; k <= n; k++ )
		{
			kk = kk + k;

			if ( z[k-1] != 0.0 )
			{
				ek = ek * Blas1.r8_sign ( -z[k-1] );
			}

			if ( ap[kk-1] < Blas1.r8_abs ( ek - z[k-1] ) )
			{
				s = ap[kk-1] / Blas1.r8_abs ( ek - z[k-1] );
				for ( i = 1; i <= n; i++ )
				{
					z[i-1] = s * z[i-1];
				}
				ek = s * ek;
			}
			wk = ek - z[k-1];
			wkm = -ek - z[k-1];
			s = Blas1.r8_abs ( wk );
			sm = Blas1.r8_abs ( wkm );
			wk = wk / ap[kk-1];
			wkm = wkm / ap[kk-1];
			kj = kk + k;

			if ( k + 1 <= n )
			{
				for ( j = k + 1; j <= n; j++ )
				{
					sm = sm + Blas1.r8_abs ( z[j-1] + wkm * ap[kj-1] );
					z[j-1] = z[j-1] + wk * ap[kj-1];
					s = s + Blas1.r8_abs ( z[j-1] );
					kj = kj + j;
				}

				if ( s < sm )
				{
					t = wkm - wk;
					wk = wkm;
					kj = kk + k;

					for ( j = k+1; j <= n; j++ )
					{
						z[j-1] = z[j-1] + t * ap[kj-1];
						kj = kj + j;
					}
				}
			}
			z[k-1] = wk;
		}

		s = Blas1.dasum ( n, z, 0, 1 );
		for ( i = 1; i <= n; i++ )
		{
			z[i-1] = z[i-1] / s;
		}
		//
		//  Solve R * Y = W.
		//
		for ( k = n; 1 <= k; k-- )
		{
			if ( ap[kk-1] < Blas1.r8_abs ( z[k-1] ) )
			{
				s = ap[kk-1] / Blas1.r8_abs ( z[k-1] );
				for ( i = 1; i <= n; i++ )
				{
					z[i-1] = s * z[i-1];
				}
			}
			z[k-1] = z[k-1] / ap[kk-1];
			kk = kk - k;
			t = -z[k-1];
			Blas1.daxpy ( k-1, t, ap, kk, 1, z, 0, 1 );
		}
		s = Blas1.dasum ( n, z, 0, 1 );
		for ( i = 1; i <= n; i++ )
		{
			z[i-1] = z[i-1] / s;
		}
		ynorm = 1.0;
		//
		//  Solve R' * V = Y.
		//
		for ( k = 1; k <= n; k++ )
		{
			z[k-1] = z[k-1] - Blas1.ddot ( k-1, ap, kk, 1, z, 0, 1 );
			kk = kk + k;

			if ( ap[kk-1] < Blas1.r8_abs ( z[k-1] ) )
			{
				s = ap[kk-1] / Blas1.r8_abs ( z[k-1] );
				for ( i = 1; i <= n; i++ )
				{
					z[i-1] = s * z[i-1];
				}
				ynorm = s * ynorm;
			}

			z[k-1] = z[k-1] / ap[kk-1];
		}
		s = 1.0 / Blas1.dasum ( n, z, 0, 1 );
		for ( i = 1; i <= n; i++ )
		{
			z[i-1] = s * z[i-1];
		}
		ynorm = s * ynorm;
		//
		//  Solve R * Z = V.
		//
		for ( k = n; 1 <= k; k-- )
		{
			if ( ap[kk-1] < Blas1.r8_abs ( z[k-1] ) )
			{
				s = ap[kk-1] / Blas1.r8_abs ( z[k-1] );
				for ( i = 1; i <= n; i++ )
				{
					z[i-1] = s * z[i-1];
				}
				ynorm = s * ynorm;
			}
			z[k-1] = z[k-1] / ap[kk-1];
			kk = kk - k;
			t = -z[k-1];
			Blas1.daxpy ( k-1, t, ap, kk, 1, z, 0, 1 );
		}
		//
		//  Make ZNORM = 1.0.
		//
		s = 1.0 / Blas1.dasum ( n, z, 0, 1 );
		for ( i = 1; i <= n; i++ )
		{
			z[i-1] = s * z[i-1];
		}
		ynorm = s * ynorm;

		if ( anorm != 0.0 )
		{
			rcond = ynorm / anorm;
		}
		else
		{
			rcond = 0.0;
		}

		return rcond;
	}
	//****************************************************************************80

	/*void dppdi ( double ap[], int n, double det[2], int job )
	{
		int i;
		int ii;
		int j;
		int j1;
		int jj;
		int k;
		int k1;
		int kj;
		int kk;
		double s;
		double t;
		//
		//  Compute the determinant.
		//
		if ( job / 10 != 0 )
		{
			det[0] = 1.0;
			det[1] = 0.0;
			s = 10.0;
			ii = 0;

			for ( i = 1; i <= n; i++ )
			{
				ii = ii + i;

				det[0] = det[0] * ap[ii-1] * ap[ii-1];

				if ( det[0] == 0.0 )
				{
					break;
				}

				while ( det[0] < 1.0 )
				{
					det[0] = det[0] * s;
					det[1] = det[1] - 1.0;
				}

				while ( s <= det[0] )
				{
					det[0] = det[0] / s;
					det[1] = det[1] + 1.0;
				}
			}
		}
		//
		//  Compute inverse(R).
		//
		if ( ( job % 10 ) != 0 )
		{
			kk = 0;

			for ( k = 1; k <= n; k++ )
			{
				k1 = kk + 1;
				kk = kk + k;
				ap[kk-1] = 1.0 / ap[kk-1];
				t = -ap[kk-1];
				Blas1.dscal  ( k-1, t, ap, k1-1, 1 );
				j1 = kk + 1;
				kj = kk + k;

				for ( j = k + 1; j <= n; j++ )
				{
					t = ap[kj-1];
					ap[kj-1] = 0.0;
					Blas1.daxpy ( k, t, ap, k1-1, 1, ap, j1-1, 1 );
					j1 = j1 + j;
					kj = kj + j;
				}
			}
			//
			//  Form inverse(R) * (inverse(R))'.
			//
			jj = 0;

			for ( j = 1; j <= n; j++ )
			{
				j1 = jj + 1;
				jj = jj + j;
				k1 = 1;
				kj = j1;

				for ( k = 1; k <= j-1; k++ )
				{
					t = ap[kj-1];
					Blas1.daxpy ( k, t, ap, j1-1, 1, ap, k1-1, 1 );
					k1 = k1 + k;
					kj = kj + 1;
				}
				t = ap[jj-1];
				Blas1.dscal  ( j, t, ap, j1-1, 1 );
			}
		}
		return;
	}*/
	//****************************************************************************80

	int dppfa ( double ap[], int n )
	{
		int info;
		int j;
		int jj;
		int k;
		int kj;
		int kk;
		double s;
		double t;

		info = 0;
		jj = 0;

		for ( j = 1; j <= n; j++ )
		{
			s = 0.0;
			kj = jj;
			kk = 0;

			for ( k = 1; k <= j-1; k++ )
			{
				kj = kj + 1;
				t = ap[kj-1] - Blas1.ddot ( k-1, ap, kk, 1, ap, jj, 1 );
				kk = kk + k;
				t = t / ap[kk-1];
				ap[kj-1] = t;
				s = s + t * t;
			}

			jj = jj + j;
			s = ap[jj-1] - s;

			if ( s <= 0.0 )
			{
				info = j;
				return info;
			}

			ap[jj-1] = Math.sqrt ( s );

		}
		return info;
	}
	//****************************************************************************80

	void dppsl ( double ap[], int n, double b[] )
	{
		int k;
		int kk;
		double t;

		kk = 0;

		for ( k = 1; k <= n; k++ )
		{
			t = Blas1.ddot ( k-1, ap, kk, 1, b, 0, 1 );
			kk = kk + k;
			b[k-1] = ( b[k-1] - t ) / ap[kk-1];
		}

		for ( k = n; 1 <= k; k-- )
		{
			b[k-1] = b[k-1] / ap[kk-1];
			kk = kk - k;
			t = -b[k-1];
			Blas1.daxpy ( k-1, t, ap, kk, 1, b, 0, 1 );
		}

		return;
	}
	//****************************************************************************80

	void dptsl ( int n, double d[], double e[], double b[] )
	{
		int k;
		int kbm1;
		int ke;
		int kf;
		int kp1;
		int nm1d2;
		double t1;
		double t2;
		//
		//  Check for 1 x 1 case.
		//
		if ( n == 1 )
		{
			b[0] = b[0] / d[0];
			return;
		}

		nm1d2 = ( n - 1 ) / 2;

		if ( 2 < n )
		{
			kbm1 = n - 1;
			//
			//  Zero top half of subdiagonal and bottom half of superdiagonal.
			//
			for ( k = 1; k <= nm1d2; k++ )
			{
				t1 = e[k-1] / d[k-1];
				d[k] = d[k] - t1 * e[k-1];
				b[k] = b[k] - t1 * b[k-1];
				t2 = e[kbm1-1] / d[kbm1];
				d[kbm1-1] = d[kbm1-1] - t2 * e[kbm1-1];
				b[kbm1-1] = b[kbm1-1] - t2 * b[kbm1];
				kbm1 = kbm1 - 1;
			}
		}

		kp1 = nm1d2 + 1;
		//
		//  Clean up for possible 2 x 2 block at center.
		//
		if ( ( n % 2 ) == 0 )
		{
			t1 = e[kp1-1] / d[kp1-1];
			d[kp1] = d[kp1] - t1 * e[kp1-1];
			b[kp1] = b[kp1] - t1 * b[kp1-1];
			kp1 = kp1 + 1;
		}
		//
		//  Back solve starting at the center, going towards the top and bottom.
		//
		b[kp1-1] = b[kp1-1] / d[kp1-1];

		if ( 2 < n )
		{
			k = kp1 - 1;
			ke = kp1 + nm1d2 - 1;

			for ( kf = kp1; kf <= ke; kf++ )
			{
				b[k-1] = ( b[k-1] - e[k-1] * b[k] ) / d[k-1];
				b[kf] = ( b[kf] - e[kf-1] * b[kf-1] ) / d[kf];
				k = k - 1;
			}
		}

		if ( ( n % 2 ) == 0 )
		{
			b[0] = ( b[0] - e[0] * b[1] ) / d[0];
		}

		return;
	}
	//****************************************************************************80

	void dqrdc ( double a[], int lda, int n, int p, double qraux[], int jpvt[], 
			double work[], int job )
	{
		int j;
		int jp;
		int l;
		int lup;
		int maxj;
		double maxnrm;
		double nrmxl;
		int pl;
		int pu;
		boolean swapj;
		double t;
		double tt;

		pl = 1;
		pu = 0;
		//
		//  If pivoting is requested, rearrange the columns.
		//
		if ( job != 0 )
		{
			for ( j = 1; j <= p; j++ )
			{
				swapj = ( 0 < jpvt[j-1] );

				if ( jpvt[j-1] < 0 )
				{
					jpvt[j-1] = -j;
				}
				else
				{
					jpvt[j-1] = j;
				}

				if ( swapj )
				{
					if ( j != pl )
					{
						Blas1.dswap  ( n, a, 0+(pl-1)*lda, 1, a, 0+(j-1), 1 );
					}
					jpvt[j-1] = jpvt[pl-1];
					jpvt[pl-1] = j;
					pl = pl + 1;
				}
			}
			pu = p;

			for ( j = p; 1 <= j; j-- )
			{
				if ( jpvt[j-1] < 0 )
				{
					jpvt[j-1] = -jpvt[j-1];

					if ( j != pu )
					{
						Blas1.dswap  ( n, a, 0+(pu-1)*lda, 1, a, 0+(j-1)*lda, 1 );
						jp = jpvt[pu-1];
						jpvt[pu-1] = jpvt[j-1];
						jpvt[j-1] = jp;
					}
					pu = pu - 1;
				}
			}
		}
		//
		//  Compute the norms of the free columns.
		//
		for ( j = pl; j <= pu; j++ )
		{
			qraux[j-1] = Blas1.dnrm2 ( n, a, 0+(j-1)*lda, 1 );
		}

		for ( j = pl; j <= pu; j++ )
		{
			work[j-1] = qraux[j-1];
		}
		//
		//  Perform the Householder reduction of A.
		//
		lup = Blas1.i4_min ( n, p );

		for ( l = 1; l <= lup; l++ )
		{
			//
			//  Bring the column of largest norm into the pivot position.
			//
			if ( pl <= l && l < pu )
			{
				maxnrm = 0.0;
				maxj = l;
				for ( j = l; j <= pu; j++ )
				{
					if ( maxnrm < qraux[j-1] )
					{
						maxnrm = qraux[j-1];
						maxj = j;
					}
				}

				if ( maxj != l )
				{
					Blas1.dswap  ( n, a, 0+(l-1)*lda, 1, a, 0+(maxj-1)*lda, 1 );
					qraux[maxj-1] = qraux[l-1];
					work[maxj-1] = work[l-1];
					jp = jpvt[maxj-1];
					jpvt[maxj-1] = jpvt[l-1];
					jpvt[l-1] = jp;
				}
			}
			//
			//  Compute the Householder transformation for column L.
			//
			qraux[l-1] = 0.0;

			if ( l != n )
			{
				nrmxl = Blas1.dnrm2 ( n-l+1, a, l-1+(l-1)*lda, 1 );

				if ( nrmxl != 0.0 )
				{
					if ( a[l-1+(l-1)*lda] != 0.0 )
					{
						nrmxl = nrmxl * Blas1.r8_sign ( a[l-1+(l-1)*lda] );
					}

					Blas1.dscal  ( n-l+1, 1.0 / nrmxl, a, l-1+(l-1)*lda, 1 );
					a[l-1+(l-1)*lda] = 1.0 + a[l-1+(l-1)*lda];
					//
					//  Apply the transformation to the remaining columns, updating the norms.
					//
					for ( j = l + 1; j <= p; j++ )
					{
						t = -Blas1.ddot ( n-l+1, a, l-1+(l-1)*lda, 1, a, l-1+(j-1)*lda, 1 ) 
								/ a[l-1+(l-1)*lda];
						Blas1.daxpy ( n-l+1, t, a, l-1+(l-1)*lda, 1, a, l-1+(j-1)*lda, 1 );

						if ( pl <= j && j <= pu )
						{
							if ( qraux[j-1] != 0.0 )
							{
								tt = 1.0 - Math.pow ( Blas1.r8_abs ( a[l-1+(j-1)*lda] ) / qraux[j-1], 2 );
								tt = Blas1.r8_max  ( tt, 0.0 );
								t = tt;
								tt = 1.0 + 0.05 * tt * Math.pow ( qraux[j-1] / work[j-1], 2 );

								if ( tt != 1.0 )
								{
									qraux[j-1] = qraux[j-1] * Math.sqrt ( t );
								}
								else
								{
									qraux[j-1] = Blas1.dnrm2 ( n-l, a, l+(j-1)*lda, 1 );
									work[j-1] = qraux[j-1];
								}
							}
						}
					}
					//
					//  Save the transformation.
					//
					qraux[l-1] = a[l-1+(l-1)*lda];
					a[l-1+(l-1)*lda] = -nrmxl;
				}
			}
		}
		return;
	}
	//****************************************************************************80

	int dqrsl ( double a[], int lda, int n, int k, double qraux[], double y[], 
			double qy[], double qty[], double b[], double rsd[], double ab[], int job )
	{
		boolean cab;
		boolean cb;
		boolean cqty;
		boolean cqy;
		boolean cr;
		int i;
		int info;
		int j;
		int jj;
		int ju;
		double t;
		double temp;
		//
		//  set info flag.
		//
		info = 0;
		//
		//  Determine what is to be computed.
		//
		cqy =  (   job / 10000          != 0 );
		cqty = ( ( job %  10000 )       != 0 );
		cb =   ( ( job %   1000 ) / 100 != 0 );
		cr =   ( ( job %    100 ) /  10 != 0 );
		cab =  ( ( job %     10 )       != 0 );

		ju = Blas1.i4_min ( k, n-1 );
		//
		//  Special action when N = 1.
		//
		if ( ju == 0 )
		{
			if ( cqy )
			{
				qy[0] = y[0];
			}

			if ( cqty )
			{
				qty[0] = y[0];
			}

			if ( cab )
			{
				ab[0] = y[0];
			}

			if ( cb )
			{
				if ( a[0+0*lda] == 0.0 )
				{
					info = 1;
				}
				else
				{
					b[0] = y[0] / a[0+0*lda];
				}
			}

			if ( cr )
			{
				rsd[0] = 0.0;
			}
			return info;
		}
		//
		//  Set up to compute QY or QTY.
		//
		if ( cqy )
		{
			for ( i = 1; i <= n; i++ )
			{
				qy[i-1] = y[i-1];
			}
		}

		if ( cqty )
		{
			for ( i = 1; i <= n; i++ )
			{
				qty[i-1] = y[i-1];
			}
		}
		//
		//  Compute QY.
		//
		if ( cqy )
		{
			for ( jj = 1; jj <= ju; jj++ )
			{
				j = ju - jj + 1;

				if ( qraux[j-1] != 0.0 )
				{
					temp = a[j-1+(j-1)*lda];
					a[j-1+(j-1)*lda] = qraux[j-1];
					t = -Blas1.ddot ( n-j+1, a, j-1+(j-1)*lda, 1, qy, j-1, 1 ) / a[j-1+(j-1)*lda];
					Blas1.daxpy ( n-j+1, t, a, j-1+(j-1)*lda, 1, qy, j-1, 1 );
					a[j-1+(j-1)*lda] = temp;
				}
			}
		}
		//
		//  Compute Q'*Y.
		//
		if ( cqty )
		{
			for ( j = 1; j <= ju; j++ )
			{
				if ( qraux[j-1] != 0.0 )
				{
					temp = a[j-1+(j-1)*lda];
					a[j-1+(j-1)*lda] = qraux[j-1];
					t = -Blas1.ddot ( n-j+1, a, j-1+(j-1)*lda, 1, qty, j-1, 1 ) / a[j-1+(j-1)*lda];
					Blas1.daxpy ( n-j+1, t, a, j-1+(j-1)*lda, 1, qty, j-1, 1 );
					a[j-1+(j-1)*lda] = temp;
				}
			}
		}
		//
		//  Set up to compute B, RSD, or AB.
		//
		if ( cb )
		{
			for ( i = 1; i <= k; i++ )
			{
				b[i-1] = qty[i-1];
			}
		}

		if ( cab )
		{
			for ( i = 1; i <= k; i++ )
			{
				ab[i-1] = qty[i-1];
			}
		}

		if ( cr && k < n )
		{
			for ( i = k+1; i <= n; i++ )
			{
				rsd[i-1] = qty[i-1];
			}
		}

		if ( cab && k+1 <= n )
		{
			for ( i = k+1; i <= n; i++ )
			{
				ab[i-1] = 0.0;
			}
		}

		if ( cr )
		{
			for ( i = 1; i <= k; i++ )
			{
				rsd[i-1] = 0.0;
			}
		}
		//
		//  Compute B.
		//
		if ( cb )
		{
			for ( jj = 1; jj <= k; jj++ )
			{
				j = k - jj + 1;

				if ( a[j-1+(j-1)*lda] == 0.0 )
				{
					info = j;
					break;
				}

				b[j-1] = b[j-1] / a[j-1+(j-1)*lda];

				if ( j != 1 )
				{
					t = -b[j-1];
					Blas1.daxpy ( j-1, t, a, 0+(j-1)*lda, 1, b, 0, 1 );
				}
			}
		}
		//
		//  Compute RSD or AB as required.
		//
		if ( cr || cab )
		{
			for ( jj = 1; jj <= ju; jj++ )
			{
				j = ju - jj + 1;

				if ( qraux[j-1] != 0.0 )
				{
					temp = a[j-1+(j-1)*lda];
					a[j-1+(j-1)*lda] = qraux[j-1];

					if ( cr )
					{
						t = -Blas1.ddot ( n-j+1, a, j-1+(j-1)*lda, 1, rsd, j-1, 1 ) 
								/ a[j-1+(j-1)*lda];
						Blas1.daxpy ( n-j+1, t, a, j-1+(j-1)*lda, 1, rsd, j-1, 1 );
					}

					if ( cab )
					{
						t = -Blas1.ddot ( n-j+1, a, j-1+(j-1)*lda, 1, ab, j-1, 1 ) 
								/ a[j-1+(j-1)*lda];
						Blas1.daxpy ( n-j+1, t, a, j-1+(j-1)*lda, 1, ab, j-1, 1 );
					}
					a[j-1+(j-1)*lda] = temp;
				}
			}
		}

		return info;
	}
	//****************************************************************************80

	double dsico ( double a[], int lda, int n, int kpvt[], double z[] )
	{
		double ak;
		double akm1;
		double anorm;
		double bk;
		double bkm1;
		double denom;
		double ek;
		int i;
		int info;
		int j;
		int k;
		int kp;
		int kps;
		int ks;
		double rcond;
		double s;
		double t;
		double ynorm;
		//
		//  Find the norm of A, using only entries in the upper half of the matrix.
		//
		for ( j = 1; j <= n; j++ )
		{
			z[j-1] = Blas1.dasum ( j, a, 0+(j-1)*lda, 1 );
			for ( i = 1; i <= j-1; i++ )
			{
				z[i-1] = z[i-1] + Blas1.r8_abs ( a[i-1+(j-1)*lda] );
			}
		}

		anorm = 0.0;
		for ( i = 1; i <= n; i++ )
		{
			anorm = Blas1.r8_max  ( anorm, z[i-1] );
		}
		//
		//  Factor.
		//
		info = dsifa ( a, lda, n, kpvt );
		//
		//  RCOND = 1/(norm(A)*(estimate of norm(inverse(A)))).
		//
		//  Estimate = norm(Z)/norm(Y) where A*Z = Y and A*Y = E.
		//
		//  The components of E are chosen to cause maximum local
		//  growth in the elements of W where U*D*W = E.
		//
		//  The vectors are frequently rescaled to avoid overflow.
		//
		//  Solve U * D * W = E.
		//
		ek = 1.0;
		for ( i = 1; i <= n; i++ )
		{
			z[i-1] = 0.0;
		}
		k = n;

		while ( k != 0 )
		{
			if ( kpvt[k-1] < 0 )
			{
				ks = 2;
			}
			else
			{
				ks = 1;
			}
			kp = Math.abs ( kpvt[k-1] );
			kps = k + 1 - ks;

			if ( kp != kps )
			{
				t = z[kps-1];
				z[kps-1] = z[kp-1];
				z[kp-1] = t;
			}

			if ( z[k-1] != 0.0 )
			{
				ek = ek * Blas1.r8_sign ( z[k-1] );
			}

			z[k-1] = z[k-1] + ek;
			Blas1.daxpy ( k-ks, z[k-2], a, 0+(k-1)*lda, 1, z, 0, 1 );

			if ( ks != 1 )
			{
				if ( z[k-2] != 0.0 )
				{
					ek = ek * Blas1.r8_sign ( z[k-2] );
				}
				z[k-2] = z[k-2] + ek;
				Blas1.daxpy ( k-ks, z[k-2], a, 0+(k-2)*lda, 1, z, 0, 1 );
			}

			if ( ks != 2 )
			{
				if ( Blas1.r8_abs ( a[k-1+(k-1)*lda] ) < Blas1.r8_abs ( z[k-1] ) )
				{
					s = Blas1.r8_abs ( a[k-1+(k-1)*lda] ) / Blas1.r8_abs ( z[k-1] );
					for ( i = 1; i <= n; i++ )
					{
						z[i-1] = s * z[i-1];
					}
					ek = s * ek;
				}

				if ( a[k-1+(k-1)*lda] != 0.0 )
				{
					z[k-1] = z[k-1] / a[k-1+(k-1)*lda];
				}
				else
				{
					z[k-1] = 1.0;
				}
			}
			else
			{
				ak = a[k-1+(k-1)*lda] / a[k-2+(k-1)*lda];
				akm1 = a[k-2+(k-2)*lda] / a[k-2+(k-1)*lda];
				bk = z[k-1] / a[k-2+(k-1)*lda];
				bkm1 = z[k-2] / a[k-2+(k-1)*lda];
				denom = ak * akm1 - 1.0;
				z[k-1] = ( akm1 * bk - bkm1 ) / denom;
				z[k-2] = ( ak * bkm1 - bk ) / denom;
			}
			k = k - ks;
		}
		s = Blas1.dasum ( n, z, 0, 1 );
		for ( i = 1; i <= n; i++ )
		{
			z[i-1] = z[i-1] / s;
		}
		//
		//  Solve U' * Y = W.
		//
		k = 1;

		while ( k <= n )
		{
			if ( kpvt[k-1] < 0 )
			{
				ks = 2;
			}
			else
			{
				ks = 1;
			}

			if ( k != 1 )
			{
				z[k-1] = z[k-1] + Blas1.ddot ( k-1, a, 0+(k-1)*lda, 1, z, 0, 1 );

				if ( ks == 2 )
				{
					z[k] = z[k] + Blas1.ddot ( k-1, a, 0+k*lda, 1, z, 0, 1 );
				}

				kp = Math.abs ( kpvt[k-1] );

				if ( kp != k )
				{
					t = z[k-1];
					z[k-1] = z[kp-1];
					z[kp-1] = t;
				}
			}
			k = k + ks;
		}
		s = Blas1.dasum ( n, z, 0, 1 );
		for ( i = 1; i <= n; i++ )
		{
			z[i-1] = z[i-1] / s;
		}
		ynorm = 1.0;
		//
		//  Solve U * D * V = Y.
		//
		k = n;

		while ( k != 0 )
		{
			if ( kpvt[k-1] < 0 )
			{
				ks = 2;
			}
			else
			{
				ks = 1;
			}

			if ( k != ks )
			{
				kp = Math.abs ( kpvt[k-1] );
				kps = k + 1 - ks;

				if ( kp != kps )
				{
					t = z[kps-1];
					z[kps-1] = z[kp-1];
					z[kp-1] = t;
				}

				Blas1.daxpy ( k-ks, z[k-1], a, 0+(k-1)*lda, 1, z, 0, 1 );

				if ( ks == 2 )
				{
					Blas1.daxpy ( k-ks, z[k-2], a, 0+(k-2)*lda, 1, z, 0, 1 );
				}
			}

			if ( ks != 2 )
			{
				if ( Blas1.r8_abs ( a[k-1+(k-1)*lda] ) < Blas1.r8_abs ( z[k-1] ) )
				{
					s = Blas1.r8_abs ( a[k-1+(k-1)*lda] ) / Blas1.r8_abs ( z[k-1] );
					for ( i = 1; i <= n; i++ )
					{
						z[i-1] = s * z[i-1];
					}
					ynorm = s * ynorm;
				}

				if ( a[k-1+(k-1)*lda] != 0.0 )
				{
					z[k-1] = z[k-1] / a[k-1+(k-1)*lda];
				}
				else
				{
					z[k-1] = 1.0;
				}
			}
			else
			{
				ak = a[k-1+(k-1)*lda] / a[k-2+(k-1)*lda];
				akm1 = a[k-2+(k-2)*lda] / a[k-2+(k-1)*lda];
				bk = z[k-1] / a[k-2+(k-1)*lda];
				bkm1 = z[k-2] / a[k-2+(k-1)*lda];
				denom = ak * akm1 - 1.0;
				z[k-1] = ( akm1 * bk - bkm1 ) / denom;
				z[k-2] = ( ak * bkm1 - bk ) / denom;
			}
			k = k - ks;
		}

		s = 1.0 / Blas1.dasum ( n, z, 0, 1 );
		for ( i = 1; i <= n; i++ )
		{
			z[i-1] = s * z[i-1];
		}
		ynorm = s * ynorm;
		//
		//  Solve U' * Z = V.
		//
		k = 1;

		while ( k <= n )
		{
			if ( kpvt[k-1] < 0 )
			{
				ks = 2;
			}
			else
			{
				ks = 1;
			}

			if ( k != 1 )
			{
				z[k-1] = z[k-1] + Blas1.ddot ( k-1, a, 0+(k-1)*lda, 1, z, 0, 1 );
				if ( ks == 2 )
				{
					z[k] = z[k] + Blas1.ddot ( k-1, a, 0+k*lda, 1, z, 0, 1 );
				}
				kp = Math.abs ( kpvt[k-1] );

				if ( kp != k )
				{
					t = z[k-1];
					z[k-1] = z[kp-1];
					z[kp-1] = t;
				}
			}
			k = k + ks;
		}
		//
		//  Make ZNORM = 1.0.
		//
		s = 1.0 / Blas1.dasum ( n, z, 0, 1 );
		for ( i = 1; i <= n; i++ )
		{
			z[i-1] = s * z[i-1];
		}
		ynorm = s * ynorm;

		if ( anorm != 0.0 )
		{
			rcond = ynorm / anorm;
		}
		else
		{
			rcond = 0.0;
		}

		return rcond;
	}
	//****************************************************************************80

/*	void dsidi ( double a[], int lda, int n, int kpvt[], double det[2], 
			int inert[3], double work[], int job )
	{
		double ak;
		double akkp1;
		double akp1;
		double d;
		boolean dodet;
		boolean doert;
		boolean doinv;
		int j;
		int jb;
		int k;
		int ks;
		int kstep;
		double t;
		double temp;

		doinv = ( job %   10 )       != 0;
		dodet = ( job %  100 ) /  10 != 0;
		doert = ( job % 1000 ) / 100 != 0;

		if ( dodet || doert )
		{
			if ( doert )
			{
				inert[0] = 0;
				inert[1] = 0;
				inert[2] = 0;
			}

			if ( dodet )
			{
				det[0] = 1.0;
				det[1] = 0.0;
			}

			t = 0.0;

			for ( k = 1; k <= n; k++ )
			{
				d = a[k-1+(k-1)*lda];
				//
				//  2 by 2 block.
				//
				//  use det (d  s)  =  (d/t * c - t) * t,  t = abs ( s )
				//	          (s  c)
				//  to avoid underflow/overflow troubles.
				//
				//  Take two passes through scaling.  Use T for flag.
				//
				if ( kpvt[k-1] <= 0 )
				{
					if ( t == 0.0 )
					{
						t = Blas1.r8_abs ( a[k-1+k*lda] );
						d = ( d / t ) * a[k+k*lda] - t;
					}
					else
					{
						d = t;
						t = 0.0;
					}
				}

				if ( doert )
				{
					if ( 0.0 < d )
					{
						inert[0] = inert[0] + 1;
					}
					else if ( d < 0.0 )
					{
						inert[1] = inert[1] + 1;
					}
					else if ( d == 0.0 )
					{
						inert[2] = inert[2] + 1;
					}
				}

				if ( dodet )
				{
					det[0] = det[0] * d;

					if ( det[0] != 0.0 )
					{
						while ( Blas1.r8_abs ( det[0] ) < 1.0 )
						{
							det[0] = det[0] * 10.0;
							det[1] = det[1] - 1.0;
						}

						while ( 10.0 <= Blas1.r8_abs ( det[0] ) )
						{
							det[0] = det[0] / 10.0;
							det[1] = det[1] + 1.0;
						}
					}
				}
			}
		}
		//
		//  Compute inverse(A).
		//
		if ( doinv )
		{
			k = 1;

			while ( k <= n )
			{
				if ( 0 <= kpvt[k-1] )
				{
					//
					//  1 by 1.
					//
					a[k-1+(k-1)*lda] = 1.0 / a[k-1+(k-1)*lda];

					if ( 2 <= k )
					{
						Blas1.dcopy ( k-1, a, 0+(k-1)*lda, 1, work, 1 );

						for ( j = 1; j <= k-1; j++ )
						{
							a[j-1+(k-1)*lda] = Blas1.ddot ( j, a, 0+(j-1)*lda, 1, work, 1 );
							Blas1.daxpy ( j-1, work[j-1], a, 0+(j-1)*lda, 1, a, 0+(k-1)*lda, 1 );
						}
						a[k-1+(k-1)*lda] = a[k-1+(k-1)*lda] 
								+ Blas1.ddot ( k-1, work, 1, a, 0+(k-1)*lda, 1 );
					}
					kstep = 1;
				}
				//
				//  2 by 2.
				//
				else
				{
					t = Blas1.r8_abs ( a[k-1+k*lda] );
					ak = a[k-1+(k-1)*lda] / t;
					akp1 = a[k+k*lda] / t;
					akkp1 = a[k-1+k*lda] / t;
					d = t * ( ak * akp1 - 1.0 );
					a[k-1+(k-1)*lda] = akp1 / d;
					a[k+k*lda] = ak / d;
					a[k-1+k*lda] = -akkp1 / d;

					if ( 2 <= k )
					{
						Blas1.dcopy ( k-1, a, 0+k*lda, 1, work, 1 );

						for ( j = 1; j <= k-1; j++ )
						{
							a[j-1+k*lda] = Blas1.ddot ( j, a, 0+(j-1)*lda, 1, work, 1 );
							Blas1.daxpy ( j-1, work[j-1], a, 0+(j-1)*lda, 1, a, 0+k*lda, 1 );
						}
						a[k+k*lda] = a[k+k*lda] + Blas1.ddot ( k-1, work, 1, a, 0+k*lda, 1 );
						a[k-1+k*lda] = a[k-1+k*lda] 
								+ Blas1.ddot ( k-1, a, 0+(k-1)*lda, 1, a, 0+k*lda, 1 );
						Blas1.dcopy ( k-1, a, 0+(k-1)*lda, 1, work, 1 );

						for ( j = 1; j <= k-1; j++ )
						{
							a[j-1+(k-1)*lda] = Blas1.ddot ( j, a, 0+(j-1)*lda, 1, work, 1 );
							Blas1.daxpy ( j-1, work[j-1], a, 0+(j-1)*lda, 1, a, 0+(k-1)*lda, 1 );
						}
						a[k-1+(k-1)*lda] = a[k-1+(k-1)*lda] 
								+ Blas1.ddot ( k-1, work, 1, a, 0+(k-1)*lda, 1 );
					}
					kstep = 2;
				}
				//
				//  Swap.
				//
				ks = abs ( kpvt[k-1] );

				if ( ks != k )
				{
					Blas1.dswap  ( ks, a, 0+(ks-1)*lda, 1, a, 0+(k-1)*lda, 1 );

					for ( jb = ks; jb <= k; jb++ )
					{
						j = k + ks - jb;
						temp = a[j-1+(k-1)*lda];
						a[j-1+(k-1)*lda] = a[ks-1+(j-1)*lda];
						a[ks-1+(j-1)*lda] = temp;
					}

					if ( kstep != 1 )
					{
						temp = a[ks-1+k*lda];
						a[ks-1+k*lda] = a[k-1+k*lda];
						a[k-1+k*lda] = temp;
					}
				}
				k = k + kstep;
			}
		}
		return;
	}*/
	//****************************************************************************80

	int dsifa ( double a[], int lda, int n, int kpvt[] )
	{
		double absakk;
		double ak;
		double akm1;
		double alpha;
		double bk;
		double bkm1;
		double colmax;
		double denom;
		int imax;
		int imaxp1;
		int info;
		int j;
		int jj;
		int jmax;
		int k;
		int km1;
		int kstep;
		double mulk;
		double mulkm1;
		double rowmax;
		boolean swap;
		double t;
		//
		//  ALPHA is used in choosing pivot block size.
		//
		alpha = ( 1.0 + Math.sqrt ( 17.0 ) ) / 8.0;

		info = 0;
		//
		//  Main loop on K, which goes from N to 1.
		//
		k = n;

		while ( 0 < k )
		{
			if ( k == 1 )
			{
				kpvt[0] = 1;
				if ( a[0+0*lda] == 0.0 )
				{
					info = 1;
				}
				return info;
			}
			//
			//  This section of code determines the kind of
			//  elimination to be performed.  When it is completed,
			//  KSTEP will be set to the size of the pivot block, and
			//  SWAP will be set to .true. if an interchange is required.
			//
			km1 = k - 1;
			absakk = Blas1.r8_abs ( a[k-1+(k-1)*lda] );
			//
			//  Determine the largest off-diagonal element in column K.
			//
			imax = Blas1.idamax  ( k-1, a, 0+(k-1)*lda, 1 );
			colmax = Blas1.r8_abs ( a[imax-1+(k-1)*lda] );

			if ( alpha * colmax <= absakk )
			{
				kstep = 1;
				swap = false;
			}
			//
			//  Determine the largest off-diagonal element in row IMAX.
			//
			else
			{
				rowmax = 0.0;
				imaxp1 = imax + 1;
				for ( j = imaxp1; j <= k; j++ )
				{
					rowmax = Blas1.r8_max  ( rowmax, Blas1.r8_abs ( a[imax-1+(j-1)*lda] ) );
				}

				if ( imax != 1 )
				{
					jmax = Blas1.idamax  ( imax-1, a, 0+(imax-1)*lda, 1 );
					rowmax = Blas1.r8_max  ( rowmax, Blas1.r8_abs ( a[jmax-1+(imax-1)*lda] ) );
				}

				if ( alpha * rowmax <= Blas1.r8_abs ( a[imax-1+(imax-1)*lda] ) )
				{
					kstep = 1;
					swap = true;
				}
				else if ( alpha * colmax * ( colmax / rowmax ) <= absakk ) 
				{
					kstep = 1;
					swap = false;
				}
				else
				{
					kstep = 2;
					swap = ( imax != k-1 );
				}
			}
			//
			//  Column K is zero.
			//  Set INFO and iterate the loop.
			//
			if ( Blas1.r8_max  ( absakk, colmax ) == 0.0 )
			{
				kpvt[k-1] = k;
				info = k;
			}
			//
			//  1 x 1 pivot block.
			//
			//  Perform an interchange.
			//
			else if ( kstep != 2 )
			{
				if ( swap )
				{
					Blas1.dswap  ( imax, a, 0+(imax-1)*lda, 1, a, 0+(k-1)*lda, 1 );

					for ( jj = imax; jj <= k; jj++ )
					{
						j = k + imax - jj;
						t = a[j-1+(k-1)*lda];
						a[j-1+(k-1)*lda] = a[imax-1+(j-1)*lda];
						a[imax-1+(j-1)*lda] = t;
					}
				}
				//
				//  Perform the elimination.
				//
				for ( jj = 1; jj <= k-1; jj++ )
				{
					j = k - jj;
					mulk = -a[j-1+(k-1)*lda] / a[k-1+(k-1)*lda];
					t = mulk;
					Blas1.daxpy ( j, t, a, 0+(k-1)*lda, 1, a, 0+(j-1)*lda, 1 );
					a[j-1+(k-1)*lda] = mulk;
				}
				//
				//  Set the pivot array.
				//
				if ( swap )
				{
					kpvt[k-1] = imax;
				}
				else
				{
					kpvt[k-1] = k;
				}
			}
			//
			//  2 x 2 pivot block.
			//
			//  Perform an interchange.
			//
			else
			{
				if ( swap )
				{
					Blas1.dswap  ( imax, a, 0+(imax-1)*lda, 1, a, 0+(k-2)*lda, 1 );

					for ( jj = imax; jj <= k-1; jj++ )
					{
						j = k-1 + imax - jj;
						t = a[j-1+(k-1)*lda];
						a[j-1+(k-1)*lda] = a[imax-1+(j-1)*lda];
						a[imax-1+(j-1)*lda] = t;
					}

					t = a[k-2+(k-1)*lda];
					a[k-2+(k-1)*lda] = a[imax-1+(k-1)*lda];
					a[imax-1+(k-1)*lda] = t;
				}
				//
				//  Perform the elimination.
				//
				if ( k-2 != 0 )
				{
					ak = a[k-1+(k-1)*lda] / a[k-2+(k-1)*lda];
					akm1 = a[k-2+(k-2)*lda] / a[k-2+(k-1)*lda];
					denom = 1.0 - ak * akm1;

					for ( jj = 1; jj <= k-2; jj++ )
					{
						j = k-1 - jj;
						bk = a[j-1+(k-1)*lda] / a[k-2+(k-1)*lda];
						bkm1 = a[j-1+(k-2)*lda] / a[k-2+(k-1)*lda];
						mulk = ( akm1 * bk - bkm1 ) / denom;
						mulkm1 = ( ak * bkm1 - bk ) / denom;
						t = mulk;
						Blas1.daxpy ( j, t, a, 0+(k-1)*lda, 1, a, 0+(j-1)*lda, 1 );
						t = mulkm1;
						Blas1.daxpy ( j, t, a, 0+(k-2)*lda, 1, a, 0+(j-1)*lda, 1 );
						a[j-1+(k-1)*lda] = mulk;
						a[j-1+(k-2)*lda] = mulkm1;
					}
				}
				//
				//  Set the pivot array.
				//
				if ( swap )
				{
					kpvt[k-1] = -imax;
				}
				else
				{
					kpvt[k-1] = 1 - k;
				}
				kpvt[k-2] = kpvt[k-1];
			}
			k = k - kstep;
		}
		return info;
	}
	//****************************************************************************80

	void dsisl ( double a[], int lda, int n, int kpvt[], double b[] )
	{
		double ak;
		double akm1;
		double bk;
		double bkm1;
		double denom;
		int k;
		int kp;
		double temp;
		//
		//  Loop backward applying the transformations and D inverse to B.
		//
		k = n;

		while ( 0 < k )
		{
			if ( 0 <= kpvt[k-1] )
			{
				//
				//  1 x 1 pivot block.
				//
				if ( k != 1 )
				{
					kp = kpvt[k-1];
					//
					//  Interchange.
					//
					if ( kp != k )
					{
						temp = b[k-1];
						b[k-1] = b[kp-1];
						b[kp-1] = temp;
					}
					//
					//  Apply the transformation.
					//
					Blas1.daxpy ( k-1, b[k-1], a, 0+(k-1)*lda, 1, b, 0, 1 );
				}
				//
				//  Apply D inverse.
				//
				b[k-1] = b[k-1] / a[k-1+(k-1)*lda];
				k = k - 1;
			}
			else
			{
				//
				//  2 x 2 pivot block.
				//
				if ( k != 2 )
				{
					kp = Math.abs ( kpvt[k-1] );
					//
					//  Interchange.
					//
					if ( kp != k-1 )
					{
						temp = b[k-2];
						b[k-2] = b[kp-1];
						b[kp-1] = temp;
					}
					//
					//  Apply the transformation.
					//
					Blas1.daxpy ( k-2, b[k-1], a, 0+(k-1)*lda, 1, b, 0, 1 );
					Blas1.daxpy ( k-2, b[k-2], a, 0+(k-2)*lda, 1, b, 0, 1 );
				}
				//
				//  Apply D inverse.
				//
				ak = a[k-1+(k-1)*lda] / a[k-2+(k-1)*lda];
				akm1 = a[k-2+(k-2)*lda] / a[k-2+(k-1)*lda];
				bk = b[k-1] / a[k-2+(k-1)*lda];
				bkm1 = b[k-2] / a[k-2+(k-1)*lda];
				denom = ak * akm1 - 1.0;
				b[k-1] = ( akm1 * bk - bkm1 ) / denom;
				b[k-2] = ( ak * bkm1 - bk ) / denom;
				k = k - 2;
			}
		}
		//
		//  Loop forward applying the transformations.
		//
		k = 1;

		while ( k <= n )
		{
			if ( 0 <= kpvt[k-1] )
			{
				//
				//  1 x 1 pivot block.
				//
				if ( k != 1 )
				{
					//
					//  Apply the transformation.
					//
					b[k-1] = b[k-1] + Blas1.ddot ( k-1, a, 0+(k-1)*lda, 1, b, 0, 1 );
					kp = kpvt[k-1];
					//
					//  Interchange.
					//
					if ( kp != k )
					{
						temp = b[k-1];
						b[k-1] = b[kp-1];
						b[kp-1] = temp;
					}
				}
				k = k + 1;
			}
			else
			{
				//
				//  2 x 2 pivot block.
				//
				if ( k != 1 )
				{
					//
					//  Apply the transformation.
					//
					b[k-1] = b[k-1] + Blas1.ddot ( k-1, a, 0+(k-1)*lda, 1, b, 0, 1 );
					b[k] = b[k] + Blas1.ddot ( k-1, a, 0+k*lda, 1, b, 0, 1 );
					kp = Math.abs ( kpvt[k-1] );
					//
					//  Interchange.
					//
					if ( kp != k )
					{
						temp = b[k-1];
						b[k-1] = b[kp-1];
						b[kp-1] = temp;
					}
				}
				k = k + 2;
			}
		}
		return;
	}
	//****************************************************************************80

	double dspco ( double ap[], int n, int kpvt[], double z[] )
	{
		double ak;
		double akm1;
		double anorm;
		double bk;
		double bkm1;
		double denom;
		double ek;
		int i;
		int ij;
		int ik;
		int ikm1;
		int ikp1;
		int info;
		int j;
		int j1;
		int k;
		int kk;
		int km1k;
		int km1km1;
		int kp;
		int kps;
		int ks;
		double rcond;
		double s;
		double t;
		double ynorm;
		//
		//  Find norm of A using only upper half.
		//
		j1 = 1;
		for ( j = 1; j <= n; j++ )
		{
			z[j-1] = Blas1.dasum ( j, ap, j1-1, 1 );
			ij = j1;
			j1 = j1 + j;
			for ( i = 1; i <= j-1; i++ )
			{
				z[i-1] = z[i-1] + Blas1.r8_abs ( ap[ij-1] );
				ij = ij + 1;
			}
		}
		anorm = 0.0;
		for ( i = 1; i <= n; i++ )
		{
			anorm = Blas1.r8_max  ( anorm, z[i-1] );
		}
		//
		//  Factor.
		//
		info = dspfa ( ap, n, kpvt );
		//
		//  RCOND = 1/(norm(A)*(estimate of norm(inverse(A)))).
		//
		//  Estimate = norm(Z)/norm(Y) where A*Z = Y and A*Y = E.
		//
		//  The components of E are chosen to cause maximum local
		//  growth in the elements of W where U*D*W = E.
		//
		//  The vectors are frequently rescaled to avoid overflow.
		//
		//  Solve U * D * W = E.
		//
		ek = 1.0;
		for ( i = 1; i <= n; i++ )
		{
			z[i-1] = 0.0;
		}

		k = n;
		ik = ( n * ( n - 1 ) ) / 2;

		while ( k != 0 )
		{
			kk = ik + k;
			ikm1 = ik - ( k - 1 );

			if ( kpvt[k-1] < 0 )
			{
				ks = 2;
			}
			else
			{
				ks = 1;
			}

			kp = Math.abs ( kpvt[k-1] );
			kps = k + 1 - ks;

			if ( kp != kps )
			{
				t = z[kps-1];
				z[kps-1] = z[kp-1];
				z[kp-1] = t;
			}

			if ( z[k-1] != 0.0 )
			{
				ek = ek * Blas1.r8_sign ( z[k-1] );
			}

			z[k-1] = z[k-1] + ek;
			Blas1.daxpy ( k-ks, z[k-1], ap, ik, 1, z, 0, 1 );

			if ( ks != 1 )
			{
				if ( z[k-2] != 0.0 )
				{
					ek = ek * Blas1.r8_sign ( z[k-2] );
				}
				z[k-2] = z[k-2] + ek;
				Blas1.daxpy ( k-ks, z[k-2], ap, ikm1, 1, z, 0, 1 );
			}

			if ( ks != 2 )
			{
				if ( Blas1.r8_abs ( ap[kk-1] ) < Blas1.r8_abs ( z[k-1] ) )
				{
					s = Blas1.r8_abs ( ap[kk-1] ) / Blas1.r8_abs ( z[k-1] );
					for ( i = 1; i <= n; i++ )
					{
						z[i-1] = s * z[i-1];
					}
					ek = s * ek;
				}

				if ( ap[kk-1] != 0.0 )
				{
					z[k-1] = z[k-1] / ap[kk-1];
				}
				else
				{
					z[k-1] = 1.0;
				}
			}
			else
			{
				km1k = ik + k - 1;
				km1km1 = ikm1 + k - 1;
				ak = ap[kk-1] / ap[km1k-1];
				akm1 = ap[km1km1-1] / ap[km1k-1];
				bk = z[k-1] / ap[km1k-1];
				bkm1 = z[k-2] / ap[km1k-1];
				denom = ak * akm1 - 1.0;
				z[k-1] = ( akm1 * bk - bkm1 ) / denom;
				z[k-2] = ( ak * bkm1 - bk ) / denom;
			}
			k = k - ks;
			ik = ik - k;
			if ( ks == 2 )
			{
				ik = ik - ( k + 1 );
			}
		}

		s = Blas1.dasum ( n, z, 0, 1 );
		for ( i = 1; i <= n; i++ )
		{
			z[i-1] = z[i-1] / s;
		}
		//
		//  Solve U' * Y = W.
		//
		k = 1;
		ik = 0;

		while ( k <= n )
		{
			if ( kpvt[k-1] < 0 )
			{
				ks = 2;
			}
			else
			{
				ks = 1;
			}

			if ( k != 1 )
			{
				z[k-1] = z[k-1] + Blas1.ddot ( k-1, ap, ik, 1, z, 0, 1 );
				ikp1 = ik + k;

				if ( ks == 2 )
				{
					z[k] = z[k] + Blas1.ddot ( k-1, ap, ikp1, 1, z, 0, 1 );
				}

				kp = Math.abs ( kpvt[k-1] );

				if ( kp != k )
				{
					t = z[k-1];
					z[k-1] = z[kp-1];
					z[kp-1] = t;
				}
			}

			ik = ik + k;
			if ( ks == 2 )
			{
				ik = ik + ( k + 1 );
			}
			k = k + ks;
		}
		for ( i = 1; i <= n; i++ )
		{
			z[i-1] = s * z[i-1];
		}
		ynorm = 1.0;
		//
		//  Solve U * D * V = Y.
		//
		k = n;

		ik = ( n * ( n - 1 ) ) / 2;

		while ( 0 < k )
		{
			kk = ik + k;
			ikm1 = ik - ( k - 1 );

			if ( kpvt[k-1] < 0 )
			{
				ks = 2;
			}
			else
			{
				ks = 1;
			}

			if ( k != ks )
			{
				kp = Math.abs ( kpvt[k-1] );
				kps = k + 1 - ks;

				if ( kp != kps )
				{
					t = z[kps-1];
					z[kps-1] = z[kp-1];
					z[kp-1] = t;
				}

				Blas1.daxpy ( k-ks, z[k-1], ap, ik, 1, z, 0, 1 );

				if ( ks == 2 )
				{
					Blas1.daxpy ( k-ks, z[k-2], ap, ikm1, 1, z, 0, 1 );
				}
			}

			if ( ks != 2 )
			{
				if ( Blas1.r8_abs ( ap[kk-1] ) < Blas1.r8_abs ( z[k-1] ) )
				{
					s = Blas1.r8_abs ( ap[kk-1] ) / Blas1.r8_abs ( z[k-1] );
					for ( i = 1; i <= n; i++ )
					{
						z[i-1] = s * z[i-1];
					}
					ynorm = s * ynorm;
				}

				if ( ap[kk-1] != 0.0 )
				{
					z[k-1] = z[k-1] / ap[kk-1];
				}
				else
				{
					z[k-1] = 1.0;
				}
			}
			else
			{
				km1k = ik + k - 1;
				km1km1 = ikm1 + k - 1;
				ak = ap[kk-1] / ap[km1k-1];
				akm1 = ap[km1km1-1] / ap[km1k-1];
				bk = z[k-1] / ap[km1k-1];
				bkm1 = z[k-2] / ap[km1k-1];
				denom = ak * akm1 - 1.0;
				z[k-1] = ( akm1 * bk - bkm1 ) / denom;
				z[k-2] = ( ak * bkm1 - bk ) / denom;
			}
			k = k - ks;
			ik = ik - k;
			if ( ks == 2 )
			{
				ik = ik - ( k + 1 );
			}
		}
		s = 1.0 / Blas1.dasum ( n, z, 0, 1 );
		for ( i = 1; i <= n; i++ )
		{
			z[i-1] = s * z[i-1];
		}
		ynorm = s * ynorm;
		//
		//  Solve U' * Z = V.
		//
		k = 1;
		ik = 0;

		while ( k <= n )
		{
			if ( kpvt[k-1] < 0 )
			{
				ks = 2;
			}
			else
			{
				ks = 1;
			}

			if ( k != 1 )
			{
				z[k-1] = z[k-1] + Blas1.ddot ( k-1, ap, ik, 1, z, 0, 1 );
				ikp1 = ik + k;

				if ( ks == 2 )
				{
					z[k] = z[k] + Blas1.ddot ( k-1, ap, ikp1, 1, z, 0, 1 );
				}

				kp = Math.abs ( kpvt[k-1] );

				if ( kp != k )
				{
					t = z[k-1];
					z[k-1] = z[kp-1];
					z[kp-1] = t;
				}
			}

			ik = ik + k;
			if ( ks == 2 )
			{
				ik = ik + ( k + 1 );
			}
			k = k + ks;
		}
		//
		//  Make ZNORM = 1.0.
		//
		s = 1.0 / Blas1.dasum ( n, z, 0, 1 );
		for ( i = 1; i <= n; i++ )
		{
			z[i-1] = s * z[i-1];
		}
		ynorm = s * ynorm;

		if ( anorm != 0.0 )
		{
			rcond = ynorm / anorm;
		}
		else
		{
			rcond = 0.0;
		}

		return rcond;
	}
	//****************************************************************************80

/*	void dspdi ( double ap[], int n, int kpvt[], double det[2], int inert[3], 
			double work[], int job )
	{
		double ak;
		double akkp1;
		double akp1;
		double d;
		boolean dodet;
		boolean doert;
		boolean doinv;
		int ij;
		int ik;
		int ikp1;
		int iks;
		int j;
		int jb;
		int jk;
		int jkp1;
		int k;
		int kk;
		int kkp1;
		int km1;
		int ks;
		int ksj;
		int kskp1;
		int kstep;
		double t;
		double temp;

		doinv = ( job %   10 )       != 0;
		dodet = ( job %  100 ) /  10 != 0;
		doert = ( job % 1000 ) / 100 != 0;

		if ( dodet || doert )
		{
			if ( doert )
			{
				inert[0] = 0;
				inert[1] = 0;
				inert[2] = 0;
			}

			if ( dodet )
			{
				det[0] = 1.0;
				det[1] = 0.0;
			}

			t = 0.0;
			ik = 0;

			for ( k = 1; k <= n; k++ )
			{
				kk = ik + k;
				d = ap[kk-1];
				//
				//  2 by 2 block
				//  use det (d  s)  =  (d/t * c - t) * t,  t = abs ( s )
				//	          (s  c)
				//  to avoid underflow/overflow troubles.
				//
				//  Take two passes through scaling.  Use T for flag.
				//
				if ( kpvt[k-1] <= 0 )
				{
					if ( t == 0.0 )
					{
						ikp1 = ik + k;
						kkp1 = ikp1 + k;
						t = Blas1.r8_abs ( ap[kkp1-1] );
						d = ( d / t ) * ap[kkp1] - t;
					}
					else
					{
						d = t;
						t = 0.0;
					}
				}

				if ( doert )
				{
					if ( 0.0 < d )
					{
						inert[0] = inert[0] + 1;
					}
					else if ( d < 0.0 )
					{
						inert[1] = inert[1] + 1;
					}
					else if ( d == 0.0 )
					{
						inert[2] = inert[2] + 1;
					}
				}

				if ( dodet )
				{
					det[0] = det[0] * d;

					if ( det[0] != 0.0 )
					{
						while ( Blas1.r8_abs ( det[0] ) < 1.0 )
						{
							det[0] = det[0] * 10.0;
							det[1] = det[1] - 1.0;
						}

						while ( 10.0 <= Blas1.r8_abs ( det[0] ) )
						{
							det[0] = det[0] / 10.0;
							det[1] = det[1] + 1.0;
						}
					}
				}
				ik = ik + k;
			}
		}
		//
		//  Compute inverse(A).
		//
		if ( doinv )
		{
			k = 1;
			ik = 0;

			while ( k <= n )
			{
				km1 = k - 1;
				kk = ik + k;
				ikp1 = ik + k;
				kkp1 = ikp1 + k;

				if ( 0 <= kpvt[k-1] )
				{
					//
					//  1 by 1.
					//
					ap[kk-1] = 1.0 / ap[kk-1];

					if ( 2 <= k )
					{
						Blas1.dcopy ( k-1, ap, ik, 1, work, 1 );
						ij = 0;

						for ( j = 1; j <= k-1; j++ )
						{
							jk = ik + j;
							ap[jk-1] = Blas1.ddot ( j, ap, ij, 1, work, 1 );
							Blas1.daxpy ( j-1, work[j-1], ap, ij, 1, ap, ik, 1 );
							ij = ij + j;
						}
						ap[kk-1] = ap[kk-1] + Blas1.ddot ( k-1, work, 1, ap, ik, 1 );
					}
					kstep = 1;
				}
				else
				{
					//
					//  2 by 2.
					//
					t = Blas1.r8_abs ( ap[kkp1-1] );
					ak = ap[kk-1] / t;
					akp1 = ap[kkp1] / t;
					akkp1 = ap[kkp1-1] / t;
					d = t * ( ak * akp1 - 1.0 );
					ap[kk-1] = akp1 / d;
					ap[kkp1] = ak / d;
					ap[kkp1-1] = -akkp1 / d;

					if ( 1 <= km1 )
					{
						Blas1.dcopy ( km1, ap, ikp1, 1, work, 1 );
						ij = 0;

						for ( j = 1; j <= km1; j++ )
						{
							jkp1 = ikp1 + j;
							ap[jkp1-1] = Blas1.ddot ( j, ap, ij, 1, work, 1 );
							Blas1.daxpy ( j-1, work[j-1], ap, ij, 1, ap, ikp1, 1 );
							ij = ij + j;
						}

						ap[kkp1] = ap[kkp1] + Blas1.ddot ( km1, work, 1, ap, ikp1, 1 );
						ap[kkp1-1] = ap[kkp1-1] + Blas1.ddot ( km1, ap, ik, 1, ap, ikp1, 1 );
						Blas1.dcopy ( km1, ap, ik, 1, work, 1 );
						ij = 0;

						for ( j = 1; j <= km1; j++ )
						{
							jk = ik + j;
							ap[jk-1] = Blas1.ddot ( j, ap, ij, 1, work, 1 );
							Blas1.daxpy ( j-1, work[j-1], ap, ij, 1, ap, ik, 1 );
							ij = ij + j;
						}
						ap[kk-1] = ap[kk-1] + Blas1.ddot ( km1, work, 1, ap, ik, 1 );
					}
					kstep = 2;
				}
				//
				//  Swap.
				//
				ks = abs ( kpvt[k-1] );

				if ( ks != k )
				{
					iks = ( ks * ( ks - 1 ) ) / 2;
					Blas1.dswap  ( ks, ap, iks, 1, ap, ik, 1 );
					ksj = ik + ks;

					for ( jb = ks; jb <= k; jb++ )
					{
						j = k + ks - jb;
						jk = ik + j;
						temp = ap[jk-1];
						ap[jk-1] = ap[ksj-1];
						ap[ksj-1] = temp;
						ksj = ksj - ( j - 1 );
					}

					if ( kstep != 1 )
					{
						kskp1 = ikp1 + ks;
						temp = ap[kskp1-1];
						ap[kskp1-1] = ap[kkp1-1];
						ap[kkp1-1] = temp;
					}
				}
				ik = ik + k;
				if ( kstep == 2 )
				{
					ik = ik + k + 1;
				}
				k = k + kstep;
			}
		}

		return;
	}*/
	//****************************************************************************80

	int dspfa ( double ap[], int n, int kpvt[] )
	{
		double absakk;
		double ak;
		double akm1;
		double alpha;
		double bk;
		double bkm1;
		double colmax;
		double denom;
		int ij;
		int ijj;
		int ik;
		int ikm1;
		int im=0;
		int imax;
		int imaxp1;
		int imim;
		int imj;
		int imk;
		int info;
		int j;
		int jj;
		int jk;
		int jkm1;
		int jmax;
		int jmim;
		int k;
		int kk;
		int km1;
		int km1k;
		int km1km1;
		int kstep;
		double mulk;
		double mulkm1;
		double rowmax;
		boolean swap;
		double t;
		//
		//  ALPHA is used in choosing pivot block size.
		//
		alpha = ( 1.0 + Math.sqrt ( 17.0 ) ) / 8.0;

		info = 0;
		//
		//  Main loop on K, which goes from N to 1.
		//
		k = n;
		ik = ( n * ( n - 1 ) ) / 2;

		for ( ; ; )
		{
			//
			//  Leave the loop if K = 0 or K = 1.
			//
			if ( k == 0 )
			{
				break;
			}

			if ( k == 1 )
			{
				kpvt[0] = 1;
				if ( ap[0] == 0.0 )
				{
					info = 1;
				}
				break;
			}
			//
			//  This section of code determines the kind of elimination to be performed.
			//  When it is completed, KSTEP will be set to the size of the pivot block,
			//  and SWAP will be set to .true. if an interchange is required.
			//
			km1 = k - 1;
			kk = ik + k;
			absakk = Blas1.r8_abs ( ap[kk-1] );
			//
			//  Determine the largest off-diagonal element in column K.
			//
			imax = Blas1.idamax  ( k-1, ap, ik, 1 );
			imk = ik + imax;
			colmax = Blas1.r8_abs ( ap[imk-1] );

			if ( alpha * colmax <= absakk )
			{
				kstep = 1;
				swap = false;
			}
			//
			//  Determine the largest off-diagonal element in row IMAX.
			//
			else
			{
				rowmax = 0.0;
				imaxp1 = imax + 1;
				im = ( imax * ( imax - 1 ) ) / 2;
				imj = im + 2 * imax;

				for ( j = imaxp1; j <= k; j++ )
				{
					rowmax = Blas1.r8_max  ( rowmax, Blas1.r8_abs ( ap[imj-1] ) );
					imj = imj + j;
				}

				if ( imax != 1 )
				{
					jmax = Blas1.idamax  ( imax-1, ap, im, 1 );
					jmim = jmax + im;
					rowmax = Blas1.r8_max  ( rowmax, Blas1.r8_abs ( ap[jmim-1] ) );
				}

				imim = imax + im;

				if ( alpha * rowmax <= Blas1.r8_abs ( ap[imim-1] ) )
				{
					kstep = 1;
					swap = true;
				}
				else if ( alpha * colmax * ( colmax / rowmax ) <= absakk )
				{
					kstep = 1;
					swap = false;
				}
				else
				{
					kstep = 2;
					swap = imax != km1;
				}
			}
			//
			//  Column K is zero.  Set INFO and iterate the loop.
			//
			if ( Blas1.r8_max  ( absakk, colmax ) == 0.0 )
			{
				kpvt[k-1] = k;
				info = k;
			}
			else
			{
				if ( kstep != 2 )
				{
					//
					//  1 x 1 pivot block.
					//
					if ( swap )
					{
						//
						//  Perform an interchange.
						//
						Blas1.dswap  ( imax, ap, im, 1, ap, ik, 1 );
						imj = ik + imax;

						for ( jj = imax; jj <= k; jj++ )
						{
							j = k + imax - jj;
							jk = ik + j;
							t = ap[jk-1];
							ap[jk-1] = ap[imj-1];
							ap[imj-1] = t;
							imj = imj - ( j - 1 );
						}
					}
					//
					//  Perform the elimination.
					//
					ij = ik - ( k - 1 );

					for ( jj = 1; jj <= km1; jj++ )
					{
						j = k - jj;
						jk = ik + j;
						mulk = -ap[jk-1] / ap[kk-1];
						t = mulk;
						Blas1.daxpy ( j, t, ap, ik, 1, ap, ij, 1 );
						ijj = ij + j;
						ap[jk-1] = mulk;
						ij = ij - ( j - 1 );
					}
					//
					//  Set the pivot array.
					//
					if ( swap )
					{
						kpvt[k-1] = imax;
					}
					else
					{
						kpvt[k-1] = k;
					}
				}
				else
				{
					//
					//  2 x 2 pivot block.
					//
					km1k = ik + k - 1;
					ikm1 = ik - ( k - 1 );
					//
					//  Perform an interchange.
					//
					if ( swap )
					{
						Blas1.dswap  ( imax, ap, im, 1, ap, ikm1, 1 );
						imj = ikm1 + imax;

						for ( jj = imax; jj <= km1; jj++ )
						{
							j = km1 + imax - jj;
							jkm1 = ikm1 + j;
							t = ap[jkm1-1];
							ap[jkm1-1] = ap[imj-1];
							ap[imj-1] = t;
							imj = imj - ( j - 1 );
						}
						t = ap[km1k-1];
						ap[km1k-1] = ap[imk-1];
						ap[imk-1] = t;
					}
					//
					//  Perform the elimination.
					//
					if ( k-2 != 0 )
					{
						ak = ap[kk-1] / ap[km1k-1];
						km1km1 = ikm1 + k - 1;
						akm1 = ap[km1km1-1] / ap[km1k-1];
						denom = 1.0 - ak * akm1;
						ij = ik - ( k - 1 ) - ( k - 2 );

						for ( jj = 1; jj <= k-2; jj++ )
						{
							j = km1 - jj;
							jk = ik + j;
							bk = ap[jk-1] / ap[km1k-1];
							jkm1 = ikm1 + j;
							bkm1 = ap[jkm1-1] / ap[km1k-1];
							mulk = ( akm1 * bk - bkm1 ) / denom;
							mulkm1 = ( ak * bkm1 - bk ) / denom;
							t = mulk;
							Blas1.daxpy ( j, t, ap, ik, 1, ap, ij, 1 );
							t = mulkm1;
							Blas1.daxpy ( j, t, ap, ikm1, 1, ap, ij, 1 );
							ap[jk-1] = mulk;
							ap[jkm1-1] = mulkm1;
							ijj = ij + j;
							ij = ij - ( j - 1 );
						}
					}
					//
					//  Set the pivot array.
					//
					if ( swap )
					{
						kpvt[k-1] = -imax;
					}
					else
					{
						kpvt[k-1] = 1 - k;
					}
					kpvt[k-2] = kpvt[k-1];
				}
			}

			ik = ik - ( k - 1 );
			if ( kstep == 2 )
			{
				ik = ik - ( k - 2 );
			}

			k = k - kstep;

		}

		return info;
	}
	//****************************************************************************80

	void dspsl ( double ap[], int n, int kpvt[], double b[] )
	{
		double ak;
		double akm1;
		double bk;
		double bkm1;
		double denom;
		int ik;
		int ikm1;
		int ikp1;
		int k;
		int kk;
		int km1k;
		int km1km1;
		int kp;
		double temp;
		//
		//  Loop backward applying the transformations and D inverse to B.
		//
		k = n;
		ik = ( n * ( n - 1 ) ) / 2;

		while ( 0 < k )
		{
			kk = ik + k;

			if ( 0 <= kpvt[k-1] )
			{
				//
				//  1 x 1 pivot block.
				//
				if ( k != 1 )
				{
					kp = kpvt[k-1];
					//
					//  Interchange.
					//
					if ( kp != k )
					{
						temp = b[k-1];
						b[k-1] = b[kp-1];
						b[kp-1] = temp;
					}
					//
					//  Apply the transformation.
					//
					Blas1.daxpy ( k-1, b[k-1], ap, ik, 1, b, 0, 1 );
				}
				//
				//  Apply D inverse.
				//
				b[k-1] = b[k-1] / ap[kk-1];
				k = k - 1;
				ik = ik - k;
			}
			else
			{
				//
				//  2 x 2 pivot block.
				//
				ikm1 = ik - ( k - 1 );

				if ( k != 2 )
				{
					kp = Math.abs ( kpvt[k-1] );
					//
					//  Interchange.
					//
					if ( kp != k-1 )
					{
						temp = b[k-2];
						b[k-2] = b[kp-1];
						b[kp-1] = temp;
					}
					//
					//  Apply the transformation.
					//
					Blas1.daxpy ( k-2, b[k-1], ap, ik, 1, b, 0, 1 );
					Blas1.daxpy ( k-2, b[k-2], ap, ikm1, 1, b, 0, 1 );
				}
				//
				//  Apply D inverse.
				//
				km1k = ik + k - 1;
				kk = ik + k;
				ak = ap[kk-1] / ap[km1k-1];
				km1km1 = ikm1 + k - 1;
				akm1 = ap[km1km1-1] / ap[km1k-1];
				bk = b[k-1] / ap[km1k-1];
				bkm1 = b[k-2] / ap[km1k-1];
				denom = ak * akm1 - 1.0;
				b[k-1] = ( akm1 * bk - bkm1 ) / denom;
				b[k-2] = ( ak * bkm1 - bk ) / denom;
				k = k - 2;
				ik = ik - ( k + 1 ) - k;
			}
		}
		//
		//  Loop forward applying the transformations.
		//
		k = 1;
		ik = 0;

		while ( k <= n )
		{
			if ( 0 <= kpvt[k-1] )
			{
				//
				//  1 x 1 pivot block.
				//
				if ( k != 1 )
				{
					//
					//  Apply the transformation.
					//
					b[k-1] = b[k-1] + Blas1.ddot ( k-1, ap, ik, 1, b, 0, 1 );
					kp = kpvt[k-1];
					//
					//  Interchange.
					//
					if ( kp != k )
					{
						temp = b[k-1];
						b[k-1] = b[kp-1];
						b[kp-1] = temp;
					}
				}
				ik = ik + k;
				k = k + 1;
			}
			else
			{
				//
				//  2 x 2 pivot block.
				//
				if ( k != 1 )
				{
					//
					//  Apply the transformation.
					//
					b[k-1] = b[k-1] + Blas1.ddot ( k-1, ap, ik, 1, b, 0, 1 );
					ikp1 = ik + k;
					b[k] = b[k] + Blas1.ddot ( k-1, ap, ikp1, 1, b, 0, 1 );
					kp = Math.abs ( kpvt[k-1] );
					//
					//  Interchange.
					//
					if ( kp != k )
					{
						temp = b[k-1];
						b[k-1] = b[kp-1];
						b[kp-1] = temp;
					}
				}
				ik = ik + k + k + 1;
				k = k + 2;
			}
		}

		return;
	}
	//****************************************************************************80

	static int dsvdc ( double a[], int lda, int m, int n, double s[], double e[], 
			double u[], int ldu, double v[], int ldv, double work[], int job )
	{
		double b;
		double c;
		double cs=0.0;
		double el;
		double emm1;
		double f;
		double g;
		int i;
		int info;
		int iter;
		int j;
		int jobu;
		int k;
		int kase;
		int kk;
		int l;
		int ll;
		int lls;
		int ls=0;
		int lu;
		int maxit = 30;
		int mm;
		int mm1;
		int mn;
		int mp1;
		int nct;
		int nctp1;
		int ncu;
		int nrt;
		int nrtp1;
		double scale;
		double shift;
		double sl;
		double sm;
		double smm1;
		double sn=0.0;
		double t;
		double t1;
		double test;
		boolean wantu;
		boolean wantv;
		double ztest;
		//
		//  Determine what is to be computed.
		//
		info = 0;
		wantu = false;
		wantv = false;
		jobu = ( job % 100 ) / 10;

		if ( 1 < jobu )
		{
			ncu = Blas1.i4_min ( m, n );
		}
		else
		{
			ncu = m;
		}

		if ( jobu != 0 )
		{
			wantu = true;
		}

		if ( ( job % 10 ) != 0 )
		{
			wantv = true;
		}
		//
		//  Reduce A to bidiagonal form, storing the diagonal elements
		//  in S and the super-diagonal elements in E.
		//
		nct = Blas1.i4_min ( m-1, n );
		nrt = Blas1.i4_max ( 0, Blas1.i4_min ( m, n-2 ) );
		lu = Blas1.i4_max ( nct, nrt );

		for ( l = 1; l <= lu; l++ )
		{
			//
			//  Compute the transformation for the L-th column and
			//  place the L-th diagonal in S(L).
			//
			if ( l <= nct )
			{
				s[l-1] = Blas1.dnrm2 ( m-l+1, a, l-1+(l-1)*lda, 1 );

				if ( s[l-1] != 0.0 )
				{
					if ( a[l-1+(l-1)*lda] != 0.0 )
					{
						s[l-1] = Blas1.r8_sign ( a[l-1+(l-1)*lda] ) * Blas1.r8_abs ( s[l-1] );
					}
					Blas1.dscal  ( m-l+1, 1.0 / s[l-1], a, l-1+(l-1)*lda, 1 );
					a[l-1+(l-1)*lda] = 1.0 + a[l-1+(l-1)*lda];
				}
				s[l-1] = -s[l-1];
			}

			for ( j = l+1; j <= n; j++ )
			{
				//
				//  Apply the transformation.
				//
				if ( l <= nct && s[l-1] != 0.0 )
				{
					t = - Blas1.ddot ( m-l+1, a, l-1+(l-1)*lda, 1, a, l-1+(j-1)*lda, 1 ) 
							/ a[l-1+(l-1)*lda];
					Blas1.daxpy ( m-l+1, t, a, l-1+(l-1)*lda, 1, a, l-1+(j-1)*lda, 1 );
				}
				//
				//  Place the L-th row of A into E for the
				//  subsequent calculation of the row transformation.
				//
				e[j-1] = a[l-1+(j-1)*lda];
			}
			//
			//  Place the transformation in U for subsequent back multiplication.
			//
			if ( wantu && l <= nct )
			{
				for ( i = l; i <= m; i++ )
				{
					u[i-1+(l-1)*ldu] = a[i-1+(l-1)*lda];
				}
			}

			if ( l <= nrt )
			{
				//
				//  Compute the L-th row transformation and place the
				//  L-th superdiagonal in E(L).
				//
				e[l-1] = Blas1.dnrm2 ( n-l, e, l, 1 );

				if ( e[l-1] != 0.0 )
				{
					if ( e[l] != 0.0 )
					{
						e[l-1] = Blas1.r8_sign ( e[l] ) * Blas1.r8_abs ( e[l-1] );
					}
					Blas1.dscal  ( n-l, 1.0 / e[l-1], e, l, 1 );
					e[l] = 1.0 + e[l];
				}

				e[l-1] = -e[l-1];
				//
				//  Apply the transformation.
				//
				if ( l+1 <= m && e[l-1] != 0.0 )
				{
					for ( j = l+1; j <= m; j++ )
					{
						work[j-1] = 0.0;
					}

					for ( j = l+1; j <= n; j++ )
					{
						Blas1.daxpy ( m-l, e[j-1], a, l+(j-1)*lda, 1, work, l, 1 );
					}

					for ( j = l+1; j <= n; j++ )
					{
						Blas1.daxpy ( m-l, -e[j-1]/e[l], work, l, 1, a, l+(j-1)*lda, 1 );
					}
				}
				//
				//  Place the transformation in V for subsequent back multiplication.
				//
				if ( wantv )
				{
					for ( j = l+1; j <= n; j++ )
					{
						v[j-1+(l-1)*ldv] = e[j-1];
					}
				}
			}
		}
		//
		//  Set up the final bidiagonal matrix of order MN.
		//
		mn = Blas1.i4_min ( m + 1, n );
		nctp1 = nct + 1;
		nrtp1 = nrt + 1;

		if ( nct < n )
		{
			s[nctp1-1] = a[nctp1-1+(nctp1-1)*lda];
		}

		if ( m < mn )
		{
			s[mn-1] = 0.0;
		}

		if ( nrtp1 < mn )
		{
			e[nrtp1-1] = a[nrtp1-1+(mn-1)*lda];
		}

		e[mn-1] = 0.0;
		//
		//  If required, generate U.
		//
		if ( wantu )
		{
			for ( i = 1; i <= m; i++ )
			{
				for ( j = nctp1; j <= ncu; j++ )
				{
					u[(i-1)+(j-1)*ldu] = 0.0;
				}
			}

			for ( j = nctp1; j <= ncu; j++ )
			{
				u[j-1+(j-1)*ldu] = 1.0;
			}

			for ( ll = 1; ll <= nct; ll++ )
			{
				l = nct - ll + 1;

				if ( s[l-1] != 0.0 )
				{
					for ( j = l+1; j <= ncu; j++ )
					{
						t = - Blas1.ddot ( m-l+1, u, (l-1)+(l-1)*ldu, 1, u, (l-1)+(j-1)*ldu, 1 ) 
								/ u[l-1+(l-1)*ldu];
						Blas1.daxpy ( m-l+1, t, u, (l-1)+(l-1)*ldu, 1, u, (l-1)+(j-1)*ldu, 1 );
					}

					Blas1.dscal  ( m-l+1, -1.0, u, (l-1)+(l-1)*ldu, 1 );
					u[l-1+(l-1)*ldu] = 1.0 + u[l-1+(l-1)*ldu];
					for ( i = 1; i <= l-1; i++ )
					{
						u[i-1+(l-1)*ldu] = 0.0;
					}
				}
				else
				{
					for ( i = 1; i <= m; i++ )
					{
						u[i-1+(l-1)*ldu] = 0.0;
					}
					u[l-1+(l-1)*ldu] = 1.0;
				}
			}
		}
		//
		//  If it is required, generate V.
		//
		if ( wantv )
		{
			for ( ll = 1; ll <= n; ll++ )
			{
				l = n - ll + 1;

				if ( l <= nrt && e[l-1] != 0.0 )
				{
					for ( j = l+1; j <= n; j++ )
					{
						t = - Blas1.ddot ( n-l, v, l+(l-1)*ldv, 1, v, l+(j-1)*ldv, 1 ) 
								/ v[l+(l-1)*ldv];
						Blas1.daxpy ( n-l, t, v, l+(l-1)*ldv, 1, v, l+(j-1)*ldv, 1 );
					}

				}
				for ( i = 1; i <= n; i++ )
				{
					v[i-1+(l-1)*ldv] = 0.0;
				}
				v[l-1+(l-1)*ldv] = 1.0;
			}
		}
		//
		//  Main iteration loop for the singular values.
		//
		mm = mn;
		iter = 0;

		while ( 0 < mn )
		{
			//
			//  If too many iterations have been performed, set flag and return.
			//
			if ( maxit <= iter )
			{
				info = mn;
				return info;
			}
			//
			//  This section of the program inspects for
			//  negligible elements in the S and E arrays.
			//
			//  On completion the variables KASE and L are set as follows:
			//
			//  KASE = 1     if S(MN) and E(L-1) are negligible and L < MN
			//  KASE = 2     if S(L) is negligible and L < MN
			//  KASE = 3     if E(L-1) is negligible, L < MN, and
			//	               S(L), ..., S(MN) are not negligible (QR step).
			//  KASE = 4     if E(MN-1) is negligible (convergence).
			//
			for ( ll = 1; ll <= mn; ll++ )
			{
				l = mn - ll;

				if ( l == 0 )
				{
					break;
				}

				test = Blas1.r8_abs ( s[l-1] ) + Blas1.r8_abs ( s[l] );
				ztest = test + Blas1.r8_abs ( e[l-1] );

				if ( ztest == test )
				{
					e[l-1] = 0.0;
					break;
				}
			}

			if ( l == mn - 1 )
			{
				kase = 4;
			}
			else
			{
				mp1 = mn + 1;

				for ( lls = l+1; lls <= mn+1; lls++ )
				{
					ls = mn - lls + l + 1;

					if ( ls == l )
					{
						break;
					}

					test = 0.0;
					if ( ls != mn )
					{
						test = test + Blas1.r8_abs ( e[ls-1] );
					}

					if ( ls != l + 1 )
					{
						test = test + Blas1.r8_abs ( e[ls-2] );
					}

					ztest = test + Blas1.r8_abs ( s[ls-1] );

					if ( ztest == test )
					{
						s[ls-1] = 0.0;
						break;
					}

				}

				if ( ls == l )
				{
					kase = 3;
				}
				else if ( ls == mn )
				{
					kase = 1;
				}
				else
				{
					kase = 2;
					l = ls;
				}
			}

			l = l + 1;
			//
			//  Deflate negligible S(MN).
			//
			if ( kase == 1 )
			{
				mm1 = mn - 1;
				f = e[mn-2];
				e[mn-2] = 0.0;

				for ( kk = 1; kk <= mm1; kk++ )
				{
					k = mm1 - kk + l;
					t1 = s[k-1];
					//Blas1.drotg ( &t1, &f, &cs, &sn );
					//////////////////////////////////
					double [] drotg = new double [4];
					drotg[0] = t1;
					drotg[1] = f;
					drotg[2] = cs;
					drotg[3] = sn;
					///////////////////////////////////
					Blas1.drotg (drotg);
					///////////////////////////////////
					t1 = drotg[0];
					f = drotg[1];
					cs = drotg[2];
					sn = drotg[3];
					///////////////////////////////////
					s[k-1] = t1;

					if ( k != l )
					{
						f = -sn * e[k-2];
						e[k-2] = cs * e[k-2];
					}

					if ( wantv )
					{
						Blas1.drot ( n, v, 0+(k-1)*ldv, 1, v, 0+(mn-1)*ldv, 1, cs, sn );
					}
				}
			}
			//
			//  Split at negligible S(L).
			//
			else if ( kase == 2 )
			{
				f = e[l-2];
				e[l-2] = 0.0;

				for ( k = l; k <= mn; k++ )
				{
					t1 = s[k-1];
					//Blas1.drotg ( &t1, &f, &cs, &sn );
					//////////////////////////////////
					double [] drotg = new double [4];
					drotg[0] = t1;
					drotg[1] = f;
					drotg[2] = cs;
					drotg[3] = sn;
					///////////////////////////////////
					Blas1.drotg (drotg);
					///////////////////////////////////
					t1 = drotg[0];
					f = drotg[1];
					cs = drotg[2];
					sn = drotg[3];
					///////////////////////////////////
					s[k-1] = t1;
					f = - sn * e[k-1];
					e[k-1] = cs * e[k-1];
					if ( wantu )
					{
						Blas1.drot ( m, u, 0+(k-1)*ldu, 1, u, 0+(l-2)*ldu, 1, cs, sn );
					}
				}
			}
			//
			//  Perform one QR step.
			//
			else if ( kase == 3 )
			{
				//
				//  Calculate the shift.
				//
				scale = Blas1.r8_max  ( Blas1.r8_abs ( s[mn-1] ), 
						Blas1.r8_max  ( Blas1.r8_abs ( s[mn-2] ), 
								Blas1.r8_max  ( Blas1.r8_abs ( e[mn-2] ), 
										Blas1.r8_max  ( Blas1.r8_abs ( s[l-1] ), Blas1.r8_abs ( e[l-1] ) ) ) ) );

				sm = s[mn-1] / scale;
				smm1 = s[mn-2] / scale;
				emm1 = e[mn-2] / scale;
				sl = s[l-1] / scale;
				el = e[l-1] / scale;
				b = ( ( smm1 + sm ) * ( smm1 - sm ) + emm1 * emm1 ) / 2.0;
				c = ( sm * emm1 ) * ( sm * emm1 );
				shift = 0.0;

				if ( b != 0.0 || c != 0.0 )
				{
					shift = Math.sqrt ( b * b + c );
					if ( b < 0.0 )
					{
						shift = -shift;
					}
					shift = c / ( b + shift );
				}

				f = ( sl + sm ) * ( sl - sm ) - shift;
				g = sl * el;
				//
				//  Chase zeros.
				//
				mm1 = mn - 1;

				for ( k = l; k <= mm1; k++ )
				{
					//Blas1.drotg ( &f, &g, &cs, &sn );
					//////////////////////////////////
					double [] drotg = new double [4];
					drotg[0] = f;
					drotg[1] = g;
					drotg[2] = cs;
					drotg[3] = sn;
					///////////////////////////////////
					Blas1.drotg (drotg);
					///////////////////////////////////
					f = drotg[0];
					g = drotg[1];
					cs = drotg[2];
					sn = drotg[3];
					///////////////////////////////////
					if ( k != l )
					{
						e[k-2] = f;
					}

					f = cs * s[k-1] + sn * e[k-1];
					e[k-1] = cs * e[k-1] - sn * s[k-1];
					g = sn * s[k];
					s[k] = cs * s[k];

					if ( wantv )
					{
						Blas1.drot ( n, v, 0+(k-1)*ldv, 1, v, 0+k*ldv, 1, cs, sn );
					}

					//Blas1.drotg ( &f, &g, &cs, &sn );
					//////////////////////////////////
					drotg[0] = f;
					drotg[1] = g;
					drotg[2] = cs;
					drotg[3] = sn;
					///////////////////////////////////
					Blas1.drotg (drotg);
					///////////////////////////////////
					f = drotg[0];
					g = drotg[1];
					cs = drotg[2];
					sn = drotg[3];
					///////////////////////////////////
					s[k-1] = f;
					f = cs * e[k-1] + sn * s[k];
					s[k] = -sn * e[k-1] + cs * s[k];
					g = sn * e[k];
					e[k] = cs * e[k];

					if ( wantu && k < m )
					{
						Blas1.drot ( m, u, 0+(k-1)*ldu, 1, u, 0+k*ldu, 1, cs, sn );
					}
				}
				e[mn-2] = f;
				iter = iter + 1;
			}
			//
			//  Convergence.
			//
			else if ( kase == 4 )
			{
				//
				//  Make the singular value nonnegative.
				//
				if ( s[l-1] < 0.0 )
				{
					s[l-1] = -s[l-1];
					if ( wantv )
					{
						Blas1.dscal  ( n, -1.0, v, 0+(l-1)*ldv, 1 );
					}
				}
				//
				//  Order the singular value.
				//
				for ( ; ; )
				{
					if ( l == mm )
					{
						break;
					}

					if ( s[l] <= s[l-1] )
					{
						break;
					}

					t = s[l-1];
					s[l-1] = s[l];
					s[l] = t;

					if ( wantv && l < n )
					{
						Blas1.dswap  ( n, v, 0+(l-1)*ldv, 1, v, 0+l*ldv, 1 );
					}

					if ( wantu && l < m )
					{
						Blas1.dswap  ( m, u, 0+(l-1)*ldu, 1, u, 0+l*ldu, 1 );
					}

					l = l + 1;
				}
				iter = 0;
				mn = mn - 1;
			}
		}

		return info;
	}
	//****************************************************************************80

	double dtrco ( double t[], int ldt, int n, double z[], int job )
	{
		double ek;
		int i;
		int i1;
		int j;
		int j1;
		int j2;
		int k;
		int kk;
		int l;
		boolean lower;
		double rcond;
		double s;
		double sm;
		double temp;
		double tnorm;
		double w;
		double wk;
		double wkm;
		double ynorm;

		lower = ( job == 0 );
		//
		//  Compute the 1-norm of T.
		//
		tnorm = 0.0;

		for ( j = 1; j <= n; j++ )
		{
			if ( lower )
			{
				l = n + 1 - j;
				i1 = j;
			}
			else
			{
				l = j;
				i1 = 1;
			}
			tnorm = Blas1.r8_max  ( tnorm, Blas1.dasum ( l, t, i1-1+(j-1)*ldt, 1 ) );
		}
		//
		//  RCOND = 1/(norm(T)*(estimate of norm(inverse(T)))).
		//
		//  Estimate = norm(Z)/norm(Y) where T * Z = Y and T' * Y = E.
		//
		//  T' is the transpose of T.
		//
		//  The components of E are chosen to cause maximum local
		//  growth in the elements of Y.
		//
		//  The vectors are frequently rescaled to avoid overflow.
		//
		//  Solve T' * Y = E.
		//
		ek = 1.0;
		for ( i = 1; i <= n; i++ )
		{
			z[i-1] = 0.0;
		}

		for ( kk = 1; kk <= n; kk++ )
		{
			if ( lower )
			{
				k = n + 1 - kk;
			}
			else
			{
				k = kk;
			}

			if ( z[k-1] != 0.0 )
			{
				ek = Blas1.r8_sign ( -z[k-1] ) * ek;
			}

			if ( Blas1.r8_abs ( t[k-1+(k-1)*ldt] ) < Blas1.r8_abs ( ek - z[k-1] ) )
			{
				s = Blas1.r8_abs ( t[k-1+(k-1)*ldt] ) / Blas1.r8_abs ( ek - z[k-1] );
				for ( i = 1; i <= n; i++ )
				{
					z[i-1] = s * z[i-1];
				}
				ek = s * ek;
			}

			wk = ek - z[k-1];
			wkm = -ek - z[k-1];
			s = Blas1.r8_abs ( wk );
			sm = Blas1.r8_abs ( wkm );

			if ( t[k-1+(k-1)*ldt] != 0.0 )
			{
				wk = wk / t[k-1+(k-1)*ldt];
				wkm = wkm / t[k-1+(k-1)*ldt];
			}
			else
			{
				wk = 1.0;
				wkm = 1.0;
			}

			if ( kk != n )
			{
				if ( lower )
				{
					j1 = 1;
					j2 = k - 1;
				}
				else
				{
					j1 = k + 1;
					j2 = n;
				}

				for ( j = j1; j <= j2; j++ )
				{
					sm = sm + Blas1.r8_abs ( z[j-1] + wkm * t[k-1+(j-1)*ldt] );
					z[j-1] = z[j-1] + wk * t[k-1+(j-1)*ldt];
					s = s + Blas1.r8_abs ( z[j-1] );
				}

				if ( s < sm )
				{
					w = wkm - wk;
					wk = wkm;
					for ( j = j1; j <= j2; j++ )
					{
						z[j-1] = z[j-1] + w * t[k-1+(j-1)*ldt];
					}
				}
			}
			z[k-1] = wk;
		}

		temp = Blas1.dasum ( n, z, 0, 1 );

		for ( i = 1; i <= n; i++ )
		{
			z[i-1] = z[i-1] / temp;
		}

		ynorm = 1.0;
		//
		//  Solve T * Z = Y.
		//
		for ( kk = 1; kk <= n; kk++ )
		{
			if ( lower )
			{
				k = kk;
			}
			else
			{
				k = n + 1 - kk;
			}

			if ( Blas1.r8_abs ( t[k-1+(k-1)*ldt] ) < Blas1.r8_abs ( z[k-1] ) )
			{
				s = Blas1.r8_abs ( t[k-1+(k-1)*ldt] ) / Blas1.r8_abs ( z[k-1] );
				for ( i = 1; i <= n; i++ )
				{
					z[i-1] = s * z[i-1];
				}
				ynorm = s * ynorm;
			}

			if ( t[k-1+(k-1)*ldt] != 0.0 )
			{
				z[k-1] = z[k-1] / t[k-1+(k-1)*ldt];
			}
			else
			{
				z[k-1] = 1.0;
			}

			if ( lower )
			{
				i1 = k + 1;
			}
			else
			{
				i1 = 1;
			}

			if ( kk < n )
			{
				w = -z[k-1];
				Blas1.daxpy ( n-kk, w, t, i1-1+(k-1)*ldt, 1, z, i1-1, 1 );
			}
		}
		//
		//  Make ZNORM = 1.0.
		//
		s = 1.0 / Blas1.dasum ( n, z, 0, 1 );
		for ( i = 1; i <= n; i++ )
		{
			z[i-1] = s * z[i-1];
		}
		ynorm = s * ynorm;

		if ( tnorm != 0.0 )
		{
			rcond = ynorm / tnorm;
		}
		else
		{
			rcond = 0.0;
		}

		return rcond;
	}
	//****************************************************************************80

	int dtrdi ( double t[], int ldt, int n, double det[], int job )
	{
		int i;
		int info;
		int j;
		int k;
		double temp;
		//
		//  Determinant.
		//
		info = 0;

		if ( job / 100 != 0 )
		{
			det[0] = 1.0;
			det[1] = 0.0;

			for ( i = 1; i <= n; i++ )
			{
				det[0] = det[0] * t[i-1+(i-1)*ldt];

				if ( det[0] == 0.0 )
				{
					break;
				}

				while ( Blas1.r8_abs ( det[0] ) < 1.0 )
				{
					det[0] = det[0] * 10.0;
					det[1] = det[1] - 1.0;
				}

				while ( 10.0 <= Blas1.r8_abs ( det[0] ) )
				{
					det[0] = det[0] / 10.0;
					det[1] = det[1] + 1.0;
				}
			}
		}

		if ( ( ( job / 10 ) % 10 ) == 0 )
		{
			return info;
		}
		//
		//  Inverse of an upper triangular matrix.
		//
		if ( ( job % 10 ) != 0 )
		{
			info = 0;

			for ( k = 1; k <= n; k++ )
			{
				if ( t[k-1+(k-1)*ldt] == 0.0 )
				{
					info = k;
					break;
				}

				t[k-1+(k-1)*ldt] = 1.0 / t[k-1+(k-1)*ldt];
				temp = -t[k-1+(k-1)*ldt];
				Blas1.dscal  ( k-1, temp, t, 0+(k-1)*ldt, 1 );

				for ( j = k + 1; j <= n; j++ )
				{
					temp = t[k-1+(j-1)*ldt];
					t[k-1+(j-1)*ldt] = 0.0;
					Blas1.daxpy ( k, temp, t, 0+(k-1)*ldt, 1, t, 0+(j-1)*ldt, 1 );
				}
			}
		}
		//
		//  Inverse of a lower triangular matrix.
		//
		else
		{
			info = 0;

			for ( k = n; 1 <= k; k-- )
			{
				if ( t[k-1+(k-1)*ldt] == 0.0 )
				{
					info = k;
					break;
				}

				t[k-1+(k-1)*ldt] = 1.0 / t[k-1+(k-1)*ldt];
				temp = -t[k-1+(k-1)*ldt];

				if ( k != n )
				{
					Blas1.dscal  ( n-k, temp, t, k+(k-1)*ldt, 1 );
				}

				for ( j = 1; j <= k-1; j++ )
				{
					temp = t[k-1+(j-1)*ldt];
					t[k-1+(j-1)*ldt] = 0.0;
					Blas1.daxpy ( n-k+1, temp, t, k-1+(k-1)*ldt, 1, t, k-1+(j-1)*ldt, 1 );
				}
			}
		}

		return info;
	}
	//****************************************************************************80

	int dtrsl ( double t[], int ldt, int n, double b[], int job )
	{
		int info;
		int j;
		int jj;
		int kase;
		double temp;
		//
		//  Check for zero diagonal elements.
		//
		for ( j = 1; j <= n; j++ )
		{
			if ( t[j-1+(j-1)*ldt] == 0.0 )
			{
				info = j;
				return info;
			}
		}

		info = 0;
		//
		//  Determine the task and go to it.
		//
		if ( ( job % 10 ) == 0 )
		{
			kase = 1;
		}
		else
		{
			kase = 2;
		}

		if ( ( job % 100 ) / 10 != 0 )
		{
			kase = kase + 2;
		}
		//
		//  Solve T * X = B for T lower triangular.
		//
		if ( kase == 1 )
		{
			b[0] = b[0] / t[0+0*ldt];
			for ( j = 2; j <= n; j++ )
			{
				temp = -b[j-2];
				Blas1.daxpy ( n-j+1, temp, t, (j-1)+(j-2)*ldt, 1, b, j-1, 1 );
				b[j-1] = b[j-1] / t[j-1+(j-1)*ldt];
			}
		}
		//
		//  Solve T * X = B for T upper triangular.
		//
		else if ( kase == 2 )
		{
			b[n-1] = b[n-1] / t[n-1+(n-1)*ldt];
			for ( jj = 2; jj <= n; jj++ )
			{
				j = n - jj + 1;
				temp = -b[j];
				Blas1.daxpy ( j, temp, t, 0+j*ldt, 1, b, 0, 1 );
				b[j-1] = b[j-1] / t[j-1+(j-1)*ldt];
			}
		}
		//
		//  Solve T' * X = B for T lower triangular.
		//
		else if ( kase == 3 )
		{
			b[n-1] = b[n-1] / t[n-1+(n-1)*ldt];
			for ( jj = 2; jj <= n; jj++ )
			{
				j = n - jj + 1;
				b[j-1] = b[j-1] - Blas1.ddot ( jj-1, t, j+(j-1)*ldt, 1, b, j, 1 );
				b[j-1] = b[j-1] / t[j-1+(j-1)*ldt];
			}
		}
		//
		//  Solve T' * X = B for T upper triangular.
		//
		else if ( kase == 4 )
		{
			b[0] = b[0] / t[0+0*ldt];
			for ( j = 2; j <= n; j++ )
			{
				b[j-1] = b[j-1] - Blas1.ddot ( j-1, t, 0+(j-1)*ldt, 1, b, 0, 1 );
				b[j-1] = b[j-1] / t[j-1+(j-1)*ldt];
			}
		}

		return info;
	}

}
