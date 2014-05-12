package affy.qualityControl;

public class Blas1 {
	//****************************************************************************80

	public static double dasum ( int n, double x[], int x_first, int incx )
	{
		int i;
		int j;
		double value;

		value = 0.0;
		j = 0;

		for ( i = 0; i < n; i++ )
		{
			value = value + r8_abs ( x[j+x_first] );
			j = j + incx;
		}

		return value;
	}
	//****************************************************************************80

	public static void daxpy ( int n, double da, double dx[], int dx_first, int incx, double dy[], int dy_first, int incy )
	{
		int i;
		int ix;
		int iy;
		int m;

		if ( n <= 0 )
		{
			return;
		}

		if ( da == 0.0 )
		{
			return;
		}
		//
		//  Code for unequal increments or equal increments
		//  not equal to 1.
		//
		if ( incx != 1 || incy != 1 )
		{
			if ( 0 <= incx )
			{
				ix = 0;
			}
			else
			{
				ix = ( - n + 1 ) * incx;
			}

			if ( 0 <= incy )
			{
				iy = 0;
			}
			else
			{
				iy = ( - n + 1 ) * incy;
			}

			for ( i = 0; i < n; i++ )
			{
				dy[iy+dy_first] = dy[iy+dy_first] + da * dx[ix+dx_first];
				ix = ix + incx;
				iy = iy + incy;
			}
		}
		//
		//  Code for both increments equal to 1.
		//
		else
		{
			m = n % 4;

			for ( i = 0; i < m; i++ )
			{
				dy[i+dy_first] = dy[i+dy_first] + da * dx[i+dx_first];
			}

			for ( i = m; i < n; i = i + 4 )
			{
				dy[i+dy_first  ] = dy[i+dy_first  ] + da * dx[i+dx_first  ];
				dy[i+1+dy_first] = dy[i+1+dy_first] + da * dx[i+1+dx_first];
				dy[i+2+dy_first] = dy[i+2+dy_first] + da * dx[i+2+dx_first];
				dy[i+3+dy_first] = dy[i+3+dy_first] + da * dx[i+3+dx_first];
			}

		}

		return;
	}
	//****************************************************************************80

	public static void dcopy ( int n, double dx[], int incx, double dy[], int incy )
	{
		int i;
		int ix;
		int iy;
		int m;

		if ( n <= 0 )
		{
			return;
		}

		if ( incx == 1 && incy == 1 )
		{
			m = n % 7;

			if ( m != 0 )
			{
				for ( i = 0; i < m; i++ )
				{
					dy[i] = dx[i];
				}
			}

			for ( i = m; i < n; i = i + 7 )
			{
				dy[i] = dx[i];
				dy[i + 1] = dx[i + 1];
				dy[i + 2] = dx[i + 2];
				dy[i + 3] = dx[i + 3];
				dy[i + 4] = dx[i + 4];
				dy[i + 5] = dx[i + 5];
				dy[i + 6] = dx[i + 6];
			}
		}
		else
		{
			if ( 0 <= incx )
			{
				ix = 0;
			}
			else
			{
				ix = ( -n + 1 ) * incx;
			}

			if ( 0 <= incy )
			{
				iy = 0;
			}
			else
			{
				iy = ( -n + 1 ) * incy;
			}

			for ( i = 0; i < n; i++ )
			{
				dy[iy] = dx[ix];
				ix = ix + incx;
				iy = iy + incy;
			}

		}

		return;
	}
	//****************************************************************************80

	public static double ddot ( int n, double dx[], int dx_first, int incx, double dy[], int dy_first, int incy )
	{
		double dtemp;
		int i;
		int ix;
		int iy;
		int m;

		dtemp = 0.0;

		if ( n <= 0 )
		{
			return dtemp;
		}
		//
		//  Code for unequal increments or equal increments
		//  not equal to 1.
		//
		if ( incx != 1 || incy != 1 )
		{
			if ( 0 <= incx )
			{
				ix = 0;
			}
			else
			{
				ix = ( - n + 1 ) * incx;
			}

			if ( 0 <= incy )
			{
				iy = 0;
			}
			else
			{
				iy = ( - n + 1 ) * incy;
			}

			for ( i = 0; i < n; i++ )
			{
				dtemp = dtemp + dx[ix+dx_first] * dy[iy+dy_first];
				ix = ix + incx;
				iy = iy + incy;
			}
		}
		//
		//  Code for both increments equal to 1.
		//
		else
		{
			m = n % 5;

			for ( i = 0; i < m; i++ )
			{
				dtemp = dtemp + dx[i+dx_first] * dy[i+dy_first];
			}

			for ( i = m; i < n; i = i + 5 )
			{
				dtemp = dtemp + dx[i+dx_first  ] * dy[i+dy_first  ]
						+ dx[i+1+dx_first] * dy[i+1+dy_first]
								+ dx[i+2+dx_first] * dy[i+2+dy_first]
										+ dx[i+3+dx_first] * dy[i+3+dy_first]
												+ dx[i+4+dx_first] * dy[i+4+dy_first];
			}

		}

		return dtemp;
	}
	//****************************************************************************80

	public static double dmach ( int job )
	{
		double eps;
		double huge;
		double s;
		double tiny;
		double value;

		eps = 1.0;
		for ( ; ; )
		{
			value = 1.0 + ( eps / 2.0 );
			if ( value <= 1.0 )
			{
				break;
			}
			eps = eps / 2.0;
		}

		s = 1.0;

		for ( ; ; )
		{
			tiny = s;
			s = s / 16.0;

			if ( s * 1.0 == 0.0 )
			{
				break;
			}

		}

		tiny = ( tiny / eps ) * 100.0;
		huge = 1.0 / tiny;

		if ( job == 1 )
		{
			value = eps;
		}
		else if ( job == 2 )
		{
			value = tiny;
		}
		else if ( job == 3 )
		{
			value = huge;
		}
		else
		{
			System.err.println("\n"+"DMACH - Fatal error!\n"+"  Illegal input value of JOB = "+job+"\n");
			System.exit ( 1 );
		}

		return value;
	}
	//****************************************************************************80

	public static double dnrm2 ( int n, double x[], int x_first, int incx )
	{
		double absxi;
		int i;
		int ix;
		double norm;
		double scale;
		double ssq;
		double value;

		if ( n < 1 || incx < 1 )
		{
			norm = 0.0;
		}
		else if ( n == 1 )
		{
			norm = r8_abs ( x[0+x_first] );
		}
		else
		{
			scale = 0.0;
			ssq = 1.0;
			ix = 0;

			for ( i = 0; i < n; i++ )
			{
				if ( x[ix+x_first] != 0.0 )
				{
					absxi = r8_abs ( x[ix+x_first] );
					if ( scale < absxi )
					{
						ssq = 1.0 + ssq * ( scale / absxi ) * ( scale / absxi );
						scale = absxi;
					}
					else
					{
						ssq = ssq + ( absxi / scale ) * ( absxi / scale );
					}
				}
				ix = ix + incx;
			}

			norm  = scale * Math.sqrt ( ssq );
		}

		return norm;
	}
	//****************************************************************************80

	public static void drot ( int n, double x[], int x_first, int incx, double y[], int y_first, int incy, double c,
			double s )
	{
		int i;
		int ix;
		int iy;
		double stemp;

		if ( n <= 0 )
		{
		}
		else if ( incx == 1 && incy == 1 )
		{
			for ( i = 0; i < n; i++ )
			{
				stemp = c * x[i+x_first] + s * y[i+y_first];
				y[i+y_first]  = c * y[i+y_first] - s * x[i+x_first];
				x[i+x_first]  = stemp;
			}
		}
		else
		{
			if ( 0 <= incx )
			{
				ix = 0;
			}
			else
			{
				ix = ( - n + 1 ) * incx;
			}

			if ( 0 <= incy )
			{
				iy = 0;
			}
			else
			{
				iy = ( - n + 1 ) * incy;
			}

			for ( i = 0; i < n; i++ )
			{
				stemp = c * x[ix+x_first] + s * y[iy+y_first];
				y[iy+y_first] = c * y[iy+y_first] - s * x[ix+x_first];
				x[ix+x_first] = stemp;
				ix = ix + incx;
				iy = iy + incy;
			}

		}

		return;
	}
	//****************************************************************************80
	//x[0-4] = double *sa, double *sb, double *c, double *s
	public static void drotg ( double x[] )
	{
		double r;
		double roe;
		double scale;
		double z;

		if ( r8_abs ( x[1] ) < r8_abs ( x[0] ) )
		{
			roe = x[0];
		}
		else
		{
			roe = x[1];
		}

		scale = r8_abs ( x[0] ) + r8_abs ( x[1] );

		if ( scale == 0.0 )
		{
			x[2] = 1.0;
			x[3] = 0.0;
			r = 0.0;
		}
		else
		{
			r = scale * Math.sqrt ( ( x[0] / scale ) * ( x[0] / scale )
					+ ( x[1] / scale ) * ( x[1] / scale ) );
			r = r8_sign ( roe ) * r;
			x[2] = x[0] / r;
			x[3] = x[1] / r;
		}

		if ( 0.0 < r8_abs ( x[2] ) && r8_abs ( x[2] ) <= x[3] )
		{
			z = 1.0 / x[2];
		}
		else
		{
			z = x[3];
		}

		x[0] = r;
		x[1] = z;

		return;
	}
	//****************************************************************************80

	public static void dscal ( int n, double sa, double x[], int x_first, int incx )
	{
		int i;
		int ix;
		int m;

		if ( n <= 0 )
		{
		}
		else if ( incx == 1 )
		{
			m = n % 5;

			for ( i = 0; i < m; i++ )
			{
				x[i+x_first] = sa * x[i+x_first];
			}

			for ( i = m; i < n; i = i + 5 )
			{
				x[i+x_first]   = sa * x[i+x_first];
				x[i+1+x_first] = sa * x[i+1+x_first];
				x[i+2+x_first] = sa * x[i+2+x_first];
				x[i+3+x_first] = sa * x[i+3+x_first];
				x[i+4+x_first] = sa * x[i+4+x_first];
			}
		}
		else
		{
			if ( 0 <= incx )
			{
				ix = 0;
			}
			else
			{
				ix = ( - n + 1 ) * incx;
			}

			for ( i = 0; i < n; i++ )
			{
				x[ix+x_first] = sa * x[ix+x_first];
				ix = ix + incx;
			}

		}

		return;
	}
	//****************************************************************************80

	public static void dswap ( int n, double x[], int x_first, int incx, double y[], int y_first, int incy )
	{
		int i;
		int ix;
		int iy;
		int m;
		double temp;

		if ( n <= 0 )
		{
		}
		else if ( incx == 1 && incy == 1 )
		{
			m = n % 3;

			for ( i = 0; i < m; i++ )
			{
				temp = x[i+x_first];
				x[i+x_first] = y[i+y_first];
				y[i+y_first] = temp;
			}

			for ( i = m; i < n; i = i + 3 )
			{
				temp = x[i+x_first];
				x[i+x_first] = y[i+y_first];
				y[i+y_first] = temp;

				temp = x[i+1+x_first];
				x[i+1+x_first] = y[i+1+y_first];
				y[i+1+y_first] = temp;

				temp = x[i+2+x_first];
				x[i+2+x_first] = y[i+2+y_first];
				y[i+2+y_first] = temp;
			}
		}
		else
		{
			if ( 0 <= incx )
			{
				ix = 0;
			}
			else
			{
				ix = ( - n + 1 ) * incx;
			}

			if ( 0 <= incy )
			{
				iy = 0;
			}
			else
			{
				iy = ( - n + 1 ) * incy;
			}

			for ( i = 0; i < n; i++ )
			{
				temp = x[ix+x_first];
				x[ix+x_first] = y[iy+y_first];
				y[iy+y_first] = temp;
				ix = ix + incx;
				iy = iy + incy;
			}

		}

		return;
	}
	//****************************************************************************80

	public static int i4_max ( int i1, int i2 )
	{
		int value;

		if ( i2 < i1 )
		{
			value = i1;
		}
		else
		{
			value = i2;
		}
		return value;
	}
	//****************************************************************************80

	public static int i4_min ( int i1, int i2 )
	{
		int value;

		if ( i1 < i2 )
		{
			value = i1;
		}
		else
		{
			value = i2;
		}
		return value;
	}
	//****************************************************************************80

	public static int idamax ( int n, double dx[], int dx_first, int incx )
	{
		double dmax;
		int i;
		int ix;
		int value;

		value = 0;

		if ( n < 1 || incx <= 0 )
		{
			return value;
		}

		value = 1;

		if ( n == 1 )
		{
			return value;
		}

		if ( incx == 1 )
		{
			dmax = r8_abs ( dx[0+dx_first] );

			for ( i = 1; i < n; i++ )
			{
				if ( dmax < r8_abs ( dx[i+dx_first] ) )
				{
					value = i + 1;
					dmax = r8_abs ( dx[i+dx_first] );
				}
			}
		}
		else
		{
			ix = 0;
			dmax = r8_abs ( dx[0+dx_first] );
			ix = ix + incx;

			for ( i = 1; i < n; i++ )
			{
				if ( dmax < r8_abs ( dx[ix+dx_first] ) )
				{
					value = i + 1;
					dmax = r8_abs ( dx[ix+dx_first] );
				}
				ix = ix + incx;
			}
		}

		return value;
	}
	//****************************************************************************80

	public static boolean lsame ( char ca, char cb )
	{
		if ( ca == cb )
		{
			return true;
		}

		if ( 'A' <= ca && ca <= 'Z' )
		{
			if ( ca - 'A' == cb - 'a' )
			{
				return true;
			}
		}
		else if ( 'a' <= ca && ca <= 'z' )
		{
			if ( ca - 'a' == cb - 'A' )
			{
				return true;
			}
		}

		return false;
	}
	//****************************************************************************80

	public static double r8_abs ( double x )
	{
		double value;

		if ( 0.0 <= x )
		{
			value = x;
		}
		else
		{
			value = -x;
		}
		return value;
	}
	//****************************************************************************80

	public static double r8_max ( double x, double y )
	{
		double value;

		if ( y < x )
		{
			value = x;
		}
		else
		{
			value = y;
		}
		return value;
	}
	//****************************************************************************80

	public static double r8_sign ( double x )
	{
		double value;

		if ( x < 0.0 )
		{
			value = -1.0;
		}
		else
		{
			value = 1.0;
		}
		return value;
	}
	//****************************************************************************80

	public static void xerbla ( String srname, int info )
	{
		System.err.println("\n"+"XERBLA - Fatal error!\n"+"  On entry to routine "+srname+"\n"
				+"  input parameter number " + info + " had an illegal value.\n");
		System.exit ( 1 );
	}
}
