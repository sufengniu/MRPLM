package affy.qualityControl;

import java.util.Arrays;

public class MedianSteps {
	public static double median(double [] x, int length){
		double med;
		double [] buffer = x.clone();

		//quickSort(buffer, 0, length-1);
		Arrays.sort(buffer);
		
		int half = (length + 1)/2;
		if (length % 2 == 1){
			med = buffer[half - 1];
		} else {
			med = (buffer[half] + buffer[half-1])/2.0;
		}

		return med;
	}

	public static double median_nocopy(double [] x, int length){
		int i;
		int half;
		double med;

		//quickSort(x, 0, length-1);
		Arrays.sort(x);
		
		half = (length + 1)/2;
		if (length % 2 == 1){
			med = x[half - 1];
		} else {
			med = (x[half] + x[half-1])/2.0;
		}

		return med;
	}

	public static double median_nocopy_hasNA(double []x, int length, int num_na){
		int i;
		int half;
		double med;

		//quickSort(x, 0, length-1);
		Arrays.sort(x);
		
		half = (length - num_na+ 1)/2;
		if (length % 2 == 1){
			med = x[half - 1];
		} else {
			med = (x[half] + x[half-1])/2.0;
		}

		return med;
	}

	//Q[0]: LQ, Q[1]: UQ
	public static double quartiles(double [] x, int length, double [] Q){

		double lowindex = (double)(length -1)*0.25;
		double highindex = (double)(length -1)*0.75;

		double lowfloor = Math.floor(lowindex);
		double highfloor = Math.floor(highindex);

		double lowceil = Math.ceil(lowindex);
		double highceil = Math.ceil(highindex);

		boolean low_i = lowindex > lowfloor;
		boolean high_i = highindex > highfloor;

		double qslow = x[(int)lowfloor];
		double qshigh = x[(int)highfloor];

		double low_h = lowindex - lowfloor;
		double high_h = highindex - highfloor;

		if (low_h > 1e-10){
			qslow = (1.0 - low_h)*qslow + low_h*x[(int)lowceil];
		}
		if (high_h > 1e-10){
			qshigh = (1.0 - high_h)*qshigh + high_h*x[(int)highceil];
		}


		Q[1] = qshigh;
		Q[0] = qslow;

		return qshigh - qslow;

	}
	
	/*private static void quickSort(double arr[], int left, int right) {
		int index = partition(arr, left, right);
		if (left < index - 1)
			quickSort(arr, left, index - 1);
		if (index < right)
			quickSort(arr, index, right);
	}

	private static int partition(double arr[], int left, int right){
		int i = left, j = right;
		double tmp;
		double pivot = arr[(left + right) / 2];

		while (i <= j) {
			while (arr[i] < pivot)
				i++;
			while (arr[j] > pivot)
				j--;
			if (i <= j) {
				tmp = arr[i];
				arr[i] = arr[j];
				arr[j] = tmp;
				i++;
				j--;
			}
		};
		return i;
	}*/
}
