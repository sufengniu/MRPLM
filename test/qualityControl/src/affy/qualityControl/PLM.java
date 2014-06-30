package affy.qualityControl;

public class PLM {
	
    MedianSteps medianSteps = new MedianSteps();

    void PLMsummarize(double[][] z, int numprobes, int numchips, double[] results, double[] SEresults, double[] affinities) {

        double[] weights = new double[numprobes * numchips];
        double[] resids = new double[numprobes * numchips];
        double[] y = new double[numprobes * numchips];
        double[] beta = new double[numprobes + numchips];
        double[] se = new double[numprobes + numchips];

        double[] varcov = null;

        for (int j = 0; j < numchips; j++) {
            for (int i = 0; i < numprobes; i++) {
                resids[j * numprobes + i] = y[j * numprobes + i] = z[i][j];
            }
        }

	    jniWrapper jni = new jniWrapper();

	    jni.rlm_fit_anova(y, numprobes, numchips, beta, resids, weights, 1.345, 20);
        //rlm_fit_anova(y, numprobes, numchips, beta, resids, weights, 1.345, 20);
	    jni.rlm_compute_se_anova(numprobes, numchips, resids, weights, se, varcov, 1.345);
	

/*
        rlm_fit_anova(y, numprobes, numchips, beta, resids, weights, 1.345, 20);
	
	
        rlm_compute_se_anova(numprobes, numchips, resids, weights, se, varcov, 1.345);
*/
        for (int j = 0; j < numchips; j++) {
            results[j] = beta[j];
            SEresults[j] = se[j];
        }

        for (int j = 0; j < numchips; j++) {
            for (int i = 0; i < numprobes; i++) {
                z[i][j] = resids[j * numprobes + i];
            }
        }
    }

    void rlm_fit_anova(double[] y, int y_rows, int y_cols, double[] out_beta, double[] resids,
                       double[] weights, double psi_k, int max_iter) {

        double acc = 1e-4;
        double scale = 0.0;

        double[] old_resids = new double[y_rows * y_cols];
        double[] rowmeans = new double[y_rows];

        double[] xtwx = new double[(y_rows + y_cols - 1) * (y_rows + y_cols - 1)];
        double[] xtwy = new double[y_rows + y_cols];

        int rows = y_rows * y_cols;

/* intially use equal weights */
        for (int i = 0; i < rows; i++) {
            weights[i] = 1.0;
        }


/* starting matrix */
/* for (int i = 0; i < y_rows; i++) {
for (int j = 0; j < y_cols; j++) {
resids[j * y_rows + i] = y[j * y_rows + i];
}
}*/


/* sweep columns (ie chip effects) */
        for (int j = 0; j < y_cols; j++) {
            out_beta[j] = 0.0;
            double sumweights = 0.0;
            for (int i = 0; i < y_rows; i++) {
                out_beta[j] += weights[j * y_rows + i] * resids[j * y_rows + i];
                sumweights += weights[j * y_rows + i];
            }
            out_beta[j] /= sumweights;
            for (int i = 0; i < y_rows; i++) {
                resids[j * y_rows + i] = resids[j * y_rows + i] - out_beta[j];
            }
        }

/* sweep rows (ie probe effects) */
        for (int i = 0; i < y_rows; i++) {
            rowmeans[i] = 0.0;
            double sumweights = 0.0;
            for (int j = 0; j < y_cols; j++) {
                rowmeans[i] += weights[j * y_rows + i] * resids[j * y_rows + i];
                sumweights += weights[j * y_rows + i];
            }
            rowmeans[i] /= sumweights;
            for (int j = 0; j < y_cols; j++) {
                resids[j * y_rows + i] = resids[j * y_rows + i] - rowmeans[i];
            }
        }
        for (int i = 0; i < y_rows - 1; i++) {
            out_beta[i + y_cols] = rowmeans[i];
        }

        for (int iter = 0; iter < max_iter; iter++) {

            scale = med_abs(resids, rows) / 0.6745;
            if (Math.abs(scale) < 1e-10) {
                /*printf("Scale too small \n"); */
                break;
            }
            for (int i = 0; i < rows; i++) {
                old_resids[i] = resids[i];
            }

            for (int i = 0; i < rows; i++) {
                weights[i] = PsiFunction.huber(resids[i] / scale, psi_k, 0); /* psi_huber(resids[i]/scale,k,0); */
            }

/* weighted least squares */

            for (int i = 0; i < xtwx.length; i++) {
                xtwx[i] = 0.0;
            }

            /***************** GPU offload part ***************/

        XTWX(y_rows, y_cols, weights, xtwx);


        XTWXinv(y_rows, y_cols, xtwx);///////////////////call cholesky decomposition
        XTWY(y_rows, y_cols, weights, y, xtwy);

        for (int i = 0; i < y_rows + y_cols - 1; i++) {
    out_beta[i] = 0.0;
    for (int j = 0; j < y_rows + y_cols - 1; j++) {
    out_beta[i] += xtwx[j * (y_rows + y_cols - 1) + i] * xtwy[j];
    }
}
    //        new jniWrapper().wlsAcc(weights, y, out_beta, y_rows, y_cols);


/* residuals */
            for (int i = 0; i < y_rows - 1; i++) {
                for (int j = 0; j < y_cols; j++) {
                    resids[j * y_rows + i] = y[j * y_rows + i] - (out_beta[j] + out_beta[i + y_cols]);
                }
            }

            for (int j = 0; j < y_cols; j++) {
                double endprobe = 0.0;
                for (int i = 0; i < y_rows - 1; i++) {
                    endprobe += out_beta[i + y_cols];
                }
                resids[j * y_rows + y_rows - 1] = y[j * y_rows + y_rows - 1] - (out_beta[j] - endprobe);
            }

/*check convergence based on residuals */
            double conv = irls_delta(old_resids, resids, rows);

            if (conv < acc) {
                /* printf("Converged \n");*/
                break;

            }
        }
    }
    private double med_abs(double[] x, int length) {
        double[] buffer = new double[length];

        for (int i = 0; i < length; i++)
            buffer[i] = Math.abs(x[i]);
        //med_abs =
        return MedianSteps.median(buffer, length);
    }

    private double irls_delta(double[] old, double[] newv, int length) {
        double sum = 0.0;
        double sum2 = 0.0;
        double divisor = 1e-20;

        for (int i = 0; i < length; i++) {
            sum = sum + (old[i] - newv[i]) * (old[i] - newv[i]);
            sum2 = sum2 + old[i] * old[i];
        }

        if (sum2 >= divisor) {
            divisor = sum2;
        }

        return Math.sqrt(sum / divisor);
    }

    private void XTWY(int y_rows, int y_cols, double[] wts, double[] y, double[] xtwy) {

/*sweep columns (ie chip effects)*/
        for (int j = 0; j < y_cols; j++) {
            xtwy[j] = 0.0;
            for (int i = 0; i < y_rows; i++) {
                xtwy[j] += wts[j * y_rows + i] * y[j * y_rows + i];
            }
        }

/*sweep rows (ie probe effects)*/
        for (int i = y_rows - 1; i >= 0; i--) {
            xtwy[i + y_cols] = 0.0;
            for (int j = 0; j < y_cols; j++) {
                xtwy[i + y_cols] += wts[j * y_rows + i] * y[j * y_rows + i];
            }
            if (i < y_rows - 1) {
                xtwy[i + y_cols] = xtwy[i + y_cols] - xtwy[y_cols + y_rows - 1];
            }
        }
    }

    static void XTWX(int y_rows, int y_cols, double[] wts, double[] xtwx) {

        int Msize = y_cols + y_rows - 1;

/* diagonal elements of first part of matrix ie upper partition*/
/* for (int j = 0; j < y_cols; j++) {
for (int i = 0; i < y_rows-1; i++) {
xtwx[j * Msize + j] += wts[j * y_rows + i];
xtwx[(y_cols + i) * Msize + (y_cols + i)] += wts[j * y_rows + i];
for (int k = i; k < y_rows - 1; k++) {
xtwx[(y_cols + k) * Msize + (y_cols + i)] = xtwx[(y_cols + i) * Msize + (y_cols + k)] += wts[j * y_rows + (y_rows - 1)];
}
xtwx[j * Msize + (y_cols + i)] = xtwx[(y_cols + i) * Msize + j] = wts[j * y_rows + i] - wts[j * y_rows + (y_rows - 1)];
}
xtwx[j * Msize + j] += wts[j * y_rows + y_rows-1];
}*/
        double sum = 0;
        for (int j = 0; j < y_cols; j++) {
            for (int i = 0; i < y_rows - 1; i++) {
                xtwx[j * Msize + j] += wts[j * y_rows + i];
                xtwx[(y_cols + i) * Msize + (y_cols + i)] += wts[j * y_rows + i];
                xtwx[j * Msize + (y_cols + i)] = xtwx[(y_cols + i) * Msize + j] = wts[j * y_rows + i] - wts[j * y_rows + (y_rows - 1)];

            }
            sum += wts[j * y_rows + (y_rows - 1)];
            xtwx[j * Msize + j] += wts[j * y_rows + y_rows - 1];
        }

        for (int i = 0; i < y_rows - 1; i++) {
            for (int k = i; k < y_rows - 1; k++) {
                xtwx[(y_cols + k) * Msize + (y_cols + i)] = xtwx[(y_cols + i) * Msize + (y_cols + k)] += sum;
            }
        }

    }

    static void XTWXinv(int y_rows, int y_cols, double[] xtwx) {
        int Msize = y_cols + y_rows - 1;
        double[] P = new double[y_cols];
        double[] RP = new double[y_cols * (y_rows - 1)];
        double[] RPQ = new double[(y_rows - 1) * (y_rows - 1)];
        double[] S = new double[(y_rows - 1) * (y_rows - 1)];
        double[] work = new double[(y_rows - 1) * (y_rows - 1)];

        for (int j = 0; j < y_cols; j++) {
            for (int i = 0; i < y_rows - 1; i++) {
                RP[j * (y_rows - 1) + i] = xtwx[j * Msize + (y_cols + i)] * (1.0 / xtwx[j * Msize + j]);
            }
        }

        for (int i = 0; i < y_rows - 1; i++) {
            for (int j = i; j < y_rows - 1; j++) {
                for (int k = 0; k < y_cols; k++) {
                    RPQ[j * (y_rows - 1) + i] += RP[k * (y_rows - 1) + j] * xtwx[k * Msize + (y_cols + i)];
                }
                RPQ[i * (y_rows - 1) + j] = RPQ[j * (y_rows - 1) + i];
                RPQ[j * (y_rows - 1) + i] = RPQ[i * (y_rows - 1) + j]
                        = xtwx[(y_cols + i) * Msize + (y_cols + j)] - RPQ[i * (y_rows - 1) + j];
            }
        }

    /* Lets start making the inverse */
        MatrixFunctions.Choleski_inverse(RPQ, S, work, y_rows - 1, false);
        for (int j = 0; j < y_cols; j++) {
            for (int i = 0; i < y_rows - 1; i++) {
                xtwx[j * Msize + (y_cols + i)] = 0.0;
                for (int k = 0; k < y_rows - 1; k++) {
                    xtwx[j * Msize + (y_cols + i)] += -1.0 * (S[i * (y_rows - 1) + k]) * RP[j * (y_rows - 1) + k];
                }
                xtwx[(y_cols + i) * Msize + j] = xtwx[j * Msize + (y_cols + i)];
            }

            P[j] = 1.0 / xtwx[j * Msize + j];

            for (int i = j; i < y_cols; i++) {
                xtwx[i * Msize + j] = 0.0;
                for (int k = 0; k < y_rows - 1; k++) {
                    xtwx[i * Msize + j] += RP[i * (y_rows - 1) + k] * xtwx[j * Msize + (y_cols + k)];
                }
                xtwx[i * Msize + j] *= -1.0;
                xtwx[j * Msize + i] = xtwx[i * Msize + j];
            }
            xtwx[j * Msize + j] += P[j];
        }


        for (int j = 0; j < y_rows - 1; j++) {
            for (int i = 0; i < y_rows - 1; i++) {
                xtwx[(y_cols + j) * Msize + (y_cols + i)] = S[j * (y_rows - 1) + i];
            }
        }
    }

}
