package affy.mas5;

import affymetrix.calvin.exception.UnsignedOutOfLimitsException;
import affymetrix.fusion.cdf.*;
import affymetrix.fusion.cel.FusionCELData;
import affymetrix.gcos.cdf.GeneChipProbeSetType;

import java.io.IOException;
import java.util.*;


public class ExpressionAlgorithm {

    private final int EXP_UNKNOWN_COMP_CALL_TYPE = 0;
    private final int EXP_INCREASE_CALL_TYPE = 1;
    private final int EXP_DECREASE_CALL_TYPE = 2;
    private final int EXP_INCREASE_MODERATE_CALL_TYPE = 3;
    private final int EXP_DECREASE_MODERATE_CALL_TYPE = 4;
    private final int EXP_NO_CHANGE_CALL_TYPE = 5;
    private final int EXP_NO_COMP_CALL_TYPE = 6;

    private final int EXP_PRESENT_CALL_TYPE = 0;
    private final int EXP_MARGINAL_CALL_TYPE = 1;
    private final int EXP_ABSENT_CALL_TYPE = 2;
    private final int EXP_NO_ABS_CALL_TYPE = 3;

    //this should be extern
    //private String ExpressionAbsCallString;


    //private List<ControlInformationType> ControlInformationList;

    ////////////////////////////////////////////////////////////////////

    //private String ExpressionCompCallString;

    private FusionCELData m_Cell;
    private FusionCELData m_Baseline;
    private FusionCDFData m_Cdf;
    private List<Double> m_tTable;
    private Map<String, Integer> m_ProbeSetNames;
    private String m_Error;
    private String m_LibPath;
    private int m_NumResults;
    private int m_NumFeatures;
    private boolean m_bCompDataexists = false;
    private List<AbsStatExpressionProbeSetResultType> m_AbsStatResults;
    private List<AbsStatExpressionProbeSetResultType> m_BaselineAbsStatResults;
    private List<CompStatExpressionProbeSetResultType> m_CompStatResults;
    private ExpStatAlgSettings m_Params;
    private double m_RelNF;
    private double m_RelNF2;
    private double m_RawQ;
    private double m_BaselineRawQ;
    private Map<Integer, Boolean> maskedProbes;
    ////////////////////////////////////////////////////////////////////////
    private double presentCalls = 0;
    private double allCalls = 0;
    /////////////////////////////////////////////////////////////////////////
    // For storing the average/stdv/min/max background
    private AvgStdvMinMaxType m_BgStats;

    // For storing the average/stdv/min/max noise
    private AvgStdvMinMaxType m_NoiseStats;

    // For storing computed background zone information
    private AllZonesInfoType m_ZonesInfo;

    // The list of control info.
    private List<ControlInformationType> m_ControlInfo;

    // Do we want to write a cel file
    private boolean m_WriteCell = false;

    private PTable pTable;

    public ExpressionAlgorithm() throws IOException {
        m_Cell = new FusionCELData();
        m_Baseline = new FusionCELData();
        m_Cdf = new FusionCDFData();
        m_tTable = new ArrayList<Double>();
        m_ProbeSetNames = new HashMap<String, Integer>();
        m_AbsStatResults = new ArrayList<AbsStatExpressionProbeSetResultType>();
        m_BaselineAbsStatResults = new ArrayList<AbsStatExpressionProbeSetResultType>();
        m_CompStatResults = new ArrayList<CompStatExpressionProbeSetResultType>();
        maskedProbes = new HashMap<Integer, Boolean>();
        m_ControlInfo = new ArrayList<ControlInformationType>();
        m_Params = new ExpStatAlgSettings();
        m_BgStats = new AvgStdvMinMaxType();
        m_NoiseStats = new AvgStdvMinMaxType();
        m_ZonesInfo = new AllZonesInfoType();
        pTable = new PTable();
    }




/*	// The control information.
    private List<ControlInformationType> getControlInfo() { return m_ControlInfo; }

	// The raw noise value (RawQ)
	private double getRawQ() { return m_RawQ; }

	// The raw noise value (RawQ) for the baseline file.
	private double getBaselineRawQ() { return m_BaselineRawQ; }

	// The background stats
	private AvgStdvMinMaxType getBgStats() { return m_BgStats; }

	// The noise stats
	private AvgStdvMinMaxType getNoiseStats() { return m_NoiseStats; }


	// The number of rows of features
	private int getRows() { return m_Cdf.getHeader().getRows(); }

	// The number of cols of features
	private int getCols() { return m_Cdf.getHeader().getCols(); }

	// The chip type
	private String getChipType() { return m_Cdf.getChipType(); }

	// The CEL object.
	private FusionCELData getCellData() { return m_Cell; }

	// The baseline CEL object.
	private FusionCELData getBaselineCellData() { return m_Baseline; }


*/
    //////////////////////////////////////////////////////////////////////

    double logtwo(double value) {
        return (Math.log(value) / Math.log(2.0));
    }

    //////////////////////////////////////////////////////////////////////

    double antiLog(double value) {
        return (Math.pow(2.0, value));
    }

    //////////////////////////////////////////////////////////////////////

    double getDegreeOfFreedom(int nProbeSet) {
        double df = 0.7 * (nProbeSet - 1);
        if (df < 1.0) {
            df = 1;
        }

        return df;
    }

    //////////////////////////////////////////////////////////////////////

    double vGetT(final int nAtoms, final List<Double> tTable, double level) {
        double t = 0.0;
        if (nAtoms == 0) {
            return t;
        } else if (nAtoms <= tTable.size()) {
            t = tTable.get(nAtoms - 1);
        } else {
            double degreeF = 0.7 * (nAtoms - 1);
            t = MathLib.tCDFinversed(level, degreeF);
        }
        return t;
    }

    //////////////////////////////////////////////////////////////////////

	/*	int SortdoublesAscending(final void elem1, final void elem2)
	{
		double v1 = (double) elem1;
		double v2 = (double) elem2;

		if (v1 > v2)
			return 1;
		else if (v2 > v1)
			return -1;
		else
			return 0;
	}*/

    //////////////////////////////////////////////////////////////////////

    double computeSquaredDistance(double x1, double y1, double x2, double y2) {
        double diffx = x1 - x2;
        double diffy = y1 - y2;
        return (diffx * diffx + diffy * diffy);
    }

    //////////////////////////////////////////////////////////////////////

    double ComputeWeightAtXY(double x, double y, double centerX, double centerY, double smoothFactor) {
        return 1.0 / (computeSquaredDistance(x, y, centerX, centerY) + smoothFactor);
    }

    //////////////////////////////////////////////////////////////////////

    double trimmedInterpolation(double bwGM, double bLow, double bHigh, double left, double right) {
        double x = left;
        if (bwGM >= bHigh) {
            x = right;
        } else if (bwGM > bLow) {
            double wG = (bwGM - bLow) / (bHigh - bLow);
            x = wG * right + (1.0 - wG) * left;
        }
        return x;
    }

    //////////////////////////////////////////////////////////////////////

    double OneStepBiweightAlgorithm(final List<Double> x, double c, double epsilon) {
        if (x.size() == 0)
            return 0.0;

        double medianValue = median(x);
        double MAD = medianAbsoluteDeviation(x) * c + epsilon;
        int n = x.size();
        double value = 0.0;
        double weightedSumNumer = 0.0;
        double weightedSumDenom = 0.0;

        for (int i = 0; i < n; i++) {
            double diff = x.get(i) - medianValue;
            double u = diff / MAD;
            double uSquare = u * u;
            double oneMinusUSqaure = 1.0 - uSquare;
            if (Math.abs(u) < 1.0) {
                weightedSumNumer += diff * oneMinusUSqaure * oneMinusUSqaure;
                weightedSumDenom += oneMinusUSqaure * oneMinusUSqaure;
            }
        }

        if (weightedSumDenom != 0.0) {
            value = medianValue + weightedSumNumer / weightedSumDenom;
        }

        return value;
    }

    //////////////////////////////////////////////////////////////////////

    double UncertaintyOfEstimate(final List<Double> x, double c, double epsilon) {
        if (x.size() == 0)
            return 0.0;

        double medianValue = median(x);
        double MAD = medianAbsoluteDeviation(x) * c + epsilon;
        int n = x.size();
        double value = 0.0;
        double numer = 0.0;
        double denom = 0.0;

        for (int i = 0; i < n; i++) {
            double diff = x.get(i) - medianValue;
            double u = diff / MAD;
            double uSquare = u * u;
            double oneMinusUSquare = 1.0 - uSquare;
            if (Math.abs(u) < 1.0) {
                numer += diff * diff * Math.pow(oneMinusUSquare, 4);
                denom += oneMinusUSquare * (1.0 - 5.0 * uSquare);
            }
        }
        numer = Math.sqrt(numer);
        denom = Math.abs(denom);

        if (denom != 0.0)
            value = numer / denom;

        return value;
    }

    //////////////////////////////////////////////////////////////////////

    double ComputeEstIntenDiff(final List<Double> v, double stp) {
        double so = 0;
        if (v.size() == 0)
            return so;

        if (v.size() == 1) {
            so = v.get(0);
        } else {
            double meanV = mean(v);
            double s = 0;
            int i;
            for (i = 0; i < v.size(); ++i) {
                s += (v.get(i) - meanV) * (v.get(i) - meanV);
            }
            double std = Math.sqrt(s / (v.size() - 1));
            double high = meanV + stp * std;
            double low = meanV - stp * std;
            s = 0;
            int ctr = 0;
            for (i = 0; i < v.size(); ++i) {
                if (v.get(i) <= high && v.get(i) >= low) {
                    ctr++;
                    s += v.get(i);
                }
            }
            so = s / (double) ctr;
        }
        return so;
    }

    //////////////////////////////////////////////////////////////////////

    void formPTable(final int nT, List<ArrayList<Double>> vvP) {
        List<Integer> mask = new ArrayList<Integer>(nT);

        for (int n = 0; n < nT; ++n) {
            mask.add(1 << n);
        }
        for (int n = 1; n <= nT; ++n) {
            int twoToN = 1 << n;
            List<Integer> posRanks = new ArrayList<Integer>(twoToN);
            for (int i = 0; i < twoToN; ++i) {
                int sum = 0;
                for (int j = 0; j < n; ++j) {
                    if ((i & mask.get(j)) != 0)
                        sum += j + 1;
                }
                posRanks.add(sum);
            }

            for (int i = 0; i < twoToN; ++i) {
                //List<Integer> signedRanks = new ArrayList<Integer>(n);
                int w = 0;
                int myCode = 0;
                for (int j = 0; j < n; ++j) {
                    if ((i & mask.get(j)) != 0) {
                        w += j + 1;
                        myCode += 1 << j;
                    }
                }
                double tail = 0;
                for (int j = 0; j < twoToN; ++j) {
                    if (posRanks.get(j) > w) {
                        tail++;
                    } else if (posRanks.get(j) == w) {
                        tail += 0.5;
                    }
                }
                vvP.get(n - 1).set(myCode, (tail / (double) twoToN));
            }
        }
    }

    /////////////////////////////////////////////////////////////////////////////
    //Template Functions
    /////////////////////////////////////////////////////////////////////////////
    ExpResults oneSidedSignRank2(final List<Double> x, final double alpha) throws IOException {
        // int len =  x.size(); // unused pk
        // 1. Ignore all zero differences.
        double[] newdiff = new double[x.size()];
        int n = 0;
        int i;
        for (i = 0; i < x.size(); ++i) {
            if (x.get(i) != 0.0) {
                newdiff[n] = (x.get(i));
                n++;
            }
        }
        if (n == 0) // No non-zero differences.  Output 0.5 as the one-sided p-value and detection is absent.
            return new ExpResults(0.5, 0);
        //newdiff.setSize(n);

        // 2.  Assign integer ranks to the differences.
        double[] ranks = new double[n];
        for (i = 0; i < n; ++i) {
            ranks[i] = (double) (i + 1);
        }

        // 3. Convert differences to absolute values and sort in ascending order.
        List<ExpResults> absdiff = new ArrayList<ExpResults>(n);
        for (i = 0; i < n; ++i) {
            ExpResults expResults = new ExpResults();
            expResults.p_value = Math.abs(newdiff[i]);
            expResults.call = i;
            absdiff.add(expResults);
        }
        Collections.sort(absdiff, new ExpResultsComparator());

        // 4. If there are ties among absolute differences, all differences in a tie
        //    group are assigned to a rank equal to the average of the integer ranks.
        int nTies = 0;
        // Avoid cross-platform compatibility problems from approximate doubleing point arithmetic.
        final double tiny = 2e-09;
        for (i = 0; i < n - 1; ++i) {
            if (Math.abs(absdiff.get(i).p_value - absdiff.get(i + 1).p_value) < tiny) {
                nTies++;
                break;
            }
        }
        //int tieGroup = 0;
        double doubleVarMod = 0; // modification of variance due to ties.
        if (nTies != 0) {
            i = 0;
            while (i < n - 1) {
                double initElement = absdiff.get(i).p_value;
                int tieGroupSize = 1;
                for (int j = i + 1; j < n; ++j) {
                    if (Math.abs(absdiff.get(j).p_value - initElement) < tiny) {
                        tieGroupSize++;
                        if (j == n - 1) {
                            i = j;
                            //tieGroup ++;
                            for (int m = j - tieGroupSize + 1; m <= j; ++m) {
                                ranks[m] = (2 * j - tieGroupSize + 3) / 2.0;
                            }
                            doubleVarMod += tieGroupSize *
                                    ((double) tieGroupSize * tieGroupSize - 1);
                            break;
                        }
                    } else {
                        i = j;
                        if (tieGroupSize > 1) {
                            //tieGroup ++;
                            for (int m = j - tieGroupSize; m <= j - 1; ++m) {
                                ranks[m] = (2 * j - tieGroupSize + 1) / 2.0;
                            }
                            doubleVarMod += tieGroupSize *
                                    ((double) tieGroupSize * tieGroupSize - 1);
                        }
                        break;
                    }
                }
            }
        }
        double[] invr = new double[n];

        for (i = 0; i < n; ++i) {
            invr[absdiff.get(i).call] = ranks[i];
        }

        double w = 0;
        for (i = 0; i < n; ++i) {
            if (newdiff[i] > 0)
                w += invr[i];
        }
        ExpResults ans = new ExpResults();
        if (n > 11) {
            // Use the asymptotic approximation:
            // S' = [S - n(n+1)/4]/sqrt[n(n+1)(2n+1)/24 - c]
            // where S is the sum of all positive signed ranks
            // and c = sum(b(b*b - 1)/48 is the modification of variance due to ties.
            //
            // p-value = 1 - f(S')
            // where f(S') is the cumulative distribution function of
            // standard normal distribution.
            double dw = w - ((double) n) * (n + 1) / 4.0;
            double denom2 = (((double) n) * (n + 1) * (2 * n + 1) - 0.5 * doubleVarMod) / 24.0;
            if (denom2 <= 0) {
                return new ExpResults(0, 0);
            }
            double z = dw / Math.sqrt(denom2);
            ans.p_value = 1 - MathLib.normalCDF(z);
        } else {
            if (nTies == 0) {
                int iCode = 0;
                for (i = 0; i < n; ++i) {
                    if (newdiff[i] > 0) {
                        iCode += 1 << ((int) invr[i] - 1);
                    }
                }
                ans.p_value = pTable.fGetPValue(n - 1, iCode);
            } else {
                int twoToN = 1 << n;
                int[] mask = new int[n];
                for (i = 0; i < n; ++i) {
                    mask[i] = 1 << i;
                }
                int[] posRanks = new int[twoToN];
                for (i = 0; i < twoToN; ++i) {
                    double sum = 0;
                    for (int j = 0; j < n; ++j) {
                        if ((i & mask[j]) != 0)
                            sum += ranks[j];
                    }
                    posRanks[i] = (int) sum;
                }
                double tail = 0;
                for (i = 0; i < twoToN; ++i) {
                    if (posRanks[i] > w) {
                        tail++;
                    } else if (posRanks[i] == w) {
                        tail += 0.5;
                    }
                }
                ans.p_value = tail / (double) twoToN;
            }
        }
        ans.call = (ans.p_value < alpha) ? 1 : 0;
        return ans;
    }

    //////////////////////////////////////////////////////////////////////

    ExpResults newSignRank(final List<Double> dif, final double alpha1, final double alpha2)
            throws IOException {
        ExpResults newPH = new ExpResults();
        ExpResults oldPH = oneSidedSignRank2(dif, alpha1);
        newPH.p_value = oldPH.p_value;
        if (oldPH.call == 1)
            newPH.call = 2;
        else if (oldPH.p_value < alpha2)
            newPH.call = 1;
        else if (oldPH.p_value > 1 - alpha1)
            newPH.call = -2;
        else if (oldPH.p_value > 1 - alpha2)
            newPH.call = -1;
        else
            newPH.call = 0;
        return newPH;
    }

    //////////////////////////////////////////////////////////////////////

    double mean(final List<Double> v) {
        if (v.size() == 0) {
            return -1;
        }
        double sum = 0;
        for (int i = 0; i < v.size(); ++i) {
            sum += v.get(i);
        }
        return (sum / (double) v.size());
    }

    //////////////////////////////////////////////////////////////////////

    double median(final List<Double> v) {
        List<Double> u = new ArrayList<Double>(v);

        int len = u.size();
        if (len < 1) {
            return -1;
        }
        Collections.sort(u);
        int half = len / 2;
        if (len % 2 == 1) {
            return u.get(half);
        } else {
            return (u.get(half - 1) + u.get(half)) / 2.0;
        }
    }

    //////////////////////////////////////////////////////////////////////

    public double stddev(final List<Double> v) {
        if (v.size() == 0 || v.size() == 1) {
            return -1.0;
        }
        double meanValue = mean(v);
        int n = v.size();
        double sum = 0.0;
        for (int i = 0; i < n; ++i) {
            sum += (v.get(i) - meanValue) * (v.get(i) - meanValue);
        }
        return Math.sqrt((1.0 / (n - 1)) * sum);
    }

    //////////////////////////////////////////////////////////////////////

    public double medianAbsoluteDeviation(final List<Double> x) {
        double meanValue = median(x);
        int size = x.size();
        List<Double> v = new ArrayList<Double>(size);
        for (int i = 0; i < size; i++) {
            v.add(Math.abs(x.get(i) - meanValue));
        }
        return median(v);
    }

    //////////////////////////////////////////////////////////////////////

    DoublePair trimMeanAndStd(List<Double> v, final double p1, final double p2) {
        int total = v.size();
        DoublePair fp = new DoublePair(0, 0);
        if (total > 0) {
            Collections.sort(v);
            int n1 = 0;
            int n2 = (int) Math.floor(total * p2);
            double subtotal = n2;
            if (subtotal > 2) {
                double sum = 0;
                int i;
                for (i = n1; i < n2; ++i) {
                    sum += v.get(i);
                }
                double tMean = sum / subtotal;
                sum = 0;
                for (i = n1; i < n2; ++i) {
                    sum += (v.get(i) - tMean) * (v.get(i) - tMean);
                }
                fp.value1 = tMean;
                fp.value2 = Math.sqrt(sum / (subtotal - 1));
            } else if (subtotal == 1) {
                fp.value1 = v.get(n1);
            }
        }
        return fp;
    }

    //////////////////////////////////////////////////////////////////////

    public double trimMean(final List<Double> vec, final double p1, final double p2) {
        List<Double> whole = new ArrayList<Double>(vec);
        int total = (int) whole.size();
        if (total == 0)
            return 0.0;

        Collections.sort(whole);

        double dG1 = total * p1;
        double dG2 = total * (1.0 - p2);
        int g1 = (int) Math.floor(dG1);
        int g2 = (int) Math.floor(dG2);
        double r1 = dG1 - g1;
        double r2 = dG2 - g2;
        int last = total - g2 - 1;
        if (last <= 0.0) { // it is theoretically impossible for last < 0, but
            last = 0; // we add the code here to guarantee proper bounds even if there is any numerical unstability
        }
        double sum = (1.0 - r1) * whole.get(g1) + (1.0 - r2) * whole.get(last);
        for (int i = g1 + 1; i < last; ++i) {
            sum += whole.get(i);
        }
        double subtotal = last - g1 - 1;
        if (subtotal <= 0.0) {
            subtotal = 0.0;
        }
        subtotal += 2.0 - r1 - r2;
        return sum / subtotal;
    }

    //////////////////////////////////////////////////////////////////////

	/*	CExpressionAlgorithmImplementation()
	{
		m_RawQ = 0.0;
		m_NumResults = 0;
		m_bCompDataexists = false;
		m_ZonesInfo.pZones = NULL;
		m_WriteCell = false;
		m_OutputDirectory = NULL;
	}*/

    //////////////////////////////////////////////////////////////////////

	/*	~CExpressionAlgorithmImplementation()
	{
		Clear();
	}*/

    //////////////////////////////////////////////////////////////////////

    void Clear() {
        maskedProbes.clear();
        m_ProbeSetNames.clear();
        m_NumResults = 0;
        m_bCompDataexists = false;
        //m_ControlInfo.clear();
        try {
            m_Cell.close();
            m_Baseline.close();
            //m_Cdf.close();
        } catch (Exception e) {
        }

		/*		ArrayList<AbsStatExpressionProbeSetResultType>::iterator abs;
		for(abs = m_AbsStatResults.begin(); abs != m_AbsStatResults.end(); ++abs)
			delete[] (*abs);
		m_AbsStatResults.erase(m_AbsStatResults.begin(), m_AbsStatResults.end());
		for(abs = m_BaselineAbsStatResults.begin(); abs != m_BaselineAbsStatResults.end(); ++abs)
			delete[] (*abs);
		m_BaselineAbsStatResults.erase(m_BaselineAbsStatResults.begin(), m_BaselineAbsStatResults.end());
		ArrayList<CompStatExpressionProbeSetResultType *>::iterator comp;
		for(comp = m_CompStatResults.begin(); comp != m_CompStatResults.end(); ++comp)
			delete[] (*comp);
		m_CompStatResults.erase(m_CompStatResults.begin(), m_CompStatResults.end());

		if (m_ZonesInfo.pZones)
			delete [] m_ZonesInfo.pZones;
		m_ZonesInfo.pZones = NULL;*/
    }

    //////////////////////////////////////////////////////////////////////

    String GetError() {
        return m_Error;
    }

    //////////////////////////////////////////////////////////////////////

    void SetWriteCell(final boolean writeCell) {
        m_WriteCell = writeCell;
    }

    //////////////////////////////////////////////////////////////////////

    void SetLibPath(final String libPath) {
        m_LibPath = libPath;
    }

    //////////////////////////////////////////////////////////////////////

    int GetNumResults() {
        return m_NumResults;
    }

    //////////////////////////////////////////////////////////////////////

    boolean DoesCompDataexists() {
        return m_bCompDataexists;
    }

    //////////////////////////////////////////////////////////////////////

    String getProbeSetName(int index) {
        return m_Cdf.getProbeSetName(index);
    }

    //////////////////////////////////////////////////////////////////////

    AbsStatExpressionProbeSetResultType GetAbsStatResult(int index) {
        return m_AbsStatResults.get(index);
    }

    //////////////////////////////////////////////////////////////////////

    AbsStatExpressionProbeSetResultType GetBaselineAbsStatResult(int index) {
        return m_BaselineAbsStatResults.get(index);
    }

    //////////////////////////////////////////////////////////////////////

    CompStatExpressionProbeSetResultType GetCompStatResult(int index) {
        return (m_bCompDataexists == false ? null : m_CompStatResults.get(index));
    }

    //////////////////////////////////////////////////////////////////////

    public ExpStatAlgSettings GetParameters() {
        return m_Params;
    }

    //////////////////////////////////////////////////////////////////////

    AllZonesInfoType GetBackgroundZoneInfo() {
        return m_ZonesInfo;
    }

    //////////////////////////////////////////////////////////////////////

    boolean IsMasked(FusionCELData cel, int x, int y) {
        if (cel.isMasked(x, y) == true)
            return true;

        if (maskedProbes.containsKey(cel.xyToIndex(x, y)))
            return true;

        return false;
    }

    //////////////////////////////////////////////////////////////////////

	/*	boolean readNormMSKFile()
	{
		if (m_Params.NormMethod == ExpStatAlgSettings.NormalizationOptionsEnum.NORM_TO_SELECTED_PROBE_SETS)
		{
			// read the MSK file.
			m_Params.NormGenes.clear();
			CMSKFileData msk;
			msk.setFileName(m_Params.NormMaskFile);
			if (msk.read() == false)
			{
				m_Error = "Unable to read the MSK file.";
				Clear();
				return false;
			}

			// Get the probe setes.
			ProbeSetListConstIt begin;
			ProbeSetListConstIt end;
			ProbeSetListConstIt it;
			msk.GetProbeSetIterators(begin, end);
			for (it=begin; it!=end; ++it)
			{
				m_Params.NormGenes.push_back(*it);
			}
		}
		return true;
	}*/

    //////////////////////////////////////////////////////////////////////

	/*	boolean readScaleMSKFile()
	{
		if (m_Params.SFMethod == ExpStatAlgSettings.ScalingOptionsEnum.SCALE_TO_SELECTED_PROBE_SETS)
		{
			// read the MSK file.
			m_Params.ScaleGenes.clear();
			CMSKFileData msk;
			msk.setFileName(m_Params.ScaleMaskFile);
			if (msk.read() == false)
			{
				m_Error = "Unable to read the MSK file.";
				Clear();
				return false;
			}

			// Get the probe setes.
			ProbeSetListfinalIt begin;
			ProbeSetListfinalIt end;
			ProbeSetListfinalIt it;
			msk.GetProbeSetIterators(begin, end);
			for (it=begin; it!=end; ++it)
			{
				m_Params.ScaleGenes.push_back(*it);
			}
		}
		return true;
	}*/

    //////////////////////////////////////////////////////////////////////

	/*	boolean readProbeMSKFile()
	{
		// Check if a MSK file was given.
		maskedProbes.clear();
		if (m_Params.ProbeMaskFile == 0)
			return true;

		// read the mask file.
		CMSKFileData msk;
		msk.setFileName(m_Params.ProbeMaskFile);
		if (msk.read() == false)
		{
			m_Error = "Unable to read the MSK file.";
			Clear();
			return false;
		}

		// Get each item from the probe pair mask section.
		FusionCDFProbeSetInformation unit;
		FusionCDFProbeGroupInformation blk;
		FusionCDFProbeInformation cel;
		int index;
		int xyindex;
		ProbeSetIndiciesListfinalIt begin;
		ProbeSetIndiciesListfinalIt end;
		ProbeSetIndiciesListfinalIt it;
		msk.GetProbeSetIndiciesIterators(begin, end);
		for (it=begin; it!=end; ++it)
		{
			// Find the corresponding probe set and mask the pair (both PM and MM probe)
			index = m_ProbeSetNames.get(it.probeSetName);
			m_Cdf.getProbeSetInformation(index, unit);
			unit.getGroupInformation(0, blk);
			for (list<int>::final_iterator pairIt=it.indicies.begin(); pairIt!=it.indicies.end(); ++pairIt)
			{
				index = 2*(*pairIt);
				if (index+1 >= blk.getNumCells())
					continue;

				// Mask the PM probe
				blk.getCell(index, cel);
				xyindex = m_Cell.XYToIndex(cel.getX(), cel.getY());
				maskedProbes.insert(make_pair(xyindex, true));

				// Mask the MM probe.
				blk.getCell(index+1, cel);
				xyindex = m_Cell.XYToIndex(cel.getX(), cel.getY());
				maskedProbes.insert(make_pair(xyindex, true));
			}
		}
		return true;
	}*/

    //////////////////////////////////////////////////////////////////////

    boolean readCelFile(final String celFile) {
        // read the CEL file
        m_Cell.setFileName(celFile);
        if (m_Cell.exists() == false) {
            m_Error = "The input CEL file does not exist.";
            Clear();
            return false;
        }
        if (m_Cell.read() == false) {
            m_Error = "Unable to read the input CEL file.";
            Clear();
            return false;
        }

        // Modify the saturated intensity if the CEL file was created from the HP scanner.
        //		if (IntensityFileType.FromHP(m_Cell) == true)
        //			m_Params.SaturatedIntensity = m_Params.HPSaturatedIntensity;

        return true;
    }

    //////////////////////////////////////////////////////////////////////

    boolean readBaselineCelFile(final String baseline) {
        // read the baseline CEL file
        if (!baseline.isEmpty()) {
            m_bCompDataexists = true;
            m_Baseline.setFileName(baseline);
            if (m_Baseline.exists() == false) {
                m_Error = "The baseline CEL file does not exist.";
                Clear();
                return false;
            }
            if (m_Baseline.read() == false) {
                m_Error = "Unable to read the baseline CEL file.";
                Clear();
                return false;
            }
        }
        return true;
    }

    //////////////////////////////////////////////////////////////////////

    boolean readCdfFile(final String cdfFile) {
        // read the CDF file.
        m_Cdf.setFileName(cdfFile);
        if (m_Cdf.exists() == false) {
            m_Error = "The associated library file (CDF) does not exist.";
            Clear();
            return false;
        }
        if (m_Cdf.read() == false) {
            m_Error = "Unable to read the library file.";
            Clear();
            return false;
        }

        // Check CDF file for all expression units.
        int iunit;
        m_NumResults = m_Cdf.getHeader().getNumProbeSets();
        for (iunit = 0; iunit < m_NumResults; iunit++) {
            if (m_Cdf.getProbeSetType(iunit) != GeneChipProbeSetType.ExpressionProbeSetType) {
                m_Error = "Unable to run algorithm on non-Expression arrays.";
                Clear();
                return false;
            }
        }
        m_NumFeatures = m_Cdf.getHeader().getCols() * m_Cdf.getHeader().getRows();
        return true;
    }

    //////////////////////////////////////////////////////////////////////

    boolean readCdfFile(FusionCELData cel) {
        // read the CDF file.
        String cdfFile = m_LibPath + cel.getChipType() + ".CDF";
        m_Cdf.setFileName(cdfFile);
        if (m_Cdf.exists() == false) {
            m_Error = "The associated library file (CDF) does not exist.";
            Clear();
            return false;
        }
        if (m_Cdf.read() == false) {
            m_Error = "Unable to read the library file.";
            Clear();
            return false;
        }

        // Check CDF file for all expression units.
        int iunit;
        m_NumResults = m_Cdf.getHeader().getNumProbeSets();
        for (iunit = 0; iunit < m_NumResults; iunit++) {
            if (m_Cdf.getProbeSetType(iunit) != GeneChipProbeSetType.ExpressionProbeSetType) {
                m_Error = "Unable to run algorithm on non-Expression arrays.";
                Clear();
                return false;
            }
        }
        m_NumFeatures = m_Cdf.getHeader().getCols() * m_Cdf.getHeader().getRows();
        return true;
    }

    //////////////////////////////////////////////////////////////////////

    void AllocateMemoryForResults() {
        // Allocate memory for the results
        //m_AbsStatResults.setSize(m_NumResults);
        //if (m_bCompDataexists == true)
        //m_CompStatResults.setSize(m_NumResults);

        // Copy the probe set names and allocate memory.
        String name;
        for (int iunit = 0; iunit < m_NumResults; iunit++) {
            name = m_Cdf.getProbeSetName(iunit);
            m_ProbeSetNames.put(name, iunit);
            m_AbsStatResults.add(new AbsStatExpressionProbeSetResultType());
            if (m_bCompDataexists == true)
                m_CompStatResults.add(new CompStatExpressionProbeSetResultType());
        }
    }

    //////////////////////////////////////////////////////////////////////

    boolean CheckIntensityData() throws UnsignedOutOfLimitsException, IOException {
        // Check if there is any data left.
        if (AllMaskedOut(m_Cell)) {
            m_Error = "Unable to compute expression results. The data has all been masked.";
            Clear();
            return false;
        }
        if (IsBlankCellFile(m_Cell)) {
            m_Error = "Unable to compute expression results. The data is all zeros.";
            Clear();
            return false;
        }

        // Now check the baseline data.
        if (m_bCompDataexists == true) {
            if (AllMaskedOut(m_Baseline)) {
                m_Error = "Unable to compute expression results. The baseline data has all been masked.";
                Clear();
                return false;
            }
            if (IsBlankCellFile(m_Baseline)) {
                m_Error = "Unable to compute expression results. The baseline data is all zeros.";
                Clear();
                return false;
            }
        }
        return true;
    }

    //////////////////////////////////////////////////////////////////////

    public boolean RunStat(final String celFile, final String baseline, final String cdfFile) throws UnsignedOutOfLimitsException, IOException {
        // Clear old results.
        Clear();
        //m_Error.erase();

        // read the CEL file
        if (readCelFile(celFile) == false) {
            System.err.println("the cel file " + celFile + " is wrong!");
            return false;
        }
        // read the baseline CEL file
        if (readBaselineCelFile(baseline) == false) {
            System.err.println("the baseline cel file " + celFile + " is wrong!");
            return false;
        }

        // read the CDF file.
        if (readCdfFile(cdfFile) == false) {
            System.err.println("the cdf file " + cdfFile + " is wrong!");
            return false;
        }

        // Allocate memory for the results
        AllocateMemoryForResults();

        // read the MSK file
        //		if (readProbeMSKFile() == false)
        //			return false;

        // Check the intensity data
        if (CheckIntensityData() == false)
            return false;

        // read the MSK files for normalization and scaling probe sets.
        //		if (readNormMSKFile() == false || readScaleMSKFile() == false)
        //			return false;

        // Call the algorithms

        ComputeStat();
        return true;
    }

    //////////////////////////////////////////////////////////////////////

    void ComputeStat() throws UnsignedOutOfLimitsException, IOException {
        // Make the expression call.	
        for (int iunit = 0; iunit < m_NumResults; iunit++) {
            ComputeExpressionStat(m_Cell, iunit, m_AbsStatResults.get(iunit));
        }
        // Compute the scaled and adjusted intensities for perfect match and mis match cells
        // based on new expression algorithm by Earl Hubbell, Mark Durst.
        List<ArrayList<Boolean>> UseAtomE = new ArrayList<ArrayList<Boolean>>(m_NumResults);
        List<ArrayList<Double>> PMe = new ArrayList<ArrayList<Double>>(m_NumResults);        // ScaledAdjusted Intensity, SA
        List<ArrayList<Double>> MMe = new ArrayList<ArrayList<Double>>(m_NumResults);        // ScaledAdjusted Intensity, SA
        List<ArrayList<DoublePair>> BGe = new ArrayList<ArrayList<DoublePair>>(m_NumResults);        // Smoothing Background. (PM in value1, MM in value2)
        List<ArrayList<DoublePair>> NoiseE = new ArrayList<ArrayList<DoublePair>>(m_NumResults);    // Noise. (PM in value1, MM in value2)
        List<Double> featureIntene = new ArrayList<Double>();        // ScaledAdjusted Intensity by feature

        List<Double> avgMeasurementE = new ArrayList<Double>(m_NumResults);    // Average measurement.
        List<ArrayList<Double>> PVe = new ArrayList<ArrayList<Double>>(m_NumResults);


        for (int i = 0; i < m_NumResults; i++) {
            UseAtomE.add(new ArrayList<Boolean>());
            PMe.add(new ArrayList<Double>());
            MMe.add(new ArrayList<Double>());
            BGe.add(new ArrayList<DoublePair>());
            NoiseE.add(new ArrayList<DoublePair>());
            PVe.add(new ArrayList<Double>());
        }

        // Pass in member zone info variable for experiment
        ComputeScaledAdjustedIntensity(m_Cell, PMe, MMe, UseAtomE, BGe, NoiseE, m_ZonesInfo, featureIntene);


        // Compute the measurement and confidence.
        ComputeMeasurement(PMe, MMe, UseAtomE, avgMeasurementE, PVe);

        // Set the measurement to the chip file.
        SetMeasurement(BGe, avgMeasurementE, m_AbsStatResults);

        // Determine the scale factor based on the average measurement.
        m_Params.ScaleFactor = DetermineScaleFactor(m_AbsStatResults);
        if (m_Params.ScaleFactor != 1.0)
            ModifyIntensities(m_AbsStatResults, m_Params.ScaleFactor);

        // Compute the absolute statistics for the baseline CEL file.
        if (m_bCompDataexists == true) {
            // Make the expression call.	
            //ArrayList<AbsStatExpressionProbeSetResultType> baselineAbsStatResults;
            //m_BaselineAbsStatResults.setSize(m_NumResults);
            for (int iunit = 0; iunit < m_NumResults; iunit++) {
                m_BaselineAbsStatResults.add(new AbsStatExpressionProbeSetResultType());
                ComputeExpressionStat(m_Baseline, iunit, m_BaselineAbsStatResults.get(iunit));
            }

            // Compute the scaled and adjusted intensities for perfect match and mis match cells
            // based on new expression algorithm by Earl Hubbell, Mark Durst.
            List<ArrayList<Boolean>> UseAtomB = new ArrayList<ArrayList<Boolean>>(m_NumResults);
            List<ArrayList<Double>> PMb = new ArrayList<ArrayList<Double>>(m_NumResults);        // ScaledAdjusted Intensity, SA
            List<ArrayList<Double>> MMb = new ArrayList<ArrayList<Double>>(m_NumResults);        // ScaledAdjusted Intensity, SA
            List<ArrayList<DoublePair>> BGb = new ArrayList<ArrayList<DoublePair>>(m_NumResults);        // Smoothing Background. (PM in value1, MM in value2)
            List<ArrayList<DoublePair>> NoiseB = new ArrayList<ArrayList<DoublePair>>(m_NumResults);    // Noise. (PM in value1, MM in value2)
            List<Double> featureIntenb = new ArrayList<Double>();        // ScaledAdjusted Intensity by feature

            // Pass in local zone info variable for baseline
            AllZonesInfoType zonesinfo = new AllZonesInfoType();
            ComputeScaledAdjustedIntensity(m_Baseline, PMb, MMb, UseAtomB, BGb, NoiseB, zonesinfo, featureIntenb);
            //if (zonesinfo.pZones) delete[] zonesinfo.pZones;

            // Compute the measurement and confidence.
            List<Double> avgMeasurementB = new ArrayList<Double>(m_NumResults);    // Average measurement.
            List<ArrayList<Double>> PVb = new ArrayList<ArrayList<Double>>(m_NumResults);
            for (int i = 0; i < m_NumResults; i++) {
                PVb.add(new ArrayList<Double>());
            }
            ComputeMeasurement(PMb, MMb, UseAtomB, avgMeasurementB, PVb);

            // Set the measurement to the chip file.
            SetMeasurement(BGb, avgMeasurementB, m_BaselineAbsStatResults);

            // Determine the scale factor based on the average measurement.
            m_Params.BaseScaleFactor = DetermineScaleFactor(m_BaselineAbsStatResults);
            if (m_Params.BaseScaleFactor != 1.0)
                ModifyIntensities(m_BaselineAbsStatResults, m_Params.BaseScaleFactor);

            // Compute the comparison call.
            CallComparison(m_Params.ScaleFactor, m_Params.BaseScaleFactor, m_BaselineAbsStatResults, BGe, BGb);


            // Determine the normalization factor
            DetermineNormFactor(m_AbsStatResults, m_BaselineAbsStatResults);

            // Normalize the intensities and signal values.
            if (m_Params.NormFactor != 1.0)
                ModifyIntensities(m_AbsStatResults, m_Params.NormFactor);

            // Scale and Normalize the PV before computing fold change.
            for (int i = 0; i < m_NumResults; i++) {
                int nProbePair = PVe.get(i).size();
                for (int j = 0; j < nProbePair; j++) {
                    PVe.get(i).set(j, PVe.get(i).get(j) + (double) logtwo(m_Params.ScaleFactor) + (double) logtwo(m_Params.NormFactor));  // PV is log scale.
                    PVb.get(i).set(j, PVb.get(i).get(j) + (double) logtwo(m_Params.BaseScaleFactor) + (double) logtwo(1.0));
                }
            }


            // Compute the fold change.
            ComputeFoldChange(PVb, PVe, UseAtomB, UseAtomE);
        }

        // Compute some summary statistics
        ComputeAvgMaxMinStdBgNoise(BGe, NoiseE);

        RawQ rawQ = new RawQ();
        //rawQ.setDefaults();
        rawQ.setVerticalZones(m_Params.NumberVertZones);
        rawQ.setHorizontalZones(m_Params.NumberHorZones);
        rawQ.setPercentBG(m_Params.NumberBGCells);
        m_RawQ = rawQ.computeRawQ(m_Cell, m_Cdf);
        if (m_bCompDataexists == true)
            m_BaselineRawQ = rawQ.computeRawQ(m_Baseline, m_Cdf);

        ReportCornerControls(FusionGeneChipQCProbeSetType.CheckerboardPositiveQCProbeSetType);
        ReportCornerControls(FusionGeneChipQCProbeSetType.CheckerboardNegativeQCProbeSetType);
        ReportCornerControls(FusionGeneChipQCProbeSetType.CentralCrossPositiveQCProbeSetType);
        ReportCornerControls(FusionGeneChipQCProbeSetType.CentralCrossNegativeQCProbeSetType);

        // Optionally generate a cel file containing background corrected
        // experiment intensities
		/*		if (m_WriteCell)
		{
			CELFileWriter writeCell;
			writeCell.setFileName (m_Cell.getFileName());
			if (! writeCell.read())
			{
				System.out.println("Unable to read cel file.");
				return;
			}
			// ensure not mapped read-only
			writeCell.EnsureNotMmapped();
			for (int celIx = 0; celIx < m_NumFeatures; ++celIx)
				writeCell.SetIntensity (celIx, featureIntene.get(celIx));
			if (m_OutputDirectory.size()() == 0)
			{
				System.out.println("An output directory is required to write a cel file.");
				return;
			}
			String outCell = m_OutputDirectory+m_Cell.getFileName();
			writeCell.setFileName (outCell);
			if (! writeCell.WriteTextCel())
			{
				System.out.println("Error writing cel file.");
				System.out.println(writeCell.GetError());
			}
		}*/

    }

    public String print() {
        String output = "\tsfs: " + m_Params.ScaleFactor
                + "\tavbg: " + m_BgStats.avg
                + "\tpps: " + (presentCalls / allCalls) * 100;
        return output;
    }

    public String isLowQuality() {
        double pps = (presentCalls / allCalls) * 100;

        if ((m_Params.ScaleFactor < -2) || (m_Params.ScaleFactor > 2))
            return "sfs";
        if ((m_BgStats.avg < 20) || (m_BgStats.avg > 100))
            return "avbg";
        if(pps < 35)
            return "pps";
        return "good"; //:)
    }

    //////////////////////////////////////////////////////////////////////

    boolean DetermineRelCallNormFactor(double sfe, double sfb) throws UnsignedOutOfLimitsException, IOException {
        FusionCDFProbeSetInformation unit = new FusionCDFProbeSetInformation();
        FusionCDFProbeGroupInformation blk = new FusionCDFProbeGroupInformation();
        FusionCDFProbeInformation pmCell = new FusionCDFProbeInformation();
        FusionCDFProbeInformation mmCell = new FusionCDFProbeInformation();

        m_RelNF = 1.0;
        m_RelNF2 = 1.0;

        double absNFDiff = (double) Math.abs(m_Params.NormFactor - 1.0);
        if ((m_Params.NormMethod == ExpStatAlgSettings.NormalizationOptionsEnum.DEFINED_NORMALIZATION_FACTOR) && (absNFDiff > 0.001f)) {
            double baselineSF = sfb;
            double expSF = sfe;
            m_RelNF = m_Params.NormFactor * expSF / baselineSF;
            m_RelNF2 = m_RelNF;
            return true;
        }

        int numAtomsUsed;
        int numBaseAtomsUsed;
        int numUnitUsed = 0;
        //int numBaseUnitUsed = 0;

        // double factorE = sfe; // unused pk
        // double factorB = sfb; // unused pk

        // Loop over all of the units.
        int UnitsPerChip = m_Cdf.getHeader().getNumProbeSets();
        List<ArrayList<Double>> diffE = new ArrayList<ArrayList<Double>>(UnitsPerChip);
        List<ArrayList<Double>> diffB = new ArrayList<ArrayList<Double>>(UnitsPerChip);
        List<Double> estIntenDiffE = new ArrayList<Double>(UnitsPerChip);
        List<Double> estIntenDiffB = new ArrayList<Double>(UnitsPerChip);
        List<Double> estIntenPosDiffE = new ArrayList<Double>(UnitsPerChip);
        List<Double> estIntenPosDiffB = new ArrayList<Double>(UnitsPerChip);
        List<ArrayList<Double>> PMe = new ArrayList<ArrayList<Double>>(UnitsPerChip);
        List<ArrayList<Double>> PMb = new ArrayList<ArrayList<Double>>(UnitsPerChip);
        List<Double> estPMIntenE = new ArrayList<Double>(UnitsPerChip);
        List<Double> estPMIntenB = new ArrayList<Double>(UnitsPerChip);

        double stp = m_Params.STP;

        //String probeSetName;
        for (int iUnit = 0; iUnit < UnitsPerChip; iUnit++) {
            //probeSetName = m_Cdf.getProbeSetName(iUnit);
            if (m_Params.NormMethod == ExpStatAlgSettings.NormalizationOptionsEnum.NORM_TO_SELECTED_PROBE_SETS)// &&
                //this.UseUnitInNormFactor(probeSetName, m_Params.NormGenes) == false)
                continue;

            m_Cdf.getProbeSetInformation(iUnit, unit);
            unit.getGroup(0, blk);
            int numCells = blk.getNumCells();
            //int nAtoms = blk.getNumLists();
            //diffE.get(iUnit).setSize(nAtoms);
            //diffB.get(iUnit).setSize(nAtoms);
            //PMe.get(iUnit).setSize(nAtoms);
            //PMb.get(iUnit).setSize(nAtoms);
            numAtomsUsed = 0;
            numBaseAtomsUsed = 0;

            for (int iCell = 0; iCell < numCells; iCell += 2) {
                blk.getCell(iCell, pmCell);
                blk.getCell(iCell + 1, mmCell);

                // Check if atom is used in the exp chip.
                if (IsMasked(m_Cell, pmCell.getX(), pmCell.getY()) == false &&
                        IsMasked(m_Cell, mmCell.getX(), mmCell.getY()) == false) {
                    diffE.get(iUnit).add((double) (
                            m_Cell.getIntensity(pmCell.getX(), pmCell.getY()) -
                                    m_Cell.getIntensity(mmCell.getX(), mmCell.getY())));
                    PMe.get(iUnit).add((double) m_Cell.getIntensity(pmCell.getX(), pmCell.getY()));
                    numAtomsUsed++;
                }

                // Check if atom is used in the baseline chip.
                if (IsMasked(m_Baseline, pmCell.getX(), pmCell.getY()) == false &&
                        IsMasked(m_Baseline, mmCell.getX(), mmCell.getY()) == false) {
                    diffB.get(iUnit).add((double) (
                            m_Baseline.getIntensity(pmCell.getX(), pmCell.getY()) -
                                    m_Baseline.getIntensity(mmCell.getX(), mmCell.getY())));
                    PMb.get(iUnit).add((double) m_Baseline.getIntensity(pmCell.getX(), pmCell.getY()));
                    numBaseAtomsUsed++;
                }
            }

            //diffE.get(iUnit).setSize(numAtomsUsed);
            //diffB.get(iUnit).setSize(numBaseAtomsUsed);
            //PMe.get(iUnit).setSize(numAtomsUsed);
            //PMb.get(iUnit).setSize(numBaseAtomsUsed);
            if (numAtomsUsed > 0) {
                estIntenDiffE.add(ComputeEstIntenDiff(diffE.get(iUnit), stp));
                estIntenPosDiffE.add((estIntenDiffE.get(numUnitUsed) > 0.0) ? estIntenDiffE.get(numUnitUsed) : 0.0);
                estPMIntenE.add(ComputeEstIntenDiff(PMe.get(iUnit), stp));
                numUnitUsed++;
            }

            if (numBaseAtomsUsed > 0) {
                estIntenDiffB.add(ComputeEstIntenDiff(diffB.get(iUnit), stp));
                estIntenPosDiffB.add((estIntenDiffB.get(numUnitUsed) > 0.0) ? estIntenDiffB.get(numUnitUsed) : 0.0);
                estPMIntenB.add(ComputeEstIntenDiff(PMb.get(iUnit), stp));
                //numBaseUnitUsed++;
            }
        }

        //estIntenDiffE.setSize(numUnitUsed);
        //estIntenPosDiffE.setSize(numUnitUsed);
        //estIntenDiffB.setSize(numBaseUnitUsed);
        //estIntenPosDiffB.setSize(numBaseUnitUsed);
        //estPMIntenE.setSize(numUnitUsed);
        //estPMIntenB.setSize(numBaseUnitUsed);

        double p1 = m_Params.IntensityLowPercent / 100;
        double p2 = 1.0 - m_Params.IntensityHighPercent / 100;

        // Calculate the primary normalization factor for PM-MM
        double TMe = trimMean(estIntenDiffE, p1, p2);
        if (TMe <= 0.0) {
            TMe = trimMean(estIntenPosDiffE, p1, p2);
            if (TMe <= 0.0)
                return false;
        }

        double TMb = trimMean(estIntenDiffB, p1, p2);
        if (TMb <= 0.0) {
            TMb = trimMean(estIntenPosDiffB, p1, p2);
            if (TMb <= 0.0)
                return false;
        }

        this.m_RelNF = TMb / TMe;

        // Calculate the primary normalization factor for PM-B
        TMe = trimMean(estPMIntenE, p1, p2);
        TMb = trimMean(estPMIntenB, p1, p2);

        if (TMe > 0.0)
            this.m_RelNF2 = TMb / TMe;
        else
            return false;

        return true;
    }

    //////////////////////////////////////////////////////////////////////

    void CallComparison
            (double sfe, double sfb, List<AbsStatExpressionProbeSetResultType> baselineAbsStatResults,
             List<ArrayList<DoublePair>> BGe,
             List<ArrayList<DoublePair>> BGb
            ) throws UnsignedOutOfLimitsException, IOException {
        //int numAtomsUsed;
        //int numCommonAtoms;
        FusionCDFProbeSetInformation unit = new FusionCDFProbeSetInformation();
        FusionCDFProbeGroupInformation blk = new FusionCDFProbeGroupInformation();
        FusionCDFProbeInformation pmCell = new FusionCDFProbeInformation();
        FusionCDFProbeInformation mmCell = new FusionCDFProbeInformation();

        boolean bRelCallNormFactor = DetermineRelCallNormFactor(sfe, sfb);

        // double factorE = sfe; // unused pk
        // double factorB = sfb; // unused pk

        // Loop over all of the units.
        int UnitsPerChip = m_Cdf.getHeader().getNumProbeSets();
        for (int iUnit = 0; iUnit < UnitsPerChip; iUnit++) {
            CompStatExpressionProbeSetResultType pUnitCompResult = GetCompStatResult(iUnit);

            AbsStatExpressionProbeSetResultType pUnitAbsResult = GetAbsStatResult(iUnit);

            // AbsStatExpressionProbeSetResultType *pBaseUnitResult = baselineAbsStatResults[iUnit]; // unused pk

            // Initialize the number of atoms used.
            //numAtomsUsed = 0;
            //numCommonAtoms = 0;

            // Loop over all atoms to determine the number of probe pairs
            // for which the difference and ratio are significantly greater
            // than the baseline.
            int nAtoms = pUnitAbsResult.NumPairs;
            List<Double> pmB = new ArrayList<Double>(nAtoms);
            List<Double> pmE = new ArrayList<Double>(nAtoms);
            //List<Double> mmB = new ArrayList<Double> (nAtoms);
            //List<Double> mmE = new ArrayList<Double> (nAtoms);
            List<Double> PMinusBgB = new ArrayList<Double>(nAtoms);
            List<Double> PMinusBgE = new ArrayList<Double>(nAtoms);
            List<Double> diffB = new ArrayList<Double>(nAtoms);
            List<Double> diffE = new ArrayList<Double>(nAtoms);


            m_Cdf.getProbeSetInformation(iUnit, unit);
            unit.getGroup(0, blk);
            int numCells = blk.getNumCells();
            int ctr = 0;
            for (int iCell = 0; iCell < numCells; iCell += 2) {
                blk.getCell(iCell, pmCell);
                blk.getCell(iCell + 1, mmCell);

                if (IsMasked(m_Cell, pmCell.getX(), pmCell.getY()) == false &&
                        IsMasked(m_Cell, mmCell.getX(), mmCell.getY()) == false &&
                        IsMasked(m_Baseline, pmCell.getX(), pmCell.getY()) == false &&
                        IsMasked(m_Baseline, mmCell.getX(), mmCell.getY()) == false) {
                    //++numAtomsUsed;

                    double expMatchCellI = m_Cell.getIntensity(pmCell.getX(), pmCell.getY());
                    double expMismatchCellI = m_Cell.getIntensity(mmCell.getX(), mmCell.getY());
                    double baselineMatchCellI = m_Baseline.getIntensity(pmCell.getX(), pmCell.getY());
                    double baselineMismatchCellI = m_Baseline.getIntensity(mmCell.getX(), mmCell.getY());

                    // Scale back the cell intensities for comparison.
                    if (expMatchCellI < m_Params.SaturatedIntensity &&
                            expMismatchCellI < m_Params.SaturatedIntensity &&
                            baselineMatchCellI < m_Params.SaturatedIntensity &&
                            baselineMismatchCellI < m_Params.SaturatedIntensity) {
                        pmB.add(baselineMatchCellI);
                        double mm = baselineMismatchCellI;
                        // No need to scale down the background like what we did in MAS since
                        // the background is calculated on fly.
                        double backGround = BGb.get(iUnit).get(iCell / 2).value1;
                        PMinusBgB.add(pmB.get(ctr) - backGround);
                        diffB.add(pmB.get(ctr) - mm);

                        pmE.add(expMatchCellI);
                        mm = expMismatchCellI;
                        backGround = BGe.get(iUnit).get(iCell / 2).value1;
                        PMinusBgE.add(pmE.get(ctr) - backGround);
                        diffE.add(pmE.get(ctr) - mm);
                        ctr++;
                    }
                }
            }

            pUnitCompResult.NumCommonPairs = 0;

            if (ctr <= 0 || bRelCallNormFactor == false) {
                // All intensities are greater than the saturated intensity value.
                // No Call will be returned.
                pUnitCompResult.Change = EXP_NO_COMP_CALL_TYPE;
                pUnitCompResult.ChangePValue = 0;
            } else {
                //pmB.setSize(ctr);
                //mmB.setSize(ctr);
                //PMinusBgB.setSize(ctr);
                //diffB.setSize(ctr);
                //pmE.setSize(ctr);
                //	mmE.setSize(ctr);
                //PMinusBgE.setSize(ctr);
                //diffE.setSize(ctr);

                double gamma1 = m_Params.Gamma1H;
                double gamma2 = m_Params.Gamma2H;

                ComputeIntensityDependentSignificances(gamma1, gamma2, pmE, pmB);

                DetermineComparativeCall(pUnitCompResult, diffE, diffB, PMinusBgE, PMinusBgB, m_Params.Gamma1H, m_Params.Gamma2H);

                pUnitCompResult.NumCommonPairs = ctr;
            }
        }
    }

    //////////////////////////////////////////////////////////////////////

    void DetermineComparativeCall(CompStatExpressionProbeSetResultType pCompData,
                                  List<Double> diffE, List<Double> diffB,
                                  List<Double> PMinusBgE, List<Double> PMinusBgB,
                                  double gamma1, double gamma2) throws IOException {
        double oneMinusGamma1 = 1.0 - gamma1;
        double oneMinusGamma2 = 1.0 - gamma2;

        List<Double> multipliers = new ArrayList<Double>(3);
        multipliers.add(1.0 / m_Params.Perturbation);
        multipliers.add(1.0);
        multipliers.add(m_Params.Perturbation);

        List<Double> pValues = new ArrayList<Double>(3);
        int size = diffE.size();
        for (int i = 0; i < 3; i++) {
            double nf = multipliers.get(i) * m_RelNF;
            double nf2 = multipliers.get(i) * m_RelNF2;
            int size2 = 2 * size;
            List<Double> nDiff = new ArrayList<Double>(size2);
            for (int j = 0; j < size; j++) {
                nDiff.add(nf * diffE.get(j) - diffB.get(j));
            }
            for (int j = size; j < size2; j++) {
                nDiff.add(m_Params.CMultiplier *
                        (nf2 * PMinusBgE.get(j - size) - PMinusBgB.get(j - size)));
            }

            ExpResults result = newSignRank(nDiff, gamma1, gamma2);
            pValues.add(result.p_value);
        }


        int Decision = EXP_NO_CHANGE_CALL_TYPE;

        if (pValues.get(0) < gamma1 && pValues.get(1) < gamma1 && pValues.get(2) < gamma1) {
            Decision = EXP_INCREASE_CALL_TYPE;
        } else if (pValues.get(0) > oneMinusGamma1 && pValues.get(1) > oneMinusGamma1 && pValues.get(2) > oneMinusGamma1) {
            Decision = EXP_DECREASE_CALL_TYPE;
        } else if (pValues.get(0) < gamma2 && pValues.get(1) < gamma2 && pValues.get(2) < gamma2) {
            Decision = EXP_INCREASE_MODERATE_CALL_TYPE;
        } else if (pValues.get(0) > oneMinusGamma2 && pValues.get(1) > oneMinusGamma2 && pValues.get(2) > oneMinusGamma2) {
            Decision = EXP_DECREASE_MODERATE_CALL_TYPE;
        }

        pCompData.Change = (char) Decision;
        double fCriticalPValue = 0.5;
        if (pValues.get(0) < 0.5 && pValues.get(1) < 0.5 && pValues.get(2) < 0.5)
            fCriticalPValue = Math.max(pValues.get(0), Math.max(pValues.get(1), pValues.get(2)));
        else if (pValues.get(0) > 0.5 && pValues.get(1) > 0.5 && pValues.get(2) > 0.5)
            fCriticalPValue = Math.min(pValues.get(0), Math.min(pValues.get(1), pValues.get(2)));
        pCompData.ChangePValue = fCriticalPValue;
    }

    //////////////////////////////////////////////////////////////////////

    void ComputeIntensityDependentSignificances(double gamma1, double gamma2,
                                                List<Double> PMe, List<Double> PMb) {
        double p1 = m_Params.IntensityLowPercent / 100;
        double p2 = 1.0 - m_Params.IntensityHighPercent / 100;

        double TPb = trimMean(PMb, p1, p2);
        double TPe = trimMean(PMe, p1, p2);

        double bc = Math.sqrt(TPb * TPe);
        double bLow = bc * m_Params.BLCoef;
        double bHigh = bc * m_Params.BHCoef;

        int size = PMe.size();
        List<Double> gMean = new ArrayList<Double>(size);
        for (int i = 0; i < size; i++) {
            gMean.add(Math.sqrt(PMb.get(i) * PMe.get(i)));
        }
        double c = m_Params.TuningConstantCGammas;
        double epsilon = m_Params.EpsilonGammas;

        double bwGM = OneStepBiweightAlgorithm(gMean, c, epsilon);

        gamma1 = trimmedInterpolation(bwGM, bLow, bHigh, m_Params.Gamma1H, m_Params.Gamma1L);
        gamma2 = trimmedInterpolation(bwGM, bLow, bHigh, m_Params.Gamma2H, m_Params.Gamma2L);
    }

    //////////////////////////////////////////////////////////////////////

    void ComputeProbeLogRatio(List<ArrayList<Double>> PVb,
                              List<ArrayList<Double>> PVe, List<ArrayList<Boolean>> UseAtomB,
                              List<ArrayList<Boolean>> UseAtomE, List<ArrayList<Double>> PLR) {
        int nProbeSet = PVe.size();
        for (int i = 0; i < nProbeSet; i++) {
            int nProbePair = PVe.get(i).size();
            //PLR.get(i).setSize(nProbePair);
            //int nProbePairUsed = 0;
            for (int j = 0; j < nProbePair; j++) {
                if (UseAtomB.get(i).get(j) && UseAtomE.get(i).get(j)) {
                    PLR.get(i).add(PVe.get(i).get(j) - PVb.get(i).get(j));
                    //nProbePairUsed++;
                }
            }
            //PLR.get(i).setSize(nProbePairUsed);
        }
    }

    //////////////////////////////////////////////////////////////////////

    void ComputeAvgLogRatio
            (
                    List<ArrayList<Double>> PLR,
                    List<Double> avgLogRatio
            ) {
        int nProbeSet = PLR.size();
        double c = m_Params.TuningConstantCAvgLogRatio;
        double epsilon = m_Params.EpsilonAvgLogRatio;
        for (int i = 0; i < nProbeSet; i++) {
            avgLogRatio.add(OneStepBiweightAlgorithm(PLR.get(i), c, epsilon));
        }
    }

    //////////////////////////////////////////////////////////////////////

    void BuildTTable(double level) {
        final int T_TABLE_SIZE = 50;
        m_tTable.clear();
        //m_tTable.setSize(T_TABLE_SIZE);
        //m_tTable = new double [T_TABLE_SIZE];
        for (int i = 0; i < T_TABLE_SIZE; ++i) {
            double df = 0.7 * i;
            if (df < 1.0) {
                df = 1;
            }
            m_tTable.add(MathLib.tCDFinversed(level, df));
        }
    }

    //////////////////////////////////////////////////////////////////////

    void ComputeFoldChange(List<ArrayList<Double>> PVb, List<ArrayList<Double>> PVe,
                           List<ArrayList<Boolean>> UseAtomB, List<ArrayList<Boolean>> UseAtomE) {
        // Compute fold change
        int nProbeSet = PVe.size();
        List<ArrayList<Double>> PLR = new ArrayList<ArrayList<Double>>(nProbeSet);
        List<Double> avgLogRatio = new ArrayList<Double>(nProbeSet);

        ComputeProbeLogRatio(PVb, PVe, UseAtomB, UseAtomE, PLR);
        ComputeAvgLogRatio(PLR, avgLogRatio);

        List<Double> uncertainty = new ArrayList<Double>(nProbeSet);
        double avgLogRatioLow;
        double avgLogRatioHigh;
        double c = m_Params.TuningConstantCAvgLogRatio;
        double epsilon = m_Params.EpsilonAvgLogRatio;

        // Store the fold changes and confidence intervals to the chip file.
        int UnitsPerChip = m_Cdf.getHeader().getNumProbeSets();
        BuildTTable(m_Params.RelConfInterval);

        for (int iUnit = 0; iUnit < UnitsPerChip; iUnit++) {
            CompStatExpressionProbeSetResultType pUnitResult =
                    (CompStatExpressionProbeSetResultType) m_CompStatResults.get(iUnit);

            pUnitResult.SignalLogRatio = avgLogRatio.get(iUnit);

            int size = PLR.get(iUnit).size();
            double tvalue = vGetT(size, m_tTable, m_Params.RelConfInterval);

            uncertainty.add(UncertaintyOfEstimate(PLR.get(iUnit), c, epsilon));
            double confidence = tvalue * uncertainty.get(iUnit);
            avgLogRatioLow = avgLogRatio.get(iUnit) - confidence;
            avgLogRatioHigh = avgLogRatio.get(iUnit) + confidence;

            pUnitResult.SignalLogRatioLow = avgLogRatioLow;
            pUnitResult.SignalLogRatioHigh = avgLogRatioHigh;
        }
    }

    //////////////////////////////////////////////////////////////////////

    void DetermineNormFactor
            (List<AbsStatExpressionProbeSetResultType> expStatResults,
             List<AbsStatExpressionProbeSetResultType> baselineStatResults) {
        int iUnit;
        double normFactor = 0.0;

        // User defined norm factor.
        if (m_Params.NormMethod == ExpStatAlgSettings.NormalizationOptionsEnum.DEFINED_NORMALIZATION_FACTOR)
            return;


        // Loop over all of the units.
        int UnitsPerChip = m_Cdf.getHeader().getNumProbeSets();
        List<Double> expList = new ArrayList<Double>(UnitsPerChip);
        List<Double> baseList = new ArrayList<Double>(UnitsPerChip);
        int unitCount = 0;
        //int baseunitCount=0;
        //String probeSetName;
        for (iUnit = 0; iUnit < UnitsPerChip; iUnit++) {
            // Determine if unit should be used.
            //probeSetName = m_Cdf.getProbeSetName(iUnit);
            if (m_Params.NormMethod == ExpStatAlgSettings.NormalizationOptionsEnum.NORM_TO_ALL_PROBE_SETS)// ||
            //UseUnitInNormFactor(probeSetName, m_Params.NormGenes))
            {
                if (expStatResults.get(iUnit).NumUsedPairs != 0) {
                    expList.add(expStatResults.get(iUnit).Signal);
                    unitCount++;
                }

                if (baselineStatResults.get(iUnit).NumUsedPairs != 0) {
                    baseList.add(baselineStatResults.get(iUnit).Signal);
                    //baseunitCount++;
                }
            }
        }

        //expList.setSize(unitCount);
        //baseList.setSize(baseunitCount);

        // Compute the ratio of the intensities.
        double p1 = m_Params.IntensityLowPercent / 100;
        double p2 = 1.0 - m_Params.IntensityHighPercent / 100;
        double avgB = trimMean(baseList, p1, p2);
        double avgE = trimMean(expList, p1, p2);

        if (avgE != 0.0) {
            normFactor = avgB / avgE;
            double diffNF = Math.abs(normFactor - 1.0);
            if (diffNF < 0.000001f)
                normFactor = 1.0;
        } else
            normFactor = 1.0;

        // Store the norm factor
        if (unitCount == 0)
            normFactor = 1.0;

        m_Params.NormFactor = normFactor;
    }

    //////////////////////////////////////////////////////////////////////

    void ComputeExpressionStat(FusionCELData pCell, int iUnit, AbsStatExpressionProbeSetResultType pUnit)
            throws UnsignedOutOfLimitsException, IOException {
        FusionCDFProbeSetInformation unit = new FusionCDFProbeSetInformation();
        FusionCDFProbeGroupInformation blk = new FusionCDFProbeGroupInformation();
        FusionCDFProbeInformation pmCell = new FusionCDFProbeInformation();
        FusionCDFProbeInformation mmCell = new FusionCDFProbeInformation();

        m_Cdf.getProbeSetInformation(iUnit, unit);
        unit.getGroup(0, blk);
        int numCells = blk.getNumCells();

        // Loop over the cells in the unit

        List<Double> discMinusTau = new ArrayList<Double>(blk.getNumLists());
        int ctr = 0;
        int iCell;
        for (iCell = 0; iCell < numCells; iCell += 2) {
            blk.getCell(iCell, pmCell);
            blk.getCell(iCell + 1, mmCell);

            if (IsMasked(pCell, pmCell.getX(), pmCell.getY()) == false &&
                    IsMasked(pCell, mmCell.getX(), mmCell.getY()) == false) {
                double pmI = pCell.getIntensity(pmCell.getX(), pmCell.getY());
                double mmI = pCell.getIntensity(mmCell.getX(), mmCell.getY());

                // Exclude saturated probe pair.
                if (mmI < m_Params.SaturatedIntensity) {
                    double sum = pmI + mmI;
                    double tau1 = m_Params.Tau;
                    if (sum > 0.0)
                        discMinusTau.add(((pmI - mmI) / sum) - tau1);
                    else
                        discMinusTau.add(-tau1);
                    ctr++;
                }
            }
        }
        pUnit.NumPairs = blk.getNumLists();

        // Compute the absolute call.
        if (ctr <= 0) {
            pUnit.Detection = EXP_NO_ABS_CALL_TYPE;
            pUnit.Signal = -1.0;
            pUnit.DetectionPValue = 0.0;
            pUnit.NumUsedPairs = 0;
        } else {
            //discMinusTau.setSize(ctr);
            ExpResults result = new ExpResults();
            //double medianValue = 0.0;
            if (discMinusTau.size() > 0) {
                //medianValue = m_Params.Tau + median(discMinusTau);
                result = newSignRank(discMinusTau, m_Params.Alpha1, m_Params.Alpha2);
            } else {
                result.call = 2;
                result.p_value = 0.0;
            }
            pUnit.DetectionPValue = result.p_value;
            ComputeAbsoluteCall(pUnit, result);
            pUnit.NumUsedPairs = ctr;
        }
    }

    //////////////////////////////////////////////////////////////////////

    boolean AllMaskedOut(FusionCELData pCell) {
        int nUnits;
        int iUnit;
        int nCells;
        int iCell;
        FusionCDFProbeSetInformation unit = new FusionCDFProbeSetInformation();
        FusionCDFProbeGroupInformation blk = new FusionCDFProbeGroupInformation();
        FusionCDFProbeInformation cell = new FusionCDFProbeInformation();

        // Get number of units.
        nUnits = m_Cdf.getHeader().getNumProbeSets();

        // Loop over all units.
        for (iUnit = 0; iUnit < nUnits; ++iUnit) {
            if (m_Cdf.getProbeSetType(iUnit) != GeneChipProbeSetType.ExpressionProbeSetType)
                continue;

            m_Cdf.getProbeSetInformation(iUnit, unit);
            unit.getGroup(0, blk);
            nCells = blk.getNumCells();
            for (iCell = 0; iCell < nCells; iCell++) {
                blk.getCell(iCell, cell);
                if (IsMasked(pCell, cell.getX(), cell.getY()) == false)
                    return false;
            }
        }
        return true;
    }

    //////////////////////////////////////////////////////////////////////

    boolean IsBlankCellFile(FusionCELData pCell) throws UnsignedOutOfLimitsException, IOException {
        final double MINIMUM_CELL_INTENSITY = 0.00001f;
        boolean bResult = true;

        int iUnit;
        int nCells;
        int iCell;
        FusionCDFProbeSetInformation unit = new FusionCDFProbeSetInformation();
        FusionCDFProbeGroupInformation blk = new FusionCDFProbeGroupInformation();
        FusionCDFProbeInformation cell = new FusionCDFProbeInformation();


        // Loop over all of the units.
        int UnitsPerChip = m_Cdf.getHeader().getNumProbeSets();
        for (iUnit = 0; iUnit < UnitsPerChip; iUnit++) {
            if (m_Cdf.getProbeSetType(iUnit) != GeneChipProbeSetType.ExpressionProbeSetType)
                continue;

            m_Cdf.getProbeSetInformation(iUnit, unit);
            unit.getGroup(0, blk);
            nCells = blk.getNumCells();
            for (iCell = 0; iCell < nCells; iCell++) {
                blk.getCell(iCell, cell);
                if (pCell.getIntensity(cell.getX(), cell.getY()) > MINIMUM_CELL_INTENSITY) {
                    bResult = false;
                    break;
                }
            }
        }

        return bResult;
    }

    //////////////////////////////////////////////////////////////////////

    void ComputeScaledAdjustedIntensity(
            FusionCELData pCell,
            List<ArrayList<Double>> PM,
            List<ArrayList<Double>> MM,
            List<ArrayList<Boolean>> UseAtom,
            List<ArrayList<DoublePair>> BG,
            List<ArrayList<DoublePair>> Noise,
            AllZonesInfoType ZonesInfo,
            List<Double> FeatureIntensity) throws UnsignedOutOfLimitsException, IOException {
        int iUnit;
        int CellsRemaining;
        int zonex;
        int zoney;
        int NumberZones;
        // int iInten=0; // unused pk

        //ArrayList<Integer> NumberCellsPerZone = null;

        // Determine the number of remaining cells in the vertical direction.
        CellsRemaining = m_Cdf.getHeader().getCols() % m_Params.NumberVertZones;
        if (CellsRemaining != 0) {
            zonex = (m_Cdf.getHeader().getCols() +
                    (int) m_Params.NumberVertZones - CellsRemaining) /
                    (int) m_Params.NumberVertZones;
        } else {
            zonex = m_Cdf.getHeader().getCols() / m_Params.NumberVertZones;
        }

        // Determine the number of remaining cells in the horizontal direction.
        CellsRemaining = m_Cdf.getHeader().getRows() % m_Params.NumberHorZones;
        if (CellsRemaining != 0) {
            zoney = (m_Cdf.getHeader().getRows() +
                    m_Params.NumberHorZones - CellsRemaining) /
                    m_Params.NumberHorZones;
        } else {
            zoney = m_Cdf.getHeader().getRows() / m_Params.NumberHorZones;
        }

        // Ensure that there are a match and mismatch cell in the same zone.
        zoney += zoney % 2; //EXPRESSION_ATOMS_PER_CELL;

        // Determine the total number of zones.
        NumberZones = (int) m_Params.NumberVertZones * m_Params.NumberHorZones;

        // Get number of units.
        int NumUnits = m_Cdf.getHeader().getNumProbeSets();

        // Allocate space for all atoms intensities and ID's.
        int[] NumberCellsPerZone = new int[NumberZones];

        // Clear arrays.
        //memset(NumberCellsPerZone, 0, sizeof(int)*NumberZones);

        // Loop over all units to determine the zone ID's and intensities.
        List<Vector<Double>> ZoneCells = new ArrayList<Vector<Double>>(NumberZones);
        for (int i = 0; i < NumberZones; i++)
            ZoneCells.add(new Vector<Double>());
        FusionCDFProbeSetInformation unit = new FusionCDFProbeSetInformation();
        FusionCDFProbeGroupInformation blk = new FusionCDFProbeGroupInformation();
        FusionCDFProbeInformation cell = new FusionCDFProbeInformation();
        boolean bMasked;
        for (iUnit = 0; iUnit < NumUnits; ++iUnit) {
            m_Cdf.getProbeSetInformation(iUnit, unit);

            // Only process expression units.
            if (unit.getProbeSetType() == GeneChipProbeSetType.ExpressionProbeSetType) {
                unit.getGroup(0, blk);
                int numCells = blk.getNumCells();

                // Loop over the atoms in the unit
                for (int iCell = 0; iCell < numCells; iCell++) {
                    blk.getCell(iCell, cell);
                    bMasked = IsMasked(pCell, cell.getX(), cell.getY());
                    if (bMasked == false) {
                        int nZone = DetermineZone(cell.getX(), cell.getY(), zonex, zoney);
                        ZoneCells.get(nZone).setSize(ZoneCells.get(nZone).size() + 1);
                        ZoneCells.get(nZone).add(0.0);
                        ZoneCells.get(nZone).set(NumberCellsPerZone[nZone],
                                (double) pCell.getIntensity(cell.getX(), cell.getY()));

                        if (nZone >= 0 && nZone < NumberZones)
                            NumberCellsPerZone[nZone]++;
                    }
                }
            }
        }

        // Allocate zones, set smooth factor and set num zones
        ZonesInfo.pZones = new ArrayList<ZoneInfo>(NumberZones);
        for (int i = 0; i < NumberZones; i++)
            ZonesInfo.pZones.add(new ZoneInfo());
        ZonesInfo.number_zones = NumberZones;
        ZonesInfo.smooth_factor = m_Params.SmoothFactorBG;

        // compute background for each zone
        for (int iZone = 0; iZone < NumberZones; iZone++) {
            // Compute the center coordinates of each zone.
            // (x1,y1) is the upper left corner
            // (x2,y2) is the lower right corner
            double x1 = ((int) (iZone % m_Params.NumberVertZones)) * zonex;
            double y1 = ((int) (iZone / m_Params.NumberVertZones)) * zoney;
            double x2 = x1 + zonex;
            double y2 = y1 + zoney;
            ZonesInfo.pZones.get(iZone).center.x = (x1 + x2) / 2;
            ZonesInfo.pZones.get(iZone).center.y = (y1 + y2) / 2;

            //int iCell=0;
            int numCell = NumberCellsPerZone[iZone];
            ZonesInfo.pZones.get(iZone).numCell = numCell;

            List<Double> zoneI = new ArrayList<Double>(numCell);
            //List<Integer> rank = new ArrayList<Integer>(numCell);

            for (int i = 0; i < numCell; i++) {
                double inten = ZoneCells.get(iZone).get(i);
                zoneI.add(ModifyIntensitySlightly(inten));
                //iCell++;
            }
            double lowBG = 0.0;
            double highBG = m_Params.NumberBGCells / 100.0;
            DoublePair fp = trimMeanAndStd(zoneI, lowBG, highBG);
            ZonesInfo.pZones.get(iZone).background = fp.value1;
            ZonesInfo.pZones.get(iZone).noise = fp.value2;
        }
        // End of computing background intensity and noise for each zone.
        // Carried zones and NumberZones as the required information.

        // Compute b(x,y), n(x,y), SA(x,y) which was stored in PM[i][j] and MM[i][j]
        double smoothF = m_Params.SmoothFactorBG;

        FusionCDFProbeInformation pmcell = new FusionCDFProbeInformation();
        FusionCDFProbeInformation mmcell = new FusionCDFProbeInformation();
        for (iUnit = 0; iUnit < NumUnits; ++iUnit) {
            m_Cdf.getProbeSetInformation(iUnit, unit);

            // Only process expression units.

            if (unit.getProbeSetType() == GeneChipProbeSetType.ExpressionProbeSetType) {
                unit.getGroup(0, blk);
                int numCells = blk.getNumCells();
                int nAtoms = blk.getNumLists();

                //PM.get(iUnit).setSize(nAtoms);
                //MM.get(iUnit).setSize(nAtoms);
                //UseAtom.get(iUnit).setSize(nAtoms);

                //BG.get(iUnit).setSize(nAtoms);
                //Noise.get(iUnit).setSize(nAtoms);

                for (int i = 0; i < nAtoms; i++) {
                    //PM.get(iUnit).add(0.0);
                    //MM.get(iUnit).add(0.0);
                    UseAtom.get(iUnit).add(true);
                    BG.get(iUnit).add(new DoublePair());
                    Noise.get(iUnit).add(new DoublePair());
                }

                // Loop over the atoms in the unit
                //int atomCount = 0;
                for (int iCell = 0, iAtom = 0; iCell < numCells; iCell += 2, ++iAtom) {
                    blk.getCell(iCell, pmcell);
                    blk.getCell(iCell + 1, mmcell);
                    bMasked = (IsMasked(pCell, pmcell.getX(), pmcell.getY()) == true) ||
                            (IsMasked(pCell, mmcell.getX(), mmcell.getY()) == true);
                    if (bMasked == true) {
                        UseAtom.get(iUnit).set(iAtom, false);
                    }

                    // Set the background (with smoothing adjustment) for matchcell.
                    Coordinate Cellxy = new Coordinate();
                    Cellxy.x = pmcell.getX();
                    Cellxy.y = pmcell.getY();
                    double WeightedSumBg = 0.0;
                    double WeightedSumNoise = 0.0;
                    double WeightedSumDenom = 0.0;
                    double background = 0.0;
                    double noise = 0.0;
                    int k;
                    for (k = 0; k < NumberZones; k++) {
                        WeightedSumBg += ComputeWeightAtXY(Cellxy.x, Cellxy.y, ZonesInfo.pZones.get(k).center.x, ZonesInfo.pZones.get(k).center.y, smoothF) * ZonesInfo.pZones.get(k).background;
                        WeightedSumNoise += ComputeWeightAtXY(Cellxy.x, Cellxy.y, ZonesInfo.pZones.get(k).center.x, ZonesInfo.pZones.get(k).center.y, smoothF) * ZonesInfo.pZones.get(k).noise;
                        WeightedSumDenom += ComputeWeightAtXY(Cellxy.x, Cellxy.y, ZonesInfo.pZones.get(k).center.x, ZonesInfo.pZones.get(k).center.y, smoothF);
                    }
                    if (WeightedSumDenom != 0.0) {
                        background = WeightedSumBg / WeightedSumDenom;
                        noise = WeightedSumNoise / WeightedSumDenom;
                    }

                    BG.get(iUnit).get(iAtom).value1 = background;
                    Noise.get(iUnit).get(iAtom).value1 = noise;

                    //double smoothFactor;
                    double scaledAdjustedI;
                    double inten;
                    double modifiedI;

                    inten = pCell.getIntensity((int) Cellxy.x, (int) Cellxy.y);
                    modifiedI = ModifyIntensitySlightly(inten);
                    scaledAdjustedI = ComputeAdjustedIntensity(modifiedI, background, noise);
                    PM.get(iUnit).add(scaledAdjustedI);

                    //////////////////////////////////////////////////////
                    // Compute Mis-match intensities
                    //////////////////////////////////////////////////////
                    Cellxy.x = mmcell.getX();
                    Cellxy.y = mmcell.getY();
                    WeightedSumBg = 0.0;
                    WeightedSumNoise = 0.0;
                    WeightedSumDenom = 0.0;
                    background = 0.0;
                    noise = 0.0;
                    for (k = 0; k < NumberZones; k++) {
                        WeightedSumBg += ComputeWeightAtXY(Cellxy.x, Cellxy.y, ZonesInfo.pZones.get(k).center.x, ZonesInfo.pZones.get(k).center.y, smoothF) * ZonesInfo.pZones.get(k).background;
                        WeightedSumNoise += ComputeWeightAtXY(Cellxy.x, Cellxy.y, ZonesInfo.pZones.get(k).center.x, ZonesInfo.pZones.get(k).center.y, smoothF) * ZonesInfo.pZones.get(k).noise;
                        WeightedSumDenom += ComputeWeightAtXY(Cellxy.x, Cellxy.y, ZonesInfo.pZones.get(k).center.x, ZonesInfo.pZones.get(k).center.y, smoothF);
                    }
                    if (WeightedSumDenom != 0.0) {
                        background = WeightedSumBg / WeightedSumDenom;
                        noise = WeightedSumNoise / WeightedSumDenom;
                    }

                    BG.get(iUnit).get(iAtom).value2 = background;
                    Noise.get(iUnit).get(iAtom).value2 = noise;

                    inten = pCell.getIntensity((int) Cellxy.x, (int) Cellxy.y);
                    modifiedI = ModifyIntensitySlightly(inten);
                    scaledAdjustedI = ComputeAdjustedIntensity(modifiedI, background, noise);
                    MM.get(iUnit).add(scaledAdjustedI);

                    //atomCount++;
                }
                //PM.get(iUnit).setSize(atomCount);
                //MM.get(iUnit).setSize(atomCount);
            } // if expression type 
        } // for each unit

        // If a cel file with adjusted intensities is requested, loop over all probes.
        if (m_WriteCell) {
            //FeatureIntensity.setSize(m_NumFeatures);
            for (int i = 0; i < m_NumFeatures; i++)
                FeatureIntensity.add(0.0);
            final int colCount = m_Cdf.getHeader().getCols();
            final int rowCount = m_Cdf.getHeader().getRows();
            for (int colIx = 0; colIx < colCount; ++colIx)
                for (int rowIx = 0; rowIx < rowCount; ++rowIx) {
                    double WeightedSumBg = 0.0;
                    double WeightedSumNoise = 0.0;
                    double WeightedSumDenom = 0.0;
                    double background = 0.0;
                    double noise = 0.0;
                    for (int k = 0; k < NumberZones; k++) {
                        WeightedSumBg += ComputeWeightAtXY(colIx, rowIx, ZonesInfo.pZones.get(k).center.x, ZonesInfo.pZones.get(k).center.y, smoothF) * ZonesInfo.pZones.get(k).background;
                        WeightedSumNoise += ComputeWeightAtXY(colIx, rowIx, ZonesInfo.pZones.get(k).center.x, ZonesInfo.pZones.get(k).center.y, smoothF) * ZonesInfo.pZones.get(k).noise;
                        WeightedSumDenom += ComputeWeightAtXY(colIx, rowIx, ZonesInfo.pZones.get(k).center.x, ZonesInfo.pZones.get(k).center.y, smoothF);
                    }
                    if (WeightedSumDenom != 0.0) {
                        background = WeightedSumBg / WeightedSumDenom;
                        noise = WeightedSumNoise / WeightedSumDenom;
                    }

                    double inten = pCell.getIntensity(colIx, rowIx);
                    double modifiedI = ModifyIntensitySlightly(inten);
                    double scaledAdjustedI = ComputeAdjustedIntensity(modifiedI, background, noise);
                    int xyindex = m_Cell.xyToIndex(colIx, rowIx);
                    FeatureIntensity.set(xyindex, scaledAdjustedI);
                }
        }

		/*		//Clean up
		delete [] NumberCellsPerZone;	*/
    }

    //////////////////////////////////////////////////////////////////////

    void ComputeBackGroundZones(FusionCELData pCell, AllZonesInfoType ZonesInfo,
                                List<Double> FeatureIntensity) throws UnsignedOutOfLimitsException, IOException {
        int iUnit;
        int CellsRemaining;
        int zonex;
        int zoney;
        int NumberZones;
        // int iInten=0; // unused pk

        //ArrayList<Integer> NumberCellsPerZone = new ArrayList<Integer> ();

        // Determine the number of remaining cells in the vertical direction.
        CellsRemaining = m_Cdf.getHeader().getCols() % m_Params.NumberVertZones;
        if (CellsRemaining != 0) {
            zonex = (m_Cdf.getHeader().getCols() +
                    (int) m_Params.NumberVertZones - CellsRemaining) /
                    (int) m_Params.NumberVertZones;
        } else {
            zonex = m_Cdf.getHeader().getCols() / m_Params.NumberVertZones;
        }

        // Determine the number of remaining cells in the horizontal direction.
        CellsRemaining = m_Cdf.getHeader().getRows() % m_Params.NumberHorZones;
        if (CellsRemaining != 0) {
            zoney = (m_Cdf.getHeader().getRows() +
                    m_Params.NumberHorZones - CellsRemaining) /
                    m_Params.NumberHorZones;
        } else {
            zoney = m_Cdf.getHeader().getRows() / m_Params.NumberHorZones;
        }

        // Ensure that there are a match and mismatch cell in the same zone.
        zoney += zoney % 2; //EXPRESSION_ATOMS_PER_CELL;

        // Determine the total number of zones.
        NumberZones = (int) m_Params.NumberVertZones * m_Params.NumberHorZones;

        // Get number of units.
        int NumUnits = m_Cdf.getHeader().getNumProbeSets();

        // Allocate space for all atoms intensities and ID's.
        int[] NumberCellsPerZone = new int[NumberZones];

        // Clear arrays.
        //memset(NumberCellsPerZone, 0, sizeof(int)*NumberZones);

        // Loop over all units to determine the zone ID's and intensities.
        List<Vector<Double>> ZoneCells = new ArrayList<Vector<Double>>(NumberZones);
        for (int i = 0; i < NumberZones; i++)
            ZoneCells.add(new Vector<Double>());
        FusionCDFProbeSetInformation unit = new FusionCDFProbeSetInformation();
        FusionCDFProbeGroupInformation blk = new FusionCDFProbeGroupInformation();
        FusionCDFProbeInformation cell = new FusionCDFProbeInformation();
        boolean bMasked;
        for (iUnit = 0; iUnit < NumUnits; ++iUnit) {
            m_Cdf.getProbeSetInformation(iUnit, unit);

            // Only process expression units.
            if (unit.getProbeSetType() == GeneChipProbeSetType.ExpressionProbeSetType) {
                unit.getGroup(0, blk);
                int numCells = blk.getNumCells();

                // Loop over the atoms in the unit
                for (int iCell = 0; iCell < numCells; iCell++) {
                    blk.getCell(iCell, cell);
                    bMasked = IsMasked(pCell, cell.getX(), cell.getY());
                    if (bMasked == false) {
                        int nZone;
                        nZone = DetermineZone(cell.getX(), cell.getY(), zonex, zoney);
                        ZoneCells.get(nZone).setSize(ZoneCells.get(nZone).size() + 1);
                        ZoneCells.get(nZone).add(0.0);
                        ZoneCells.get(nZone).set(NumberCellsPerZone[nZone],
                                (double) pCell.getIntensity(cell.getX(), cell.getY()));

                        if (nZone >= 0 && nZone < NumberZones)
                            NumberCellsPerZone[nZone]++;
                    }
                }
            }
        }

        // Allocate zones, set smooth factor and set num zones
        ZonesInfo.pZones = new ArrayList<ZoneInfo>(NumberZones);
        ZonesInfo.number_zones = NumberZones;
        ZonesInfo.smooth_factor = m_Params.SmoothFactorBG;

        // compute background for each zone
        for (int iZone = 0; iZone < NumberZones; iZone++) {
            // Compute the center coordinates of each zone.
            // (x1,y1) is the upper left corner
            // (x2,y2) is the lower right corner
            double x1 = ((int) (iZone % m_Params.NumberVertZones)) * zonex;
            double y1 = ((int) (iZone / m_Params.NumberVertZones)) * zoney;
            double x2 = x1 + zonex;
            double y2 = y1 + zoney;
            ZonesInfo.pZones.get(iZone).center.x = (x1 + x2) / 2;
            ZonesInfo.pZones.get(iZone).center.y = (y1 + y2) / 2;

            //int iCell=0;
            int numCell = NumberCellsPerZone[iZone];
            ZonesInfo.pZones.get(iZone).numCell = numCell;

            List<Double> zoneI = new ArrayList<Double>(numCell);
            //List<Integer> rank = new ArrayList<Integer>(numCell);

            for (int i = 0; i < numCell; i++) {
                double inten = ZoneCells.get(iZone).get(i);
                zoneI.add(ModifyIntensitySlightly(inten));
                //iCell++;
            }
            double lowBG = 0.0;
            double highBG = m_Params.NumberBGCells / 100.0;
            DoublePair fp = trimMeanAndStd(zoneI, lowBG, highBG);
            ZonesInfo.pZones.get(iZone).background = fp.value1;
            ZonesInfo.pZones.get(iZone).noise = fp.value2;
        }

        //Clean up
        //delete [] NumberCellsPerZone;	
    }


    //////////////////////////////////////////////////////////////////////

    void ComputeMeasurement(List<ArrayList<Double>> PM,
                            List<ArrayList<Double>> MM,
                            List<ArrayList<Boolean>> UseAtom,
                            List<Double> avgMeasurement,
                            List<ArrayList<Double>> PV) {
        // Compute the typical difference (for each probe set) of log(PM) over log(MM).
        int nProbeSet = PM.size(); //  Number of Probe Set

        // Compute Contrast Value
        List<ArrayList<Double>> CT = new ArrayList<ArrayList<Double>>(nProbeSet);
        List<ArrayList<Double>> PVused = new ArrayList<ArrayList<Double>>(nProbeSet);

        for (int i = 0; i < nProbeSet; i++) {
            CT.add(new ArrayList<Double>());
            PVused.add(new ArrayList<Double>());
        }
        ComputeContrastValue(PM, MM, CT);

        // Compute Probe Value
        ComputeProbeValue(PM, CT, PV);

        // Compute Average Log Intensity for probe set i and its confidence interval.
        double c = m_Params.TuningConstantCAvgLogInten;
        double epsilon = m_Params.EpsilonAvgLogInten;


        GetUsedSet(PV, UseAtom, PVused);

        //List<Double> uncertainty = new ArrayList<Double> (nProbeSet);
        for (int i = 0; i < nProbeSet; i++) {
            avgMeasurement.add(OneStepBiweightAlgorithm(PVused.get(i), c, epsilon));
            avgMeasurement.set(i, antiLog(avgMeasurement.get(i)));
        }
    }

    //////////////////////////////////////////////////////////////////////

    void GetUsedSet(List<ArrayList<Double>> InputSet,
                    List<ArrayList<Boolean>> UseAtom,
                    List<ArrayList<Double>> UsedSet) {
        int nProbeSet = InputSet.size();
        for (int i = 0; i < nProbeSet; i++) {
            int size = InputSet.get(i).size();
            //UsedSet.get(i).setSize(size);
            //int used = 0;
            for (int j = 0; j < size; j++) {
                if (UseAtom.get(i).get(j)) {
                    UsedSet.get(i).add(InputSet.get(i).get(j));
                    //used++;
                }
            }
            //UsedSet.get(i).setSize(used);
        }
    }

    //////////////////////////////////////////////////////////////////////

    void ComputeContrastValue(List<ArrayList<Double>> PM,
                              List<ArrayList<Double>> MM,
                              List<ArrayList<Double>> CT) {
        int nProbeSet = PM.size();
        List<Double> SB = new ArrayList<Double>(nProbeSet);
        ComputeTypicalDifference(PM, MM, SB);
        double ContrastTau = m_Params.ContrastTau;
        for (int i = 0; i < nProbeSet; i++) {
            //CT.add(new ArrayList<Double>());
            int nProbePair = PM.get(i).size();
            //CT.get(i).setSize(nProbePair);
            for (int j = 0; j < nProbePair; j++) {
                if (MM.get(i).get(j) < PM.get(i).get(j)) {
                    CT.get(i).add(MM.get(i).get(j));
                } else if ((MM.get(i).get(j) >= PM.get(i).get(j)) &&
                        (SB.get(i) > ContrastTau)) {
                    CT.get(i).add(PM.get(i).get(j) / antiLog(SB.get(i)));
                } else if ((MM.get(i).get(j) >= PM.get(i).get(j)) &&
                        (SB.get(i) <= ContrastTau)) {
                    CT.get(i).add(PM.get(i).get(j) / antiLog(ContrastTau / (1.0 + (ContrastTau - SB.get(i)) / m_Params.ScaleTau)));
                }
            }
        }
    }

    //////////////////////////////////////////////////////////////////////

    void ComputeProbeValue(List<ArrayList<Double>> PM,
                           List<ArrayList<Double>> CT,
                           List<ArrayList<Double>> PV) {
        int nProbeSet = PM.size();
        double delta = m_Params.Delta;
        double correction = 1.0 + m_Params.BiasCorrect;
        for (int i = 0; i < nProbeSet; i++) {
            //PV.add(new ArrayList<Double>());
            int nProbePair = PM.get(i).size();
            //PV.get(i).setSize(nProbePair);
            for (int j = 0; j < nProbePair; j++) {
                double v = PM.get(i).get(j) - CT.get(i).get(j);
                if (v < delta)
                    PV.get(i).add(correction * logtwo(delta));
                else
                    PV.get(i).add(correction * logtwo(v));
            }
        }
    }

    //////////////////////////////////////////////////////////////////////

    void ComputeTypicalDifference(List<ArrayList<Double>> PM,
                                  List<ArrayList<Double>> MM,
                                  List<Double> SB) {
        int nProbeSet = PM.size();
        List<ArrayList<Double>> logPM_minus_logMM = new ArrayList<ArrayList<Double>>(nProbeSet);
        for (int i = 0; i < nProbeSet; i++)
            logPM_minus_logMM.add(new ArrayList<Double>());
        double c = m_Params.TuningConstantCSB;
        double epsilon = m_Params.EpsilonSB;

        for (int i = 0; i < nProbeSet; i++) {
            int nProbePair = PM.get(i).size(); // Number of Probe Pair in Probe Set i
            //logPM_minus_logMM.get(i).setSize(nProbePair);
            for (int j = 0; j < nProbePair; j++) {
                logPM_minus_logMM.get(i).add(logtwo(PM.get(i).get(j)) - logtwo(MM.get(i).get(j)));
            }
            SB.add(OneStepBiweightAlgorithm(logPM_minus_logMM.get(i), c, epsilon));
        }
    }

    //////////////////////////////////////////////////////////////////////

    void ModifyIntensities
            (
                    List<AbsStatExpressionProbeSetResultType> statResults,
                    double factor
            ) {
        // Loop over all of the units.
        int UnitsPerChip = m_Cdf.getHeader().getNumProbeSets();
        for (int iUnit = 0; iUnit < UnitsPerChip; iUnit++) {
            // Get the units.
            AbsStatExpressionProbeSetResultType pUnitResult = statResults.get(iUnit);
            pUnitResult.Signal = (pUnitResult.Signal * factor);
        }
    }

    //////////////////////////////////////////////////////////////////////

    double DetermineScaleFactor(List<AbsStatExpressionProbeSetResultType> statResults) {
        int iUnit;
        double avg = 0.0;


        // User defined norm factor is already stored in the algorithm parameters.
        if (m_Params.SFMethod == ExpStatAlgSettings.ScalingOptionsEnum.DEFINED_SCALING_FACTOR) {
            return m_Params.ScaleFactor;
        }

        // Loop over all of the units.
        int UnitsPerChip = m_Cdf.getHeader().getNumProbeSets();
        int unitCount = 0;
        List<Double> intensityList = new ArrayList<Double>(UnitsPerChip);
        for (iUnit = 0; iUnit < UnitsPerChip; iUnit++) {
            // Get the units.
            AbsStatExpressionProbeSetResultType pUnitResult = statResults.get(iUnit);

            // find signal only for the used genes.
            //String probeSetName = m_Cdf.getProbeSetName(iUnit);
            if (m_Params.SFMethod == ExpStatAlgSettings.ScalingOptionsEnum.SCALE_TO_ALL_PROBE_SETS)// ||
            //UseUnitInNormFactor(probeSetName, m_Params.ScaleGenes))
            {
                // Use Measurement to do scale factor
                if (pUnitResult.NumUsedPairs != 0) {
                    intensityList.add(pUnitResult.Signal);
                    unitCount++;
                }
            }
        }

        //intensityList.setSize(unitCount);

        // Compute the trimMean
        double p1 = m_Params.IntensityLowPercent / 100.0;
        double p2 = 1.0 - m_Params.IntensityHighPercent / 100.0;
        avg = trimMean(intensityList, p1, p2);

        // Store the scale factor
        double sf = 1.0;
        if (unitCount != 0 && avg != 0.0) {
            sf = m_Params.TGT / avg;
        }

        //Check for the validity of SF
        if (sf <= 0) {
            sf = 1.0;
        }
        return sf;
    }

    //////////////////////////////////////////////////////////////////////

    boolean UseUnitInNormFactor(String probeSetName, List<String> scaleGenes) {
        boolean UseUnit = false;
        for (Iterator<String> iter = scaleGenes.iterator(); iter.hasNext(); ) {
            if (probeSetName == iter.next()) {
                UseUnit = true;
                break;
            }
        }
        return UseUnit;
    }

    //////////////////////////////////////////////////////////////////////

    void ComputeAbsoluteCall(AbsStatExpressionProbeSetResultType unit, ExpResults res) {
        // Determine if the gene is present.
        if (res.call == 2) {
            unit.Detection = EXP_PRESENT_CALL_TYPE;
            presentCalls++;
        }

        // Determine if marginal
        else if (res.call == 1)
            unit.Detection = EXP_MARGINAL_CALL_TYPE;

            // Otherwise the gene is absent.
        else
            unit.Detection = EXP_ABSENT_CALL_TYPE;

        allCalls++;
    }

    //////////////////////////////////////////////////////////////////////

    double ModifyIntensitySlightly(double intensity) {
        return Math.max(intensity, m_Params.Epsilon);
    }

    //////////////////////////////////////////////////////////////////////

    double ComputeAdjustedIntensity(double intensity, double background, double noise) {
        double factoredNoise = noise * m_Params.NoiseFrac;
        double diff = intensity - background;
        double adjustedI = Math.max(Math.max(diff, factoredNoise), 0.5);
        // AlexC - 1/22/01
        // Code Comments:
        // if too frequent substitution of the noise value, alert might generated.
        // Refer to the page 4 of the Background and Spatial Variation Adjustement spec., 
        // under eq(16), it said that 
        // "Production software should probably alert the user if the noise value is being 
        //  substituted too frequently (indicating that too much data is below the noise level), 
        //  but an appropriate threshold value is not at present known."
        return adjustedI;
    }

    //////////////////////////////////////////////////////////////////////

    int DetermineZone(int cellx, int celly, int zonex, int zoney) {
        double fZx = (double) cellx / (double) zonex;
        double fZy = (double) celly / (double) zoney;

        int Zx = (int) Math.floor(fZx);
        int Zy = (int) Math.floor(fZy);

        int zoneID = Zx + Zy * m_Params.NumberVertZones;
        return zoneID;
    }

    //////////////////////////////////////////////////////////////////////

    void SetMeasurement
            (
                    List<ArrayList<DoublePair>> BG,
                    List<Double> avgMeasurement,
                    List<AbsStatExpressionProbeSetResultType> statResults
            ) {
        // FusionCDFProbeSetInformation *unit=NULL; // unused pk
        // FusionCDFProbeGroupInformation *blk=NULL; // unused pk
        // FusionCDFProbeInformation *cell=NULL; // unused pk

        int UnitsPerChip = m_Cdf.getHeader().getNumProbeSets();
        for (int iUnit = 0; iUnit < UnitsPerChip; iUnit++) {
            AbsStatExpressionProbeSetResultType pUnitResult = statResults.get(iUnit);
            if (pUnitResult.Detection == EXP_NO_ABS_CALL_TYPE)
                pUnitResult.Signal = 0.0;
            else
                pUnitResult.Signal = avgMeasurement.get(iUnit);
        }
    }

    //////////////////////////////////////////////////////////////////////

    double ComputeCellBackground
            (
                    int x,
                    int y,
                    List<ZoneInfo> zones,
                    int zoneSize,
                    double smoothFactorBg
            ) {
        double bkg = 0.0;
        if (zones != null) {
            double weightedSumBg = 0.0;
            double weight = 0.0;
            double weightedSumDenom = 0.0;
            double background = 0.0;

            for (int k = 0; k < zoneSize; k++) {
                weight = ComputeWeightAtXY(x, y, zones.get(k).center.x, zones.get(k).center.y, smoothFactorBg);
                weightedSumBg += weight * zones.get(k).background;
                weightedSumDenom += weight;
            }
            if (weightedSumDenom != 0.0) {
                background = weightedSumBg / weightedSumDenom;
            }
            bkg = background;
        }

        return bkg;
    }

    //////////////////////////////////////////////////////////////////////

    void ComputeAvgMaxMinStdBgNoise(List<ArrayList<DoublePair>> BG, List<ArrayList<DoublePair>> Noise) {
        double maxBg = 0.0;
        double minBg = 0.0;
        double avgBg = 0.0;
        double stddevBg = 0.0;
        double maxNoise = 0.0;
        double minNoise = 0.0;
        double avgNoise = 0.0;
        double stddevNoise = 0.0;

        List<Double> bgList = new ArrayList<Double>();
        List<Double> noiseList = new ArrayList<Double>();
        //int count=0;

        int UnitsPerChip = m_Cdf.getHeader().getNumProbeSets();
        for (int iUnit = 0; iUnit < UnitsPerChip; iUnit++) {
            int nAtoms = BG.get(iUnit).size();
            //bgList.setSize(count + nAtoms * 2);
            //noiseList.setSize(count + nAtoms * 2);
            for (int iAtom = 0; iAtom < nAtoms; iAtom++) {
                bgList.add(BG.get(iUnit).get(iAtom).value1);
                noiseList.add(Noise.get(iUnit).get(iAtom).value1);
                //count++;

                bgList.add(BG.get(iUnit).get(iAtom).value2);
                noiseList.add(Noise.get(iUnit).get(iAtom).value2);
                //count++;
            }
        }

        // Compute the min/max/avg/stdev values

        //		maxBg = *max_element(bgList.begin(), bgList.end());
        //		minBg = *min_element(bgList.begin(), bgList.end());
        maxBg = Collections.max(bgList);
        minBg = Collections.min(bgList);
        //		maxNoise = *max_element(noiseList.begin(), noiseList.end());
        //		minNoise = *min_element(noiseList.begin(), noiseList.end());
        maxNoise = Collections.max(noiseList);
        minNoise = Collections.min(noiseList);
        avgBg = mean(bgList);
        avgNoise = mean(noiseList);
        stddevBg = stddev(bgList);
        stddevNoise = stddev(noiseList);


        // Store the noise stats
        m_NoiseStats.avg = avgNoise;
        m_NoiseStats.stdv = stddevNoise;
        m_NoiseStats.min = minNoise;
        m_NoiseStats.max = maxNoise;


        // Store the background stats
        m_BgStats.avg = avgBg;
        m_BgStats.stdv = stddevBg;
        m_BgStats.min = minBg;
        m_BgStats.max = maxBg;
    }

    //////////////////////////////////////////////////////////////////////

    void ReportCornerControls(int qcType) throws UnsignedOutOfLimitsException, IOException {
        // Get the desired qc probe set.
        FusionCDFQCProbeSetInformation qcUnit = new FusionCDFQCProbeSetInformation();
        m_Cdf.getQCProbeSetInformationByType(qcType, qcUnit);
        if (qcUnit.getNumCells() == 0)
            return;

        // Get the average noise.
        final double NOISE_FRAC = 0.5;
        double avgNoise = m_NoiseStats.avg * NOISE_FRAC;

        // Compute the average intensity
        int size = qcUnit.getNumCells();
        double avgIntensity = 0.0;
        FusionCDFQCProbeInformation qcCell = new FusionCDFQCProbeInformation();
        for (int iCell = 0; iCell < size; iCell++) {
            qcUnit.getCell(iCell, qcCell);
            double fQCBg = ComputeCellBackground(qcCell.getX(), qcCell.getY(),
                    m_ZonesInfo.pZones, m_ZonesInfo.number_zones, m_Params.SmoothFactorBG);
            double fQCIntensity = m_Cell.getIntensity(qcCell.getX(), qcCell.getY());
            double fDelta = fQCIntensity - fQCBg;
            if (fDelta <= avgNoise)
                avgIntensity += avgNoise;
            else
                avgIntensity += fDelta;
        }

        // Compute average intensity.
        if (size != 0) {
            avgIntensity /= size;
            ControlInformationType info = new ControlInformationType();
            info.avg = avgIntensity;
            info.count = size;
            info.qcType = qcType;
            m_ControlInfo.add(info);
        }

    }

    void getTime(String location) {
        Date date = new Date();

        // display time and date using toString()
        System.out.println(location + ": " + date.toString());
    }
}


////////////////////////////////////////////////////////////////
class ExpResults {
    double p_value;
    int call;

    public ExpResults(double p, int c) {
        p_value = p;
        call = c;
    }

    public ExpResults() {
        p_value = 0.0;
        call = 0;
    }
}

class ExpResultsComparator implements Comparator<ExpResults> {

    public int compare(ExpResults o1, ExpResults o2) {
        return o1.p_value > o2.p_value ? 1 : -1;
    }
}
//////////////////////////////////////////////////////////////////////

class AbsStatExpressionProbeSetResultType {
    double DetectionPValue;
    double Signal;
    int NumPairs;
    int NumUsedPairs;
    int Detection;
	/*	_AbsStatExpressionProbeSetResultType *operator=(_AbsStatExpressionProbeSetResultType &src)
{
memcpy(this, &src, sizeof(_AbsStatExpressionProbeSetResultType));
return this;
};*/
}


////////////////////////////////////////////////////////////////////

class CompStatExpressionProbeSetResultType {
    double ChangePValue;
    double SignalLogRatio;
    double SignalLogRatioLow;
    double SignalLogRatioHigh;
    int NumCommonPairs;
    int Change;
	/*	_CompStatExpressionProbeSetResultType *operator=(_CompStatExpressionProbeSetResultType &src)
{
memcpy(this, &src, sizeof(_CompStatExpressionProbeSetResultType));
return this;
};*/
}

////////////////////////////////////////////////////////////////////

class Coordinate {
    double x;
    double y;
}

////////////////////////////////////////////////////////////////////

class ZoneInfo {
    Coordinate center;
    int numCell;
    double background;
    double noise;

    public ZoneInfo() {
        center = new Coordinate();
    }
}

////////////////////////////////////////////////////////////////////

class AllZonesInfoType {
    int number_zones;
    double smooth_factor;
    List<ZoneInfo> pZones;
}

////////////////////////////////////////////////////////////////////

class AvgStdvMinMaxType {
    double avg;
    double stdv;
    double min;
    double max;
}

////////////////////////////////////////////////////////////////////

class ControlInformationType {
    int qcType;
    double avg;
    int count;
}



