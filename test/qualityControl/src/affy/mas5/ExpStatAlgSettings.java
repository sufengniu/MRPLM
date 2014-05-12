package affy.mas5;

import java.util.ArrayList;

public class ExpStatAlgSettings {

	enum NormalizationOptionsEnum
	{
		NORM_TO_ALL_PROBE_SETS,
		NORM_TO_SELECTED_PROBE_SETS,
		DEFINED_NORMALIZATION_FACTOR
	};

	enum ScalingOptionsEnum
	{
		SCALE_TO_ALL_PROBE_SETS,
		SCALE_TO_SELECTED_PROBE_SETS,
		DEFINED_SCALING_FACTOR
	};

	NormalizationOptionsEnum NormMethod = NormalizationOptionsEnum.NORM_TO_ALL_PROBE_SETS;
	ScalingOptionsEnum SFMethod = ScalingOptionsEnum.SCALE_TO_ALL_PROBE_SETS;

	String ProbeMaskFile;
	String ScaleMaskFile;
	String NormMaskFile;

	ArrayList<String> ScaleGenes;
	ArrayList<String> NormGenes;

	int NumberHorZones = 4;
	int NumberVertZones = 4;
	double NumberBGCells = 2.0;
	double NormFactor = 1.0;
	double ScaleFactor = 1.0;
	double TGT = 100; //500;
	double IntensityLowPercent = 2.0;
	double IntensityHighPercent = 2.0;

	/* For Absolute Call */

	// User Parameters
	double	Alpha1 = 0.05;
	double	Alpha2 = 0.065;
	double	Tau = 0.015;

	// Other algorithm parameters
	double	HPSaturatedIntensity = 48000.0;
	double	SaturatedIntensity = 65000.0;
	double	Epsilon = 0.5;
	double	ContrastTau = 0.03;
	double	TuningConstantCSB = 5.0;
	double	TuningConstantCAvgLogInten = 5.0;
	double	TuningConstantCAvgLogRatio = 5.0;
	double	TuningConstantCGammas = 5.0;
	double	EpsilonSB = 0.0001;
	double	EpsilonAvgLogInten = 0.0001;
	double	EpsilonAvgLogRatio = 0.0001;
	double	EpsilonGammas = 0.0001;
	double	SmoothFactorBG = 100.0;
	double	NoiseFrac = 0.5;
	double	ScaleTau = 10.0;
	double  log2DELTA = -20;
	double	Delta = Math.pow(2.0, (double) log2DELTA);

	/* For Comparative Call */

	// User Parameters
	double	Gamma1H = 0.004500;
	double	Gamma1L = 0.004500;
	double	Gamma2H = 0.006000;
	double	Gamma2L = 0.006000;
	double	Perturbation = 1.10;

	// Other algorithm parameters
	double	CMultiplier = 0.2;
	double	BHCoef = 7.0;
	double	BLCoef = 0.8;
	double	BiasCorrect = 0.0;
	double	RelConfInterval = 0.975;
	double	STP = 3.0;

	double	BaseScaleFactor;
}
