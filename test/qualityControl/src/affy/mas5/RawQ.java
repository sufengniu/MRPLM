package affy.mas5;

import affymetrix.calvin.exception.UnsignedOutOfLimitsException;
import affymetrix.fusion.cdf.*;
import affymetrix.fusion.cel.FusionCELData;
import affymetrix.gcos.cdf.GeneChipProbeSetType;

import java.io.IOException;

public class RawQ {
	private final int EXPRESSION_ATOMS_PER_CELL = 2;
	
	/*! The number of vertical zones for the background calculation. */
	private int m_VertZones = 4;

	/*! The number of horizontal zones for the background calculation. */
	private int m_HorzZones = 4;

	/*! The percentage of probes in a zone to use for the background calculation. */
	private double m_PercentBGCells = 2;
	
	private double LARGE_double_NUMBER = 999999999999999999.0;

	/*! Sets the number of vertical zones for the background calculation.
	 * @param v The number of vertical zones for the background calculation.
	 */
	void setVerticalZones(int v) { 
		m_VertZones = v; 
	}

	/*! Gets the number of vertical zones for the background calculation.
	 * @return The number of vertical zones for the background calculation.
	 */
	int getVerticalZones() 
	{ 
		return m_VertZones; 
	}

	/*! Sets the number of horizontal zones for the background calculation.
	 * @param h The number of horizontal zones for the background calculation.
	 */
	void setHorizontalZones(int h) { m_HorzZones = h; }

	/*! Gets the number of horizontal zones for the background calculation.
	 * @return The number of horizontal zones for the background calculation.
	 */
	int getHorizontalZones() 
	{ 
		return m_HorzZones; 
	}

	/*! Sets the percentage of probes in a zone to use for the background calculation.
	 * @param p The percentage of probes in a zone to use for the background calculation.
	 */
	void setPercentBG(double p) { 
		m_PercentBGCells = p; 
	}

	/*! Gets the percentage of probes in a zone to use for the background calculation.
	 * @return The percentage of probes in a zone to use for the background calculation.
	 */
	double getPercentBG() 
	{ 
		return m_PercentBGCells; 
	}
	
	int DetermineZone(int x, int y, int zonex, int zoney, int vZones)
	{
		double fZx = (double) x / (double) zonex;
		double fZy = (double) y / (double) zoney;
		int Zx = (int)Math.floor(fZx);
		int Zy = (int)Math.floor(fZy);
		int zoneID = Zx + Zy * vZones;
		return zoneID;
	}

	/////////////////////////////////////////////////////////////////////////////

	double computeRawQ(FusionCELData cell, FusionCDFData cdf) throws UnsignedOutOfLimitsException, IOException
	{
	
		// Store the number of background cells.
		FusionCDFHeader cdfHeader = cdf.getHeader();
		double numer = cdfHeader.getCols() * cdfHeader.getRows() * m_PercentBGCells;
		double denom = m_VertZones * m_HorzZones * 100.0;
		int bgCells = (int) (numer / denom);


		// Determine the number of remaining cells in the vertical direction.
		int CellsRemaining = cdfHeader.getCols() % m_VertZones;
		int zonex;
		int zoney;
		if(CellsRemaining != 0)
			zonex = (cdfHeader.getCols() + 
					m_VertZones - CellsRemaining) /
					m_VertZones;
		else
			zonex = cdfHeader.getCols() / m_VertZones;

		// Determine the number of remaining cells in the horizontal direction.
		CellsRemaining = cdfHeader.getRows() % m_HorzZones;
		if(CellsRemaining != 0)
			zoney = (cdfHeader.getRows() +
					m_HorzZones-CellsRemaining) /
					m_HorzZones;
		else
			zoney = cdfHeader.getRows() / m_HorzZones;

		// Ensure that there are a match and mismatch cell in the same zone.
		zoney += zoney % EXPRESSION_ATOMS_PER_CELL;


		// Determine the total number of zones.
		int NumberZones = m_VertZones * m_HorzZones;

		// Allocate memory for storing background data.
		double []zonebg = new double[NumberZones];
		double []zonenoise = new double[NumberZones];
		CellStatisticsType []bgN = new CellStatisticsType[NumberZones*bgCells];
		for (int i=0; i<NumberZones*bgCells; i++)
			bgN[i] = new CellStatisticsType();
		// Get number of units.
		int NumUnits = cdfHeader.getNumProbeSets();

		// Determine the total number of atoms in the chip.
		FusionCDFProbeSetInformation unit = new FusionCDFProbeSetInformation();
		int totalCells=0;
		for (int iUnit=0; iUnit<NumUnits; ++iUnit)
		{
			if (cdf.getProbeSetType(iUnit) == GeneChipProbeSetType.ExpressionProbeSetType)
			{
				cdf.getProbeSetInformation(iUnit, unit);
				totalCells += unit.getNumCells();
			}
		}

		// Allocate space for all atoms intensities and ID's.
		int [] entryIndex = new int[totalCells];
		int [] intenszid = new int[totalCells];

		// Loop over all units to determine the zone ID's and intensities.
		int iInten=0;
		for (int iUnit=0; iUnit<NumUnits; ++iUnit)
		{
			// Only process expression units.
			if (cdf.getProbeSetType(iUnit) != GeneChipProbeSetType.ExpressionProbeSetType)
				continue;

			// Get the PM and MM intensity objects and their zone ids.
			cdf.getProbeSetInformation(iUnit, unit);
			FusionCDFProbeGroupInformation group = new FusionCDFProbeGroupInformation();
			int nGroups = unit.getNumGroups();
			for (int iGroup=0; iGroup<nGroups; iGroup++)
			{
				unit.getGroup(iGroup, group);
				int nCells = group.getNumCells();
				FusionCDFProbeInformation probe = new FusionCDFProbeInformation();
				for (int iCell=0; iCell<nCells; iCell++)
				{
					group.getCell(iCell, probe);
					entryIndex[iInten] = cell.xyToIndex(probe.getX(), probe.getY());
					intenszid[iInten] = DetermineZone(probe.getX(), probe.getY(), zonex, zoney, m_VertZones);
					++iInten;
				}
			}
		}

		// compute background for each zone
		for (int iZone=0; iZone<NumberZones; iZone++)
		{
			// Initialize the background.
			for(int bgcnt = 0; bgcnt < bgCells; bgcnt++)
				bgN[bgcnt+(iZone*bgCells)].intensity = LARGE_double_NUMBER;

			// find the lowest N intensities in each zone
			for (iInten = 0; iInten < totalCells; iInten++)
			{
				// Only process those intensities in the current zone.
				if(intenszid[iInten] == iZone)
				{
					int index_cnt;
					int index_max;
					for (index_cnt=1, index_max=0;
						 index_cnt < bgCells; index_cnt++)
					{
						if(bgN[index_cnt+(iZone*bgCells)].intensity > bgN[index_max+(iZone*bgCells)].intensity)
							index_max = index_cnt;
					}


					// Store the low intensity.
					double intensity = Math.min(bgN[index_max+(iZone*bgCells)].intensity, cell.getIntensity(entryIndex[iInten]));
					if (intensity != bgN[index_max+(iZone*bgCells)].intensity)
					{
						bgN[index_max+(iZone*bgCells)].intensity = cell.getIntensity(entryIndex[iInten]);
						bgN[index_max+(iZone*bgCells)].pixel = cell.getPixels(entryIndex[iInten]);
						bgN[index_max+(iZone*bgCells)].stdev = cell.getStdv(entryIndex[iInten]);
					}
				}
			}

			// compute the average
			double bgSum = 0.0;
			for(int bgcnt = 0; bgcnt < bgCells; bgcnt++)
				bgSum += bgN[bgcnt+(iZone*bgCells)].intensity;
			zonebg[iZone] = bgSum / bgCells;

			// Compute the noise.
			bgSum = 0.0;
			for(int bgcnt = 0; bgcnt < bgCells; bgcnt++)
				bgSum += bgN[bgcnt+(iZone*bgCells)].stdev / (double)(Math.sqrt((double)bgN[bgcnt+(iZone*bgCells)].pixel));
			zonenoise[iZone] = bgSum / bgCells;
		}

		// Compute the average noise.
		double avgNoise=0.0;
		for (int iZone=0; iZone<NumberZones; iZone++)
			avgNoise+=zonenoise[iZone];
		double rawQ = avgNoise/NumberZones;


		// Return the RawQ value
	 	return rawQ;
	}
}

class CellStatisticsType {
	public double intensity;
	public double stdev;
	public int   pixel;
}
