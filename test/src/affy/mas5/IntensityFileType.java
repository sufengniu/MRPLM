package affy.mas5;

public class IntensityFileType {
	static boolean IsDatheaderStringFromHP(final String datheader)
	{
	    int index = (int) datheader.indexOf(":CLS=");
	    String strTypeAndID = datheader.substring(index+1);
	    strTypeAndID = strTypeAndID.substring(67);

	    // There are several sub-fields in this field. The 
	    // first sub field is the scanner ID, followed by the scanner type, 
	    // followed by three spaces. If the scanner ID is absent, 
	    // the field consists of four spaces
	    // Find the 3 spaces that are after the scanner type
	    int iStartScanType = (int) strTypeAndID.indexOf(0x14);

	    // If the start of the 3 spaces is at the begining, there is 
	    // no scannertype and scanner ID
	    String strScannerType;
	    if (iStartScanType > 0)
	    {
	        strScannerType = strTypeAndID.substring(0, iStartScanType - 3);
	        iStartScanType = (int) strScannerType.lastIndexOf(" ");
	        if(iStartScanType <= 0)
	            return true;
	        
	        strScannerType = strScannerType.substring(iStartScanType + 1, strScannerType.length() - iStartScanType - 1);
	        return (strScannerType == "HP");
	    }
	    return true;
	}

	/*
	 * Return a flag indicating if the CEL file is from an HP scanner.
	 */
/*	static boolean FromHP(FusionCELData cel)
	{ 
	    // Check if from the HT scanner. This had a upper left grid of 1,1
	    FGridCoords grid = cel.getGridCorners();
	    if (grid.getUpperLeft().getX() < 1.001f && grid.getUpperLeft().getY() < 1.001f)
	        return false;

	    GenericData gdata = cel.getGenericData();
	    if (gdata != null)
	    {
	        // Check for the scanne type parameter.
	        int nparents = gdata.getHeader().getGenericDataHdr().getParentCnt();
	        for (int iparent=0; iparent<nparents; ++iparent)
	        {
	            GenericDataheader pheader = gdata.getHeader().getGenericDataHdr().GetParent(iparent);
	            ParameterNameValueType nvt;
	            if (pheader.FindNameValParam(SCANNER_TYPE_PARAM_NAME, nvt) == true)
	            {
	                String val = nvt.ToString();
	                if (val.size() == 0 || val == "HP")
	                    return true;
	            }
	            else if (pheader.FindNameValParam(DAT_header_PARAM_NAME, nvt) == true)
	            {
	                String val = nvt.ToString();
	                return IsDatheaderStringFromHP(val);
	            }
	        }
	        return false;
	    }
	    else
	    {
	        // Check the DAT header. It it contains HP or blank for the scanner type
	        // then it came from the HP scanner.
	        String datheader = cel.getDatHeader();
	        return IsDatheaderStringFromHP(datheader);
	    }
	}
*/
}
