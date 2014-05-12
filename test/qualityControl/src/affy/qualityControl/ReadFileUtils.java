package affy.qualityControl;

import affymetrix.fusion.cdf.FusionCDFData;
import affymetrix.fusion.cel.FusionCELData;

import java.io.File;
import java.io.FilenameFilter;
import java.io.IOException;
import java.nio.file.*;

/**
 * Created with IntelliJ IDEA.
 * User: David
 * Date: 6/16/13
 * Time: 4:05 PM
 * To change this template use File | Settings | File Templates.
 */
public class ReadFileUtils {
    public static void copyFile(File src, File dest) throws IOException {
        Path FROM = Paths.get(src.toURI());
        Path TO = Paths.get(dest.toURI());
        //overwrite existing file, if exists
        CopyOption[] options = new CopyOption[]{
                StandardCopyOption.REPLACE_EXISTING,
                StandardCopyOption.COPY_ATTRIBUTES
        };
        Files.copy(FROM, TO, options);
    }

    public static void deleteFileOrDirectory(File file) {
        try {
            Files.delete(Paths.get(file.toURI()));
        } catch (NoSuchFileException x) {
            System.err.format("%s: no such file or directory%n", file.getName());
        } catch (DirectoryNotEmptyException x) {
            System.err.format("%s not empty%n", file.getName());
        } catch (IOException x) {
            // File permission problems are caught here.
            System.err.println(x);
        }
    }

    public static FusionCDFData getCdfData(String filename){
        FusionCDFData cdfData = new FusionCDFData();
        cdfData.setFileName(filename);
        if (cdfData.exists() == false) {
            System.err.println("The CDF file: "+filename+" does not exist.");
            System.exit(-1);
        }
        if (cdfData.read() == false) {
            System.err.println("Failed to read the CDF file: "+filename);
            System.exit(-1);
        }
        return cdfData;
    }

    public static String getCDFfile(String dir) {
        File directory = new File(dir);
        if (!directory.isDirectory()) {
            System.err.println("Error! Please select right a directory that contain the .cdf(.CDF) file, instead of "+dir);
            System.exit(-1);
        }

        FilenameFilter filter = new FilenameFilter() {
            public boolean accept(File directory, String fileName) {
                return fileName.endsWith(".cdf") || fileName.endsWith(".CDF");
            }
        };
        String [] list = directory.list(filter);
        if (list.length > 1){
            System.err.println("Error! The directory have more than one .cdf(.CDF) files!");
            System.exit(-1);
        }
        else if (list.length == 0){
            System.err.println("Error! No .cdf(.CDF) file in the directory: "+dir);
            System.exit(-1);
        }
        return list[0];
    }

    public static FusionCELData getCelData(String filename){
        FusionCELData celData = new FusionCELData();
        celData.setFileName(filename);
        if (celData.exists() == false) {
            System.err.println("The CEL file: "+filename+" does not exist.");
            System.exit(-1);
        }
        if (celData.read() == false) {
            System.err.println("Failed to read the CEL file: "+filename);
            System.exit(-1);
        }
        return celData;
    }

    public static String getPath(String arg) {
        String dataPath = arg;
        // make sure the path is end with '/'
        if (dataPath.toCharArray()[dataPath.length()-1] != '/') {
            dataPath += '/';
        }
        return dataPath;
    }
}
