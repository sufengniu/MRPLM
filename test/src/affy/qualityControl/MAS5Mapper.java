package affy.qualityControl;

import affy.mas5.ExpressionAlgorithm;
import affymetrix.calvin.exception.UnsignedOutOfLimitsException;
import org.apache.hadoop.io.Text;
import org.apache.hadoop.mapreduce.Mapper;

import java.io.File;
import java.io.FilenameFilter;
import java.io.IOException;

/**
 * Created with IntelliJ IDEA.
 * User: David
 * Date: 6/4/13
 * Time: 12:18 PM
 * To change this template use File | Settings | File Templates.
 */
public class MAS5Mapper extends
        Mapper<Text, Text, Text, Text> {
    String celPath = null;
    String cdfPath = null;
    String cdfFile = null;

    private Text outkey = new Text();
    private Text outvalue = new Text();

    public void setup(Context context) throws IOException,
            InterruptedException {
        celPath = context.getConfiguration().get("cel_path");
        cdfPath = context.getConfiguration().get("cdf_path");
        cdfFile = cdfPath + getCDFfile(cdfPath);
    }

    public void map(Text key, Text value, Context context) throws IOException, InterruptedException {
        try {
            File local = new File("/local_scratch/" + key);
            ReadFileUtils.copyFile(new File(celPath + key), local);
            ExpressionAlgorithm exp = new ExpressionAlgorithm();
            exp.RunStat(local.toString(), "", cdfFile);

            outvalue.set(exp.print());
            context.write(key, outvalue);
            //write the low quality cel file name to HDFS
//            if (exp.isLowQuality().equals("sfs")){
//                outvalue.set("sfs");
//                context.write(key, outvalue);
//            }
//            if (exp.isLowQuality().equals("avbg")){
//                outvalue.set("avbg");
//                context.write(key, outvalue);
//            }
//            if (exp.isLowQuality().equals("pps")){
//                outvalue.set("pps");
//                context.write(key, outvalue);
//            }

            //delete local file
            ReadFileUtils.deleteFileOrDirectory(local);
        } catch (UnsignedOutOfLimitsException e) {
            e.printStackTrace();  //To change body of catch statement use File | Settings | File Templates.
        }
    }

    private String getCDFfile(String dir) {
        File directory = new File(dir);
        if (!directory.isDirectory()) {
            System.err.println("Error! Please select right a directory that contain the .cdf(.CDF) file, instead of " + dir);
            System.exit(-1);
        }

        FilenameFilter filter = new FilenameFilter() {
            public boolean accept(File directory, String fileName) {
                return fileName.endsWith(".cdf") || fileName.endsWith(".CDF");
            }
        };
        String[] list = directory.list(filter);
        if (list.length > 1) {
            System.err.println("Error! The directory have more than one .cdf(.CDF) files!");
            System.exit(-1);
        } else if (list.length == 0) {
            System.err.println("Error! No .cdf(.CDF) file in the directory: " + dir);
            System.exit(-1);
        }

        return list[0];
    }
}
