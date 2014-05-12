package affy.qualityControl;

import org.apache.commons.lang3.StringUtils;
import org.apache.hadoop.io.DoubleWritable;
import org.apache.hadoop.io.Text;
import org.apache.hadoop.mapreduce.Mapper;

import java.io.IOException;

/**
 * Created with IntelliJ IDEA.
 * User: David
 * Date: 6/4/13
 * Time: 11:06 PM
 * To change this template use File | Settings | File Templates.
 */
public class NUSE_RLEMapper extends
        Mapper<Text, Text, Text, DoubleWritable> {
    private Text outkey = new Text();
    private DoubleWritable outvalue = new DoubleWritable();

    // key: probeset_name-[ex|se], value: cell_name:expression or cell_name:SEexpression
    protected void map(Text key, Text value, Context context)
            throws IOException, InterruptedException {
        String [] strings = key.toString().split(":");
        String flag = strings[1];

        if (flag.equals("se")) {
            nuseMapper(key, value, context);
        }
        if (flag.equals("ex")) {
            rleMapper(key, value, context);
        }
    }

    private void nuseMapper(Text key, Text value, Context context) throws IOException, InterruptedException {
        String [] strings = StringUtils.split(value.toString(), ',');
        int numArray = strings.length;
        String [] celNames = new String[numArray];
        double [] expressions = new double[numArray];

        for (int i = 0; i < numArray; i++) {
            String[] tmpStrings = strings[i].split(":");
            celNames[i] = tmpStrings[0];
            expressions[i] = Double.parseDouble(tmpStrings[1]);
        }

        double med = MedianSteps.median(expressions, numArray);
        int numNA = 0;
        if (med == 0.0) {
            med = 1.0;
        }
        if (med != -1.0) {
            for (int j = 0; j < numArray; j++) {
                expressions[j] = expressions[j] / med;
            }
        } else {
            for (int j = 0; j < numArray; j++) {
                expressions[j] = Double.NaN;
            }
            numNA++;
        }
        // output
        for (int i = 0; i < numArray; i++) {
            outkey.set("nuse:"+celNames[i] + ":" + numNA);
            outvalue.set(expressions[i]);
            context.write(outkey, outvalue);
        }
    }

    private void rleMapper(Text key, Text value, Context context) throws IOException, InterruptedException {
        String [] strings = StringUtils.split(value.toString(), ',');
        int numArray = strings.length;
        String [] celNames = new String[numArray];
        double [] expressions = new double[numArray];

        for (int i=0; i<numArray; i++){
            String [] tmpStrings = strings[i].split(":");
            celNames[i] = tmpStrings[0];
            expressions[i] = Double.parseDouble(tmpStrings[1]);
        }

        double med = MedianSteps.median(expressions, numArray);

        for (int i=0; i<numArray; i++){
            expressions[i] = expressions[i] - med;
        }
        // output
        for (int i=0; i<numArray; i++){
            outkey.set("rle:"+celNames[i]);
            outvalue.set(expressions[i]);
            context.write(outkey, outvalue);
        }
    }
}

