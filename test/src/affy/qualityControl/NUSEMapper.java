package affy.qualityControl;

import org.apache.commons.lang3.StringUtils;
import org.apache.hadoop.io.DoubleWritable;
import org.apache.hadoop.io.Text;
import org.apache.hadoop.mapreduce.Mapper;

import java.io.IOException;

/**
 * Created with IntelliJ IDEA.
 * User: David
 * Date: 6/18/13
 * Time: 11:26 PM
 * To change this template use File | Settings | File Templates.
 */
public class NUSEMapper extends
        Mapper<Text, Text, Text, DoubleWritable> {
    private Text outkey = new Text();
    private DoubleWritable outvalue = new DoubleWritable();

    @Override
    protected void map(Text key, Text value, Context context) throws IOException, InterruptedException {
        String [] strings = StringUtils.split(value.toString().trim(), ',');
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
}
