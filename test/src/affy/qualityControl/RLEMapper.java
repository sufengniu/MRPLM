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
 * Time: 11:29 PM
 * To change this template use File | Settings | File Templates.
 */
public class RLEMapper extends
        Mapper<Text, Text, Text, DoubleWritable> {
    private Text outkey = new Text();
    private DoubleWritable outvalue = new DoubleWritable();
    @Override
    protected void map(Text key, Text value, Mapper.Context context) throws IOException, InterruptedException {
        String [] strings = StringUtils.split(value.toString().trim(), ',');
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
