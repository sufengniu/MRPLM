package affy.qualityControl;

import org.apache.hadoop.io.DoubleWritable;
import org.apache.hadoop.io.Text;
import org.apache.hadoop.mapreduce.Reducer;

import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

/**
 * Created with IntelliJ IDEA.
 * User: David
 * Date: 6/4/13
 * Time: 11:22 PM
 * To change this template use File | Settings | File Templates.
 */
public class NUSE_RLEReducer extends
        Reducer<Text, DoubleWritable, Text, Text> {
    private Text outkey = new Text();
    private Text outvalue = new Text();

    @Override
    protected void reduce(Text key, Iterable<DoubleWritable> values, Context context)
            throws IOException, InterruptedException {
        String[] tokens = key.toString().split(":");
        String flag = tokens[0];

        if (flag.equals("nuse")) {
             nuseReduce(key, values, context);
        }
        if (flag.equals("rle")) {
             rleReduce(key, values, context);
        }
    }

    private void nuseReduce(Text key, Iterable<DoubleWritable> values,
                            Context context) throws IOException, InterruptedException {
        String[] tokens = key.toString().split(":");
        String celName = tokens[1];
        int numNA = Integer.parseInt(tokens[2]);

        List<Double> data = new ArrayList<Double>();
        for (DoubleWritable value : values) {
            data.add(value.get());
        }

        double[] buffer = new double[data.size()];
        for (int i = 0; i < data.size(); i++) {
            if (Double.isNaN(data.get(i))) {
                buffer[i] = Double.POSITIVE_INFINITY;
            } else {
                buffer[i] = data.get(i);
            }
        }
        double med = MedianSteps.median_nocopy_hasNA(buffer, buffer.length, numNA);
        //Q[0]: LQ, Q[1]: UQ
        double[] Q = new double[2];
        MedianSteps.quartiles(buffer, buffer.length - numNA, Q);

//    		StringBuilder sBuilder = new StringBuilder();
//            sBuilder.append(buffer[0]+" ");
//            sBuilder.append(Q[0]+" ");
//            sBuilder.append(med+" ");
//            sBuilder.append(Q[1]+" ");
//            sBuilder.append(buffer[buffer.length-numNA-1]+" ");
        outkey.set(celName);
        //outvalue.set(sBuilder.toString());
        //context.write(outkey, outvalue);
//        if (med > 1.025) {
//            outvalue.set("nuse");
//            context.write(outkey, outvalue);
//        }
        outvalue.set(med + "\tnuse");
        context.write(outkey,outvalue);
    }

    private void rleReduce(Text key, Iterable<DoubleWritable> values,
                           Context context) throws IOException, InterruptedException {
        String[] tokens = key.toString().split(":");
        String celName = tokens[1];
        List<Double> data = new ArrayList<Double>();
        for (DoubleWritable value : values) {
            data.add(value.get());
        }
        double[] buffer = new double[data.size()];

        for (int i = 0; i < data.size(); i++) {
            buffer[i] = data.get(i);
        }

        double med = MedianSteps.median_nocopy(buffer, buffer.length);
        //Q[0]: LQ, Q[1]: UQ
        double[] Q = new double[2];
        MedianSteps.quartiles(buffer, buffer.length, Q);

//		StringBuilder sBuilder = new StringBuilder();
//        sBuilder.append(buffer[0]+" ");
//		sBuilder.append(Q[0]+" ");
//		sBuilder.append(med+" ");
//		sBuilder.append(Q[1]+" ");
//		sBuilder.append(buffer[buffer.length-1]+" ");

        outkey.set(celName);
        //outvalue.set(sBuilder.toString());
        //context.write(outkey, outvalue);

//        if ((med > 0.15) || (med < -0.15)){
//            outvalue.set("rle");
//            context.write(outkey, outvalue);
//        }

        outvalue.set(med + "\trle");
        context.write(outkey,outvalue);
    }
}
