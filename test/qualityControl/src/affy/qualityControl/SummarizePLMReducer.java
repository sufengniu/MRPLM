package affy.qualityControl;

import org.apache.commons.lang3.StringUtils;
import org.apache.hadoop.io.Text;
import org.apache.hadoop.mapreduce.Reducer;
import org.apache.hadoop.mapreduce.lib.output.MultipleOutputs;

import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

public class SummarizePLMReducer extends
Reducer<Text, Text, Text, Text>{
    private MultipleOutputs<Text, Text> mos = null;
    Text outvalue = new Text();

    @Override
    protected void setup(Context context) throws IOException,
            InterruptedException {
        // Create a new MultipleOutputs using the context object
        mos = new MultipleOutputs<Text, Text>(context);
    }

	PLM plm = new PLM();
	// key: probeset name, value: celname \t intensities
	@Override
	protected void reduce(Text key, Iterable<Text> values, Context context)
			throws IOException, InterruptedException {
		double LOG2 = Math.log(2.0);

		List <String> data = new ArrayList<String>();
		// calculate the number of chips
		for (Text value : values){
			data.add(value.toString());
		}
		int numChips = data.size();
		int numProbes = data.get(0).split("\t")[1].split(",").length;

		double [][] z = new double [numProbes][numChips];

		double [] results = new double [numChips];
		double [] SEresults = new double [numChips];
		String [] celNames = new String [numChips];

		// data.size() = number of chips
		for (int i=0; i<data.size(); i++){
			//String [] strings = data.get(i).split("\t");
            String[] strings = StringUtils.split(data.get(i), '\t');
			celNames[i] = strings[0];
			//String [] probes = strings[1].split(",");
            String[] probes = StringUtils.split(strings[1], ',');
			//probes.length = number of perfect match probes in a probe group
			for (int j=0; j<probes.length; j++) {
				z[j][i] = Math.log(Double.parseDouble(probes[j])) / LOG2;
			}
		}

		double [] affinities = null;
		plm.PLMsummarize(z,numProbes, numChips, results, SEresults, affinities);
		
		StringBuilder outputExpr = new StringBuilder();
		StringBuilder outputSEexpr = new StringBuilder();
		for(int i=0; i<results.length; i++) {
			outputExpr.append(celNames[i]+":"+results[i] + ",");
			outputSEexpr.append(celNames[i]+":"+SEresults[i] + ",");
		}
        outvalue.set(outputSEexpr.toString());
        mos.write("sexpr", key, outvalue, "sexpr/part");
        outvalue.set(outputExpr.toString());
        mos.write("expr", key, outvalue, "expr/part");
		//context.write(new Text(key+":ex"), new Text(outputExpr.toString()));
		//context.write(new Text(key+":se"), new Text(outputSEexpr.toString()));
	}

    @Override
    protected void cleanup(Context context) throws IOException,
            InterruptedException {
        // Close multiple outputs!
        mos.close();
    }
}
