package affy.qualityControl;

import org.apache.commons.lang3.StringUtils;
import org.apache.hadoop.conf.Configuration;
import org.apache.hadoop.fs.FileSystem;
import org.apache.hadoop.fs.Path;
import org.apache.hadoop.io.Text;
import org.apache.hadoop.mapreduce.Mapper;

import java.io.BufferedReader;
import java.io.DataInputStream;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

public class RNAdegMapper extends
    Mapper<Text, Text, Text, Text> {
	private int N;
    private String indexFile;
    private Text outvalue = new Text();

    protected void setup(Context context) throws IOException,
            InterruptedException {
        indexFile = context.getConfiguration().get("index_file");
    }

	// key cel file name, value intensity.
    @Override
	protected void map(Text key, Text value, Context context)
			throws IOException, InterruptedException {
		Path inFile = new Path(indexFile);
		FileSystem fs = FileSystem.get(new Configuration());

		// split value into the intensity
		//String [] intensity = value.toString().split(",");
        String[] intensities = StringUtils.split(value.toString(), ',');
		//output format: celname$intensity[i]

		ArrayList<ArrayList<Double> > pm = new ArrayList<ArrayList<Double> >();
		if (!fs.exists(inFile)){// probeIndex.txt exist?
			System.err.println("Input file: "+indexFile+" not found!");
			System.exit(-1);
		}
		else {// read probeIndex.txt
			BufferedReader reader = new BufferedReader(new InputStreamReader(
                    new DataInputStream(fs.open(inFile))
            ));
			String line;
			int ps = 0; //the counter of the probe set
			while ((line = reader.readLine())!=null) {
				pm.add(new ArrayList<Double> ());
				//strings[0]: probe group name, strings[1]: index of probe in intensity
				//String [] strings = line.split("\t");
                String[] strings = StringUtils.split(line, '\t');
				//String [] indices = strings[1].split(",");
                String[] indices = StringUtils.split(strings[1], ',');
				StringBuilder output = new StringBuilder();
				output.append(key+"\t");
                for (String index : indices) {
                    pm.get(ps).add(Double.parseDouble(intensities[Integer.parseInt(index)]));
                }
				ps++;
			}
			reader.close();

			List<Double> sd = new ArrayList<Double>();
			List<Double> mean = new ArrayList<Double>();
			getMeanSd(pm, mean, sd);

			double first = mean.get(0);
            double sqrt_N = Math.sqrt(N);
			for (int j=0; j<mean.size(); j++){
				mean.set(j, mean.get(j)-first);
				mean.set(j, mean.get(j)/(sd.get(j)/sqrt_N));
			}

			ArrayList<Double> indices = new ArrayList<Double> ();
			for (int i=0; i<mean.size(); i++)
				indices.add(i+0.0);
			LinearModel lm = new LinearModel();
			lm.linearRegression(indices, mean);

//			//context.write(key, new Text(String.valueOf(lm.getSlope()+" "+lm.getSlopeProT())));
//			if (lm.getSlope()>4.5){
//                outvalue.set("rnaDeg");
//				context.write(key, outvalue);
//            }
            outvalue.set(lm.getSlope()+"\trnaDeg");
            context.write(key, outvalue);
        }
	}

	private void getMeanSd(ArrayList<ArrayList<Double> > pm, List<Double> mean, List<Double> sd){
		ArrayList<ArrayList<Double> > data1 = new ArrayList<ArrayList<Double>> ();
		Map<Integer, Integer> maxNum = new HashMap<Integer, Integer>();
		for (int i=0; i<pm.size(); i++){
			data1.add(new ArrayList<Double>());
			for (int j=0; j<pm.get(i).size(); j++){
				data1.get(i).add(log2(pm.get(i).get(j)));
			}
			int count = 0;
			if (maxNum.get(data1.get(i).size())!=null)
				count = maxNum.get(data1.get(i).size());
			// key: number of probes in a probeset, value: number of probeset have the same key
			maxNum.put(data1.get(i).size(), count + 1);
		}

		int key = 0; //the number of probes in a probeset
		int value = 0; //the count of the probesets
		for (Map.Entry<Integer, Integer> pair : maxNum.entrySet()){
			if (value < pair.getValue()){
				value = pair.getValue();
				key = pair.getKey();
			}
		}
		//new matrix to store probesets with largest number.
		ArrayList<ArrayList<Double> > data2 = new ArrayList<ArrayList<Double>> (value);
		int n = 0;
        for (ArrayList<Double> dataArray : data1) {
            if (dataArray.size() != key) {
                continue;
            }
            data2.add(new ArrayList<Double>());
            for (Double data : dataArray) {
                data2.get(n).add(data);
            }
            n++;
        }
		N = n;
		//ArrayList<Double> mns = new ArrayList<Double>(value);
		//ArrayList<Double> sds = new ArrayList<Double>(value);

		for (int i=0; i<key; i++){
			ArrayList<Double> probes = new ArrayList<Double>();
			for (int j=0; j<value; j++)
				probes.add(data2.get(j).get(i));
			mean.add(mean(probes));
			sd.add(sd(probes));
		}
	}

	private double log2(double value){
		return Math.log10(value)/Math.log10(2);
	}

	public double sd (List<Double> a){
		double sum = 0;
		double mean = mean(a);

		for (double i : a)
			sum += Math.pow((i - mean), 2);
		return Math.sqrt( sum / ( a.size() - 1 ) ); // sample
	}

	public double mean (List<Double> a){
        return sum(a) / (a.size() * 1.0);
    }

	public double sum (List<Double> a){
        double sum = 0;
		if (a.size() > 0) {
			for (double i : a) {
				sum += i;
			}
		}
		return sum;
	}

}


