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

public class SummarizePLMMapper extends
Mapper<Text, Text, Text, Text> {

    private String indexFile;
    private Text outkey = new Text();
    private Text outvalue = new Text();

    protected void setup(Context context) throws IOException,
            InterruptedException {
        indexFile = context.getConfiguration().get("index_file");
    }

	// key cel file name, value intensities.
	protected void map(Text key, Text value, Context context)
			throws IOException, InterruptedException {
		Path inFile = new Path(indexFile);
		FileSystem fs = FileSystem.get(new Configuration());

		// split value into the intensity
		//String [] intensities = value.toString().split(",");
        String[] intensities = StringUtils.split(value.toString(), ',');
		// output format: celname$intensity[i]

		if (!fs.exists(inFile)){// probeIndex.txt exist?
			System.err.println("Input file: "+indexFile+" not found!");
			System.exit(-1);
		}
		else {// read probeIndex.txt
			BufferedReader br = new BufferedReader(new InputStreamReader(
                    new DataInputStream(fs.open(inFile))
            ));
			String line = null;
			while ((line = br.readLine())!=null) {
				// strings[0]: probe group name, strings[1]: index of probe in intensity
				//String [] strings = line.split("\t");
                String[] strings = StringUtils.split(line, '\t');
				//String [] indices = strings[1].split(",");
                String[] indices = StringUtils.split(strings[1], ',');
				StringBuilder output = new StringBuilder();
				output.append(key+"\t");
                for (String indice : indices) {
                    output.append(intensities[Integer.parseInt(indice)] + ",");
                }
				// key: probe group name, value: celname \t intensities
                outkey.set(strings[0]);
                outvalue.set(output.toString());
                context.write(outkey, outvalue);
			}
			br.close();
		}
	}
}
