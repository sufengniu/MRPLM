import affy.qualityControl.*;
import org.apache.hadoop.conf.Configuration;
import org.apache.hadoop.filecache.DistributedCache;
import org.apache.hadoop.fs.FileSystem;
import org.apache.hadoop.fs.Path;
import org.apache.hadoop.io.DoubleWritable;
import org.apache.hadoop.io.NullWritable;
import org.apache.hadoop.io.Text;
import org.apache.hadoop.mapred.JobConf;
import org.apache.hadoop.mapreduce.Job;
import org.apache.hadoop.mapreduce.lib.input.KeyValueTextInputFormat;
import org.apache.hadoop.mapreduce.lib.input.MultipleInputs;
import org.apache.hadoop.mapreduce.lib.output.FileOutputFormat;
import org.apache.hadoop.mapreduce.lib.output.MultipleOutputs;
import org.apache.hadoop.mapreduce.lib.output.TextOutputFormat;
import org.apache.hadoop.util.GenericOptionsParser;

import java.io.IOException;
import java.net.URI;

/**
 * Created with IntelliJ IDEA.
 * User: David
 * Date: 6/11/13
 * Time: 3:54 PM
 * To change this template use File | Settings | File Templates.
 */
public class qualityControlDriver {
    private static int numReducers;

    private static void setInputParameters(Configuration conf, String[] args) throws IOException {
        String[] otherArgs = new GenericOptionsParser(conf, args)
                .getRemainingArgs();
        if (otherArgs.length != 1) {
            System.err.println("Usage: qualityControlJob" + " <number-of-reduce-tasks>");
            System.exit(1);
        }

        numReducers = Integer.parseInt(otherArgs[0]);
        conf.set("index_file", "part-m-00000");
    }

    public static void main(String[] args) throws Exception {
        JobConf conf = new JobConf();

        setInputParameters(conf, args);
        DistributedCache.createSymlink(conf);
        DistributedCache.addCacheFile(new URI("hdfs://user/sniu/lib/libjniWrapper.so#libjniWrapper.so"), conf);
        conf.set("mapred.reduce.child.java.opts", "-Djava.library.path=/home/sniu/hadoop-dist/app/lib/");

        int code = doPLMSummarizeOnly(conf);
        System.exit(code);
    }


    private static int doPLMSummarizeOnly(Configuration conf) throws InterruptedException, IOException, ClassNotFoundException {
        Job summarizePLMJob = submitSummarizePLMJob(conf);
        while (!summarizePLMJob.isComplete()) {
            Thread.sleep(5000);
        }

        if (summarizePLMJob.isSuccessful()) {
            System.out.println("do PLM summarize job completed successfully!");
        } else {
            System.out.println("do PLM summarize job failed!");
        }


        return (summarizePLMJob.isSuccessful() ? 0 : 1);
    }



    private static Job submitSummarizePLMJob(Configuration conf) throws IOException, ClassNotFoundException, InterruptedException {
        String inputPath = "input";
        String outputPath = "output";
        int numReduceTasks = numReducers;
        Job job = new Job(conf, "do PLM summarization");
        job.setJarByClass(qualityControlDriver.class);
        job.setMapperClass(SummarizePLMMapper.class);
        job.setReducerClass(SummarizePLMReducer.class);
        job.setNumReduceTasks(numReduceTasks);
        job.setOutputKeyClass(Text.class);
        job.setOutputValueClass(Text.class);

        job.setInputFormatClass(KeyValueTextInputFormat.class);
        KeyValueTextInputFormat.setInputPaths(job, new Path(inputPath));

        // Configure multiple outputs
        job.setOutputFormatClass(TextOutputFormat.class);
        FileOutputFormat.setOutputPath(job, new Path(outputPath));

        MultipleOutputs.addNamedOutput(job, "sexpr", TextOutputFormat.class,
                Text.class, Text.class);
        MultipleOutputs.addNamedOutput(job, "expr", TextOutputFormat.class,
                Text.class, Text.class);
        FileSystem.get(conf).delete(new Path(outputPath), true);
        job.submit();
        return job;
    }

}
