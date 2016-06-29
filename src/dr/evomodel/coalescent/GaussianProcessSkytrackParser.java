package dr.evomodel.coalescent;

import dr.stats.DiscreteStatistics;

import java.io.BufferedReader;
import java.io.FileReader;
import java.util.ArrayList;
import java.util.Arrays;

/**
 * Created by mkarcher on 6/22/16.
 */
public class GaussianProcessSkytrackParser {
    private static final int changepointsIndex = 5;
    private static final int GvaluesIndex = 6;
    private static final int precisionIndex = 8;
    private static final int tmrcaIndex = 10;

    public static void main(String[] args) {
        String filePathName = "examples/hcvNew2small.log";
        CSVstats stats = parseCSV(filePathName, 3);
        System.out.println(Arrays.toString(stats.tmrcas));
        System.out.println(DiscreteStatistics.median(stats.tmrcas));
    }

    public static CSVstats parseCSV(String filename, int skip) {
        CSVstats result = null;

        try {
            BufferedReader br = new BufferedReader(new FileReader(filename));

            String line;
            String[] fieldsArray;
            int lineNum = 0;

            ArrayList<ArrayList<Double>> changepoints = new ArrayList<ArrayList<Double>>();
            ArrayList<ArrayList<Double>> Gvalues = new ArrayList<ArrayList<Double>>();
            ArrayList<Double> signals = new ArrayList<Double>();
            ArrayList<Double> tmrcas = new ArrayList<Double>();

            while ((line = br.readLine()) != null) {
                if (lineNum >= skip) {
                    fieldsArray = line.split("\t");

                    changepoints.add(parseListStr(fieldsArray[changepointsIndex]));
                    Gvalues.add(parseListStr(fieldsArray[GvaluesIndex]));
                    signals.add(new Double(fieldsArray[precisionIndex]));
                    tmrcas.add(new Double(fieldsArray[tmrcaIndex]));
                }
                lineNum += 1;
            }

            result = new CSVstats(changepoints, Gvalues, signals, tmrcas);
        } catch (Exception e) {
            e.printStackTrace();
        }

        return result;
    }

    public static void gpPosterior(CSVstats stats) {

    }

    private static ArrayList<Double> parseListStr(String arg) {
        ArrayList<Double> vals = new ArrayList<Double>();

        String[] valStrs = arg.replaceAll("\\{|}", "").split(",");
        for (String valStr : valStrs) {
            vals.add(new Double(valStr));
        }

        return vals;
    }

    private static double sigmoidal(double x) {
        return 1/(1+Math.exp(-x));
    }

    private static class CSVstats {
        public ArrayList<double[]> changepoints, Gvalues;
        public double[] signals, tmrcas;

        public CSVstats(ArrayList<ArrayList<Double>> changepoints,
                        ArrayList<ArrayList<Double>> Gvalues,
                        ArrayList<Double> signals,
                        ArrayList<Double> tmrcas) {
            this.changepoints = new ArrayList<double[]>();
            for (ArrayList<Double> al : changepoints) {
                double[] arr = new double[al.size()];
                for (int i = 0; i < al.size(); i++) {
                    arr[i] = al.get(i);
                }
                this.changepoints.add(arr);
            }

            this.Gvalues = new ArrayList<double[]>();
            for (ArrayList<Double> al : Gvalues) {
                double[] arr = new double[al.size()];
                for (int i = 0; i < al.size(); i++) {
                    arr[i] = al.get(i);
                }
                this.Gvalues.add(arr);
            }

            this.signals = new double[signals.size()];
            for (int i = 0; i < signals.size(); i++) {
                this.signals[i] = signals.get(i);
            }

            this.tmrcas = new double[tmrcas.size()];
            for (int i = 0; i < tmrcas.size(); i++) {
                this.tmrcas[i] = tmrcas.get(i);
            }
        }
    }
}
