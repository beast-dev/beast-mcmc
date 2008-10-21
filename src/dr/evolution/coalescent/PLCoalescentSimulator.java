package dr.evolution.coalescent;

import dr.evolution.tree.Tree;
import dr.evolution.util.Date;
import dr.evolution.util.Taxa;
import dr.evolution.util.Taxon;
import dr.evolution.util.Units;

import java.io.*;
import java.util.ArrayList;
import java.util.List;

/**
 * @author Alexei Drummond
 */
public class PLCoalescentSimulator {

    public static void main(String[] arg) throws IOException {

        // READ DEMOGRAPHIC FUNCTION

        String filename = arg[0];
        BufferedReader reader = new BufferedReader(new FileReader(filename));

        double popSizeScale = 1.0;
        double generationTime = 1.0;
        if (arg.length > 2) {
            popSizeScale = Double.parseDouble(arg[2]);
        }
        if (arg.length > 3) {
            generationTime = Double.parseDouble(arg[3]);
        }

        List<Double> times = new ArrayList<Double>();

        String line = reader.readLine();
        String[] tokens = line.trim().split("[\t ]+");
        if (tokens.length < 2) throw new RuntimeException();

        List<Double>[] popSizes = new List[tokens.length - 1];
        for (int i = 0; i < tokens.length - 1; i++) {
            popSizes[i] = new ArrayList<Double>();
        }

        while (line != null) {

            double time = Double.parseDouble(tokens[0]) / generationTime;

            times.add(time);
            for (int i = 1; i < tokens.length; i++) {
                popSizes[i - 1].add(Double.parseDouble(tokens[i]));
            }
            line = reader.readLine();
            if (line != null) {
                tokens = line.trim().split("[\t ]+");
                if (tokens.length != popSizes.length + 1) throw new RuntimeException();
            }
        }

        reader.close();

        // READ SAMPLE TIMES

        String samplesFilename = arg[1];

        reader = new BufferedReader(new FileReader(samplesFilename));

        line = reader.readLine();
        Taxa taxa = new Taxa();
        int id = 0;
        while (line != null) {

            if (!line.startsWith("#")) {

                tokens = line.split("[\t ]+");

                if (tokens.length == 4) {

                    double t0 = Double.parseDouble(tokens[0]);
                    double t1 = Double.parseDouble(tokens[1]);
                    double dt = Double.parseDouble(tokens[2]);
                    int k = Integer.parseInt(tokens[3]);
                    for (double time = t0; time <= t1; time += dt) {

                        double sampleTime = time / generationTime;
                        for (int i = 0; i < k; i++) {
                            Taxon taxon = new Taxon(id + "");
                            taxon.setAttribute(dr.evolution.util.Date.DATE, new Date(sampleTime, Units.Type.GENERATIONS, true));
                            taxa.addTaxon(taxon);
                            id += 1;
                        }
                    }

                } else {

                    // sample times are in the same units as simulation
                    double sampleTime = Double.parseDouble(tokens[0]) / generationTime;
                    int count = Integer.parseInt(tokens[1]);

                    for (int i = 0; i < count; i++) {
                        Taxon taxon = new Taxon(id + "");
                        taxon.setAttribute(dr.evolution.util.Date.DATE, new Date(sampleTime, Units.Type.GENERATIONS, true));
                        taxa.addTaxon(taxon);
                        id += 1;
                    }
                }
            }
            line = reader.readLine();
        }

        double minTheta = Double.MAX_VALUE;
        double maxTheta = 0.0;

        PrintWriter out = new PrintWriter(System.out);
        if (arg.length > 4) {
            out = new PrintWriter(new FileWriter(arg[4]));
        }

        for (int i = 0; i < popSizes.length; i++) {
            double[] thetas = new double[popSizes[i].size()];
            double[] intervals = new double[times.size() - 1];

            // must reverse the direction of the model
            for (int j = intervals.length; j > 0; j--) {
                intervals[intervals.length - j] = times.get(j) - times.get(j - 1);

                double theta = popSizes[i].get(j) * popSizeScale;
                thetas[intervals.length - j] = theta;
                if (theta < minTheta) {
                    minTheta = theta;
                }
                if (theta > maxTheta) {
                    maxTheta = theta;
                }

                double t = times.get(intervals.length) - times.get(j);

                System.out.println(t + "\t" + theta);
            }

            System.out.println("min theta = " + minTheta);
            System.out.println("max theta = " + maxTheta);

            PiecewiseLinearPopulation demo = new PiecewiseLinearPopulation(intervals, thetas, Units.Type.GENERATIONS);

            CoalescentSimulator simulator = new CoalescentSimulator();
            Tree tree = simulator.simulateTree(taxa, demo);

            out.println(Tree.Utils.newick(tree));
            System.err.println(Tree.Utils.newick(tree));
        }
        out.flush();
        out.close();
    }
}
