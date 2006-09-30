/*
 * SimcoalRunner.java
 *
 * Copyright (C) 2002-2006 Alexei Drummond and Andrew Rambaut
 *
 * This file is part of BEAST.
 * See the NOTICE file distributed with this work for additional
 * information regarding copyright ownership and licensing.
 *
 * BEAST is free software; you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as
 * published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.
 *
 *  BEAST is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with BEAST; if not, write to the
 * Free Software Foundation, Inc., 51 Franklin St, Fifth Floor,
 * Boston, MA  02110-1301  USA
 */

package dr.app.simcoal;

import dr.evolution.tree.Tree;
import dr.evolution.tree.NodeRef;
import dr.evolution.io.NexusImporter;
import dr.evolution.io.Importer;
import dr.evolution.util.Taxon;
import dr.stats.DiscreteStatistics;

import java.io.*;
import java.text.NumberFormat;
import java.text.DecimalFormat;
import java.util.HashSet;

/**
 * Created by IntelliJ IDEA.
 * User: adru001
 * Date: May 30, 2006
 * Time: 5:38:44 PM
 * To change this template use File | Settings | File Templates.
 */
public class SimcoalRunner {

    static NumberFormat format = new DecimalFormat("0.000000");


    static void createParFile(
            String fileName,
            int[] popSize,
            int[] sampleSize,
            double migrationRate,
            double mutationRate,
            double transitionBias,
            int sequenceLength,
            double shapeParameter,
            int numCategories) throws IOException {

        PrintWriter writer = new PrintWriter(new FileWriter(fileName));

        writer.write("// input parameters for the coalescence simulation program: simcoal\n");
        writer.write(popSize.length + " samples to simulate\n");
        writer.write("// Deme sizes (haploid number of genes)\n");
        for (int i = 0; i < popSize.length; i++) {
            writer.write(popSize[i]+"\n");
        }
        writer.write("// Sample sizes\n");
        for (int i = 0; i < popSize.length; i++) {
            writer.write(sampleSize[i]+"\n");
        }
        writer.write("// Growth rates\n");
        for (int i = 0; i < popSize.length; i++) {
            writer.write("0\n");
        }
        writer.write("// Number of migration matrices\n");
        writer.write("1\n");


        writer.write("// Migration rates matrix 0 : A linear stepping-stone model\n");
        for (int i = 0; i < popSize.length; i++) {
            for (int j = 0; j < popSize.length; j++) {
                if (j == i + 1 || j == i-1) {
                    writer.write(format.format(migrationRate) + " ");
                } else {
                    writer.write(format.format(0.0) + " ");
                }
            }
            writer.write("\n");
        }
        writer.write("// Historical event: time, source, sink, migrants, new deme size, new growth rate, new migration matrix\n");
        writer.write("0 historical events\n");
        writer.write("// Mutation rate per generation for the whole sequence\n");
        writer.write(format.format(mutationRate) + "\n");
        writer.write("// Sequence length\n");
        writer.write(sequenceLength + "\n");
        writer.write("// Data type : either DNA, RFLP, or MICROSAT : If DNA, second term is the transition bias\n");
        writer.write("DNA " + transitionBias + "\n");
        writer.write("//Mutation rates gamma distribution shape parameter\n");
        writer.write(shapeParameter + " " + numCategories + "\n");
        writer.close();

    }

    static void createSteppingStoneParFile(
            String fileName,
            int numPopulations,
            int demeSize,
            int sampleSize,
            double migrationRate,
            double mutationRate,
            double transitionBias,
            int sequenceLength,
            double shapeParameter) throws IOException {

        int[] populationSizes = new int[numPopulations];
        int[] sampleSizes = new int[numPopulations];

        for (int i = 0; i < populationSizes.length; i++) {
            populationSizes[i] = demeSize;
            sampleSizes[i] = sampleSize;
        }

        createParFile(fileName, populationSizes, sampleSizes, migrationRate, mutationRate, transitionBias, sequenceLength, shapeParameter, 10);
    }

    static void consumeStream(final InputStream stream) {

        (new Thread() {

            public void run() {
                InputStreamReader in = new InputStreamReader(stream);
                //BufferedReader reader = new BufferedReader(in);

                try {
                    while (in.ready()) {
                        char c = (char)in.read();
                        //System.out.print(c);
                    }
                } catch (IOException e) {
                    e.printStackTrace();  //To change body of catch statement use File | Settings | File Templates.
                }

                /*String line;
                try {
                    while ((line = reader.readLine())
                              != null) {
                        System.out.println(line);
                    }
                } catch (IOException e) {
                    e.printStackTrace();  //To change body of catch statement use File | Settings | File Templates.
                }
*/
            }
        }).start();

    }

    static Tree[] runSimcoal( String simcoalLocation, String fileName, int replicates) throws Exception, InterruptedException, Importer.ImportException {

        Process process = null;
        try {
            String command = simcoalLocation + " " + fileName + " " + replicates;

        String[] cmdArray = new String[] {simcoalLocation, fileName, replicates+""};

        //System.out.println("executing command '" + cmdArray[0] + " " + cmdArray[1] + " " + cmdArray[2] + "'");

        process = Runtime.getRuntime().exec(cmdArray);

        consumeStream(process.getInputStream());
        consumeStream(process.getErrorStream());

        process.waitFor();

        //System.out.println("command '" + command + "' exited with value " + process.exitValue());

        // read in tree file

        String treeFileName = fileName.substring(0,fileName.length()-4) + "_true_trees.trees";

        NexusImporter importer = new NexusImporter(new FileReader(treeFileName));

        Tree[] trees = new Tree[replicates];
        int count = 0;
        while (importer.hasTree()) {
            trees[count] = importer.importNextTree();
            count += 1;
        }

        return trees;

        } catch (Exception e) {
            throw e;

        } finally {
            if (process != null) process.destroy();
        }

    }

    private static double[] processTree(Tree tree) {

        int n = tree.getExternalNodeCount();

        double[] x = new double[n * (n-1) / 2];
        double[] y = new double[x.length];

        int count = 0;
        for (int i = 0; i < tree.getExternalNodeCount(); i++) {
            Taxon taxon = tree.getTaxon(i);
            String id = taxon.getId();

            int location = Integer.parseInt(id.substring(id.indexOf('.')+1));

            for (int j = i+1; j < tree.getExternalNodeCount(); j++) {
                Taxon taxon2 = tree.getTaxon(j);
                String id2 = taxon2.getId();

                int location2 = Integer.parseInt(id2.substring(id2.indexOf('.')+1));

                HashSet set = new HashSet();
                set.add(id);
                set.add(id2);

                NodeRef ancestor = Tree.Utils.getCommonAncestorNode(tree, set);

                if (ancestor == null) throw new RuntimeException("what?!");

                int popDistance = Math.abs(location2 -location);

                double geneticDistance = getHeightFromBranchLengths(tree, ancestor);

                x[count] = popDistance;
                y[count] = geneticDistance;

                //System.out.println(popDistance + "\t" + geneticDistance);
                count += 1;
            }
        }

        LeastSquaresFunction lsf = new LeastSquaresFunction(x,y);
        double[] params = lsf.optimize();
        //System.out.println(format.format(params[0]) + "\t" + format.format(params[1]));
        return params;

    }

    static double getHeightFromBranchLengths(Tree tree, NodeRef node) {
        double height = 0.0;
        while (!tree.isExternal(node)) {
            node = tree.getChild(node, 0);
            height += tree.getBranchLength(node);
        }
        return height;
    }


    public static void main(String[] args) throws Exception, Importer.ImportException, InterruptedException {
        String simcoalLocation = args[0];

        double[] migrationRate = new double[] {
                0.02, 0.01, 0.005, 0.002, 0.001, 0.0005, 0.0002
        };

        int[] sampleSizes = new int[] { 1,2,5,10};

        int sampleSize = sampleSizes[0];

        int[] populationSizes = new int[] {200,500,1000};

        System.out.println("pop\tss\tmig\tm\tmlow\tmupp\tc\tclow\tcupp");

        for (int j = 0; j < populationSizes.length; j++) {
        for (int k = 0; k < migrationRate.length; k++) {

            String filename = "in.par";

            createSteppingStoneParFile(filename,10,populationSizes[j],sampleSize,migrationRate[k], 0.001, 1.0, 1000, 100.0);

            Tree[] trees = runSimcoal(simcoalLocation,filename,100);

            double[] m = new double[trees.length];
            double[] c = new double[trees.length];

            for (int i = 0; i < trees.length; i++) {
                double[] params = processTree(trees[i]);
                m[i] = params[0];
                c[i] = params[1];
            }

            System.out.print(populationSizes[j] + "\t");
            System.out.print(sampleSize + "\t");
            System.out.print(format.format(migrationRate[k]) + "\t");

            System.out.print(format.format(DiscreteStatistics.mean(m)) + "\t");
            System.out.print(format.format(DiscreteStatistics.quantile(0.025,m)) + "\t");
            System.out.print(format.format(DiscreteStatistics.quantile(0.975,m)) + "\t");

            System.out.print(format.format(DiscreteStatistics.mean(c)) + "\t");
            System.out.print(format.format(DiscreteStatistics.quantile(0.025,c)) + "\t");
            System.out.print(format.format(DiscreteStatistics.quantile(0.975,c)) + "\t");

            System.out.println();
        }
        }

    }
}
