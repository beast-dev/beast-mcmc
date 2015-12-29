/*
 * LineageCountThroughTime.java
 *
 * Copyright (c) 2002-2015 Alexei Drummond, Andrew Rambaut and Marc Suchard
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

package dr.evolution.tree;

import dr.stats.Variate;
import jebl.evolution.graphs.Node;
import jebl.evolution.io.ImportException;
import jebl.evolution.io.NewickImporter;
import jebl.evolution.io.NexusImporter;
import jebl.evolution.io.TreeImporter;
import jebl.evolution.trees.RootedTree;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

/**
 * @author Alexei Drummond
 * @author Andrew Rambaut
 */
public class LineageCountThroughTime {

    public static Variate[] getLTT(String treeFile,
                                   double minTime,
                                   double maxTime,
                                   int binCount,
                                   int skip // number of trees to skip
    ) throws IOException, ImportException {

        double delta = (maxTime - minTime) / (binCount - 1);

        BufferedReader reader = new BufferedReader(new FileReader(treeFile));

        String line = reader.readLine();

        TreeImporter importer;
        if (line.toUpperCase().startsWith("#NEXUS")) {
            importer = new NexusImporter(reader);
        } else {
            importer = new NewickImporter(reader, false);
        }

        int state = 0;

        while (importer.hasTree() && state < skip) {
            importer.importNextTree();
            state += 1;
        }

        // the age of the end of this group
        List<double[]> branchingTimes = new ArrayList<double[]>();
        state = 0;

        while (importer.hasTree()) {
            RootedTree tree = (RootedTree) importer.importNextTree();

            double[] bt = new double[tree.getInternalNodes().size()];

            int i = 0;
            for (Node node : tree.getInternalNodes()) {
                bt[i] = tree.getHeight(node);
                i++;
            }
            Arrays.sort(bt);

            branchingTimes.add(bt);
            state += 1;
        }

        Variate[] bins = new Variate[binCount];
        double height = 0;

        double n = branchingTimes.get(0).length;

        for (int k = 0; k < binCount; k++) {
            bins[k] = new Variate.D();

            if (height >= 0.0 && height <= maxTime) {
                for (state = 0; state < branchingTimes.size(); state++) {
                    int index = 0;
                    while (index < branchingTimes.get(state).length && branchingTimes.get(state)[index] < height) {
                        index += 1;
                    }

                    double lineageCount = 1;
                    if (index < branchingTimes.get(state).length) {
                        lineageCount = n - index + 1;
                    }
                    bins[k].add(lineageCount);
                }
            }
            height += delta;
        }

        Variate xData = new Variate.D();
        Variate yDataMean = new Variate.D();
        Variate yDataMedian = new Variate.D();
        Variate yDataUpper = new Variate.D();
        Variate yDataLower = new Variate.D();

        double t = minTime;
        for (Variate bin : bins) {
            xData.add(t);
            if (bin.getCount() > 0) {
                yDataMean.add(bin.getMean());
                yDataMedian.add(bin.getQuantile(0.5));
                yDataLower.add(bin.getQuantile(0.025));
                yDataUpper.add(bin.getQuantile(0.975));
            } else {
                yDataMean.add(Double.NaN);
                yDataMedian.add(Double.NaN);
                yDataLower.add(Double.NaN);
                yDataUpper.add(Double.NaN);
            }
            t += delta;
        }

        return new Variate[]{xData, yDataMean};
    }
}
