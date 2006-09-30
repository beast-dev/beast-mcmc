/*
 * ArisRejection.java
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

package dr.app.misc;

import dr.evolution.io.Importer;
import dr.evolution.io.NexusImporter;
import dr.evolution.tree.SimpleTree;
import dr.evolution.tree.Tree;
import dr.math.ExponentialDistribution;

import java.io.FileReader;
import java.io.IOException;

/**
 * Rejection sampling of KRP trees for posterior estimation of parameter values from
 * a single input tree.
 *
 *
 * @author Aris Katzourakis
 * @author Alexei Drummond
 *
 * @version $Id: ArisRejection.java,v 1.1 2005/09/26 22:14:15 rambaut Exp $
 */
public class ArisRejection {

    private static final double THRESHOLD = 0.01;

    public ArisRejection() {}


    public void sample(Tree targetTree, int reps) {

//        TreeSummaryStatistic index = new CollessIndex();
//
//        double targetIc = index.getSummaryStatistic(targetTree);
//
//        for (int i = 0; i < reps; i++) {
//
//            double B = sampleBPrior();
//            double D = sampleDPrior();
//            double I = sampleIPrior();
//            double p = Math.random();
//
//            Tree tree = simulateKRP05Tree(B, D, I, p);
//
//            double Ic = index.getSummaryStatistic(tree);
//
//            if (matches(targetIc, Ic)) {
//                // output cool parameter values
//                System.out.println(B + "\t" + D + "\t" + I + "\t" + p);
//            }
//        }
    }

    /**
     * @param ic1
     * @param ic2
     * @return true if the summary statistics are close enough to be regarded as a match
     */
    public boolean matches(double ic1, double ic2) {
        return Math.abs(ic1 - ic2) < THRESHOLD;
    }

    /**
     * @return a sample from the prior distribution shared by B, D, I
     */
    public double sampleBPrior() {

        return ExponentialDistribution.quantile(Math.random(), 1.0);
    }

    /**
     * @return a sample from the prior distribution shared by B, D, I
     */
    public double sampleDPrior() {

        return ExponentialDistribution.quantile(Math.random(), 1.0);
    }

    /**
     * @return a sample from the prior distribution shared by B, D, I
     */
    public double sampleIPrior() {

        return ExponentialDistribution.quantile(Math.random(), 1.0);
    }

    /**
     *
     * @param B
     * @param D
     * @param I
     * @param p
     * @return an KRP05 tree simulated with the given parameter values
     */
    public Tree simulateKRP05Tree(double B, double D, double I, double p) {

        // todo ARIS MUST IMPLEMENT NOW OR ELSE
        return new SimpleTree();
    }

    public static void main(String[] args) throws IOException, Importer.ImportException {

        ArisRejection ar = new ArisRejection();

        NexusImporter importer = new NexusImporter(new FileReader(args[0]));

        Tree targetTree = importer.importNextTree();

        ar.sample(targetTree, 1000000);
    }
}
