/*
 * TestTreeMetrics.java
 *
 * Copyright © 2002-2024 the BEAST Development Team
 * http://beast.community/about
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
 *
 */

package dr.evolution.tree.treemetrics;/**
 * ${CLASS_NAME}
 *
 * @author Andrew Rambaut
 */

import dr.evolution.io.Importer;
import dr.evolution.io.NewickImporter;
import dr.evolution.tree.Tree;

import java.io.IOException;

public class TestTreeMetrics {

    public static void main(String[] args) {

        try {

            //4-taxa example
            NewickImporter importer = new NewickImporter("(('A':1.2,'B':0.8):0.5,('C':0.8,'D':1.0):1.1)");
            Tree treeOne = importer.importNextTree();
            System.out.println("4-taxa tree 1: " + treeOne);

            importer = new NewickImporter("((('A':0.8,'B':1.4):0.3,'C':0.7):0.9,'D':1.0)");
            Tree treeTwo = importer.importNextTree();
            System.out.println("4-taxa tree 2: " + treeTwo);
            System.out.println();

            System.out.println("Paired trees:");
            System.out.println("BranchScore = " + new BranchScoreMetric().getMetric(treeOne, treeTwo));
            System.out.println("CladeHeight = " + new CladeHeightMetric().getMetric(treeOne, treeTwo));
            System.out.println("Kendall-Colijn lambda (0.0) = " + new KendallColijnPathDifferenceMetric(0.0).getMetric(treeOne, treeTwo));
            System.out.println("Kendall-Colijn lambda (0.5) = " + new KendallColijnPathDifferenceMetric(0.5).getMetric(treeOne, treeTwo));
            System.out.println("Kendall-Colijn lambda (1.0) = " + new KendallColijnPathDifferenceMetric(1.0).getMetric(treeOne, treeTwo));
            System.out.println("RobinsonFoulds = " + new RobinsonFouldsMetric().getMetric(treeOne, treeTwo));
            System.out.println("RootedBranchScore = " + new RootedBranchScoreMetric().getMetric(treeOne, treeTwo));
            System.out.println("Steel-Penny = " + new SteelPennyPathDifferenceMetric().getMetric(treeOne, treeTwo));
            System.out.println();

            //5-taxa example
            importer = new NewickImporter("(((('A':0.6,'B':0.6):0.1,'C':0.5):0.4,'D':0.7):0.1,'E':1.3)");
            treeOne = importer.importNextTree();
            System.out.println("5-taxa tree 1: " + treeOne);

            importer = new NewickImporter("((('A':0.8,'B':1.4):0.1,'C':0.7):0.2,('D':1.0,'E':0.9):1.3)");
            treeTwo = importer.importNextTree();
            System.out.println("5-taxa tree 2: " + treeTwo);
            System.out.println();

            //lambda = 0.0 should yield: sqrt(7) = 2.6457513110645907162
            //lambda = 1.0 should yield: sqrt(2.96) = 1.7204650534085252911

            System.out.println("Paired trees:");
            System.out.println("BranchScore = " + new BranchScoreMetric().getMetric(treeOne, treeTwo));
            System.out.println("CladeHeight = " + new CladeHeightMetric().getMetric(treeOne, treeTwo));
            System.out.println("Kendall-Colijn lambda (0.0) = " + new KendallColijnPathDifferenceMetric(0.0).getMetric(treeOne, treeTwo));
            System.out.println("Kendall-Colijn lambda (0.5) = " + new KendallColijnPathDifferenceMetric(0.5).getMetric(treeOne, treeTwo));
            System.out.println("Kendall-Colijn lambda (1.0) = " + new KendallColijnPathDifferenceMetric(1.0).getMetric(treeOne, treeTwo));
            System.out.println("RobinsonFoulds = " + new RobinsonFouldsMetric().getMetric(treeOne, treeTwo));
            System.out.println("RootedBranchScore = " + new RootedBranchScoreMetric().getMetric(treeOne, treeTwo));
            System.out.println("Steel-Penny = " + new SteelPennyPathDifferenceMetric().getMetric(treeOne, treeTwo));
            System.out.println();


        } catch (Importer.ImportException ie) {
            System.err.println(ie);
        } catch (IOException ioe) {
            System.err.println(ioe);
        }

    }

}
