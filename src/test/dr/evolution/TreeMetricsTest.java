/*
 * TreeMetricsTest.java
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

package test.dr.evolution;

import dr.evolution.io.Importer;
import dr.evolution.io.NewickImporter;
import dr.evolution.tree.*;
import dr.evolution.tree.treemetrics.BranchScoreMetric;
import dr.evolution.tree.treemetrics.SteelPennyPathDifferenceMetric;
import dr.evolution.tree.treemetrics.RobinsonFouldsMetric;
import junit.framework.TestCase;

import java.io.*;

/**
 * @author Luiz Carvalho
 */
public class TreeMetricsTest extends TestCase {

    public static void main(String[] args) {

        try {

            NewickImporter importer = new NewickImporter("((A:0.1,B:0.1):0.1,(C:0.1,D:0.1):0.1)");
            Tree treeOne = importer.importNextTree();
            System.out.println("tree 1: " + treeOne);

            importer = new NewickImporter("(((A:0.1,B:0.1):0.5,C:0.1):0.1,D:0.1)");
            Tree treeTwo = importer.importNextTree();
            System.out.println("tree 2: " + treeTwo + "\n");
            
            /* Billera et al., 2001 */
            
//            double billera = (new BilleraMetric().getMetric(TreeUtils.asJeblTree(treeOne),
//            		TreeUtils.asJeblTree(treeTwo)));
//
//            System.out.println("Billera distance = " + billera);
//            assertEquals(billera, 0.2236068);
            
            /* Robinson & Foulds, 1981*/
            double RF = (new RobinsonFouldsMetric().getMetric(treeOne, treeTwo)*2.0);
            System.out.println("Robinson-Foulds = " + RF);
            assertEquals(RF, 2.0, 0.0000001);
            
            /* Penny and Hendy, 1993*/
            double path = (new SteelPennyPathDifferenceMetric().getMetric(treeOne, treeTwo));
            System.out.println("path difference = " + path);
            assertEquals(path, 0.7141428, 0.0000001);
            
            /* Branch Score*/
            double bl = (new BranchScoreMetric().getMetric(treeOne, treeTwo));
            System.out.println("bl score = " + bl);
            assertEquals(bl, Math.sqrt(Math.pow(0.5-0.1, 2) + Math.pow(.1, 2) + Math.pow(.1, 2)), 0.0000001);
            

        } catch(Importer.ImportException ie) {
            System.err.println(ie);
        } catch(IOException ioe) {
            System.err.println(ioe);
        }

    }

}

