/*
 * GhostTreeModelTest.java
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

package test.dr.evomodel.bigfasttree;

import dr.evolution.io.NewickImporter;
import dr.evolution.tree.Tree;
import dr.evolution.util.Taxa;
import dr.evolution.util.Taxon;
import dr.evomodel.bigfasttree.ghosttree.GhostTreeModel;
import dr.evomodel.operators.SubtreeLeapOperator;
import dr.evomodel.tree.TreeModel;
import dr.inference.operators.AdaptationMode;
import dr.math.MathUtils;
import junit.framework.TestCase;

public class GhostTreeModelTest extends TestCase {
    private GhostTreeModel ghostTree;
    private TreeModel corporealTree;
    public void setUp() throws Exception {
        super.setUp();
        Taxa ghostTaxa = new Taxa("ghost");
        ghostTaxa.addTaxon(new Taxon("g1"));
        ghostTaxa.addTaxon(new Taxon("g2"));

        Taxa corporealTaxa = new Taxa("corporeal");
        corporealTaxa.addTaxon(new Taxon("0"));
        corporealTaxa.addTaxon(new Taxon("1"));
        corporealTaxa.addTaxon(new Taxon("2"));

        Taxa allTaxa = new Taxa(ghostTaxa);
        allTaxa.addTaxa(corporealTaxa);

        NewickImporter importer = new NewickImporter("(((g1:0.5,0:1):0.5,1:1):1,(g2:1,2:1):0.5)");
        Tree superTree = importer.importTree(allTaxa);

        ghostTree = new GhostTreeModel(superTree, ghostTaxa);
        corporealTree = ghostTree.getCorporealTreeModel();
    }
    private void checkState(){
        System.out.println(ghostTree);
        System.out.println(corporealTree);
//        String currentText = corporealTree.toString();
//        ghostTree.makeDirty();
//        String textFromScratch = corporealTree.toString();
//        assertEquals(textFromScratch,currentText);
    }

    public void testSetup(){
        checkState();

    }
    public void testOperation(){
        MathUtils.setSeed(1);
        checkState();
        SubtreeLeapOperator op = new SubtreeLeapOperator(ghostTree,1,0.5, SubtreeLeapOperator.DistanceKernelType.CAUCHY, AdaptationMode.ADAPTATION_OFF,0.5);

        for (int i = 0; i < 100; i++) {
            op.doOperation();
            checkState();
        }

    }

}


