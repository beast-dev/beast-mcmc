/*
 * ThorneyTreeLikelihoodTest.java
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

package test.dr.evomodel.bigfasttree.thorney;

import java.io.IOException;

import dr.evolution.io.Importer;
import dr.evolution.io.NewickImporter;
import dr.evolution.tree.NodeRef;
import dr.evolution.tree.Tree;
import dr.evolution.tree.TreeUtils;
import dr.evomodel.bigfasttree.BigFastTreeModel;
import dr.evomodel.bigfasttree.thorney.BranchLengthLikelihoodDelegate;
import dr.evomodel.bigfasttree.thorney.ConstrainedTreeBranchLengthProvider;
import dr.evomodel.bigfasttree.thorney.ConstrainedTreeModel;
import dr.evomodel.bigfasttree.thorney.ConstrainedTreeOperator;
import dr.evomodel.bigfasttree.thorney.MutationBranchMap;
import dr.evomodel.bigfasttree.thorney.PoissonBranchLengthLikelihoodDelegate;
import dr.evomodel.bigfasttree.thorney.ThorneyDataLikelihoodDelegate;
import dr.evomodel.operators.ExchangeOperator;
import dr.evomodel.branchratemodel.BranchRateModel;
import dr.evomodel.branchratemodel.StrictClockBranchRates;
import dr.evomodel.tree.TreeModel;
import dr.evomodel.treedatalikelihood.TreeDataLikelihood;
import dr.inference.model.Parameter;
import dr.inference.operators.AdaptationMode;
import dr.math.distributions.PoissonDistribution;
import junit.framework.TestCase;

public class ThorneyTreeLikelihoodTest extends TestCase {
    public void setUp() throws Exception {
        super.setUp();
        rateParameter = new Parameter.Default(1.0);
        branchRateModel = new StrictClockBranchRates(rateParameter);

        NewickImporter importer = new NewickImporter("((1:1.0,2:1.0,3:2.0):1.0,4:1.0);"); //mutations
        NewickImporter importer2 = new NewickImporter("(((1:1.0,2:1.0):0.1,3:1.1):1.0,4:1.0);");

        tree = importer.importTree(null);
        TreeModel baseTreeModel = new BigFastTreeModel(importer2.importTree(null));

        constrainedTreeModel = new ConstrainedTreeModel("testTree",baseTreeModel,tree);
        MutationBranchMap constrainedTreeBranchLengthProvider = new ConstrainedTreeBranchLengthProvider(constrainedTreeModel,tree);
        BranchLengthLikelihoodDelegate branchLengthLikelihoodDelegate = new PoissonBranchLengthLikelihoodDelegate("poissonDelgate",1.0);
        thorneyDataLikelihoodDelegate = new ThorneyDataLikelihoodDelegate(constrainedTreeModel,constrainedTreeBranchLengthProvider, branchLengthLikelihoodDelegate);

        treeDataLikelihood = new TreeDataLikelihood(thorneyDataLikelihoodDelegate, constrainedTreeModel, branchRateModel);

        expectedLL= 0;
        double[] expectations = {1d,1d,1.1,1.0,1.0,0.1};
        double[] mutations = {1d, 1d, 2.0, 1.0,1.0, 0}; 
        for (int i = 0; i < expectations.length; i++) {
            PoissonDistribution p = new PoissonDistribution(expectations[i]);
            expectedLL += p.logPdf(mutations[i]);
        }
        treeDataLikelihood.getLogLikelihood();
    }

    public void testLikelihood() {

        assertEquals(expectedLL, treeDataLikelihood.getLogLikelihood(),1E-13);
    }

    // TODO broken
//   public void testAfterTopologyChange(){
//       TreeModel subtree = constrainedTreeModel.getSubtree(constrainedTreeModel.getNode(0));
//
//        NodeRef subtreeRoot = subtree.getRoot();// should be inserted node parent
//        NodeRef tip2 = subtree.getExternalNode(1);
//        NodeRef insertedNode = subtree.getParent(tip2);
//        NodeRef tip3 = subtree.getExternalNode(2);
//        assert subtree.getParent(insertedNode) == subtreeRoot;
//
//        constrainedTreeModel.beginTreeEdit();
//        subtree.removeChild(insertedNode,tip2);
//
//        subtree.addChild(subtreeRoot, tip2);
//        subtree.addChild(insertedNode, tip3);
//        constrainedTreeModel.endTreeEdit();
//
//        double ll= 0;
//        double[] expectations = {1d,1.1,1d,1.0,1.0,0.1};//swap tips 2 and 3
//        double[] mutations = {1d, 1d, 2.0, 1.0,1.0, 0};
//       for (int i = 0; i < expectations.length; i++) {
//           PoissonDistribution p = new PoissonDistribution(expectations[i]);
//           ll += p.logPdf(mutations[i]);
//       }
//       assertEquals(ll, treeDataLikelihood.getLogLikelihood(),1E-13);
//   }


    public void testAfterHeightChange(){
        NodeRef insertedNode = constrainedTreeModel.getParent(constrainedTreeModel.getNode(0));
       constrainedTreeModel.setNodeHeight(insertedNode,0.9);
       double ll= 0;
        double[] expectations = {0.9,0.9,1.1,1.0,1.0,0.2};
        double[] mutations = {1d, 1d, 2.0, 1.0,1.0, 0.0}; 
       for (int i = 0; i < expectations.length; i++) {
           PoissonDistribution p = new PoissonDistribution(expectations[i]);
           ll += p.logPdf(mutations[i]);
       }
       assertEquals(ll, treeDataLikelihood.getLogLikelihood(),1E-13);

    }
    public void testAfterRateChange(){
        rateParameter.setValue(0, 2.0);
        double ll= 0;
        double[] expectations = {1d,1d,1.1,1.0,1.0,0.1};
        double[] mutations = {1d, 1d, 2.0, 1.0,1.0, 0}; 
        for (int i = 0; i < expectations.length; i++) {
            PoissonDistribution p = new PoissonDistribution(expectations[i]*2);
            ll += p.logPdf(mutations[i]);
        }
       assertEquals(ll, treeDataLikelihood.getLogLikelihood(),1E-13);

    }
    public void testAfterOperator(){

        ExchangeOperator narrow = new ExchangeOperator(0, null, 10);
        ConstrainedTreeOperator op = new ConstrainedTreeOperator(constrainedTreeModel,10,narrow,1.0,1, AdaptationMode.ADAPTATION_OFF,0.2);

        op.doOperation();
        System.out.println(constrainedTreeModel.toString());
        double LL = treeDataLikelihood.getLogLikelihood();
        treeDataLikelihood.makeDirty();
        double newLL = treeDataLikelihood.getLogLikelihood();
        // NO error?
        assertEquals(newLL,LL,1E-13 );
        
    }


   public void testRootPolytomy() throws IOException, Importer.ImportException, TreeUtils.MissingTaxonException {
       branchRateModel = new StrictClockBranchRates(new Parameter.Default(1.0));

       NewickImporter importer = new NewickImporter("((1:1.0,2:1.0):1.0,3:1.0,4:1.0);");
       NewickImporter importer2 = new NewickImporter("(((1:2.0,2:1.0):1.0,3:1.1):0.1,4:1.0);");

      tree = importer.importTree(null);
        TreeModel baseTreeModel = new BigFastTreeModel(importer2.importTree(null));

        constrainedTreeModel = new ConstrainedTreeModel("testTree",baseTreeModel,tree);
        MutationBranchMap constrainedTreeBranchLengthProvider = new ConstrainedTreeBranchLengthProvider(constrainedTreeModel,tree);
        BranchLengthLikelihoodDelegate branchLengthLikelihoodDelegate = new PoissonBranchLengthLikelihoodDelegate("poissonDelgate",1.0);
        StrictClockBranchRates branchRateModel = new StrictClockBranchRates(new Parameter.Default(1.0));
        thorneyDataLikelihoodDelegate = new ThorneyDataLikelihoodDelegate(constrainedTreeModel,constrainedTreeBranchLengthProvider, branchLengthLikelihoodDelegate);

        treeDataLikelihood = new TreeDataLikelihood(thorneyDataLikelihoodDelegate, constrainedTreeModel, branchRateModel);

       expectedLL= 0;
       double[] expectations = {2d,1d,1d,1d,1.1,0.1};
       double[] mutations = {1d, 1d, 1d,1.0, 1,0}; // time
       for (int i = 0; i < expectations.length; i++) {
           PoissonDistribution p = new PoissonDistribution(expectations[i]);
           expectedLL += p.logPdf(mutations[i]);
       }

       double LL = treeDataLikelihood.getLogLikelihood();


       assertEquals(LL,expectedLL,1E-13);

   }

    public void testRootPolytomyHeightChange() throws IOException, Importer.ImportException, TreeUtils.MissingTaxonException {
        branchRateModel = new StrictClockBranchRates(new Parameter.Default(1.0));

        NewickImporter importer = new NewickImporter("((1:1.0,2:1.0):1.0,3:1.0,4:1.0);");
        NewickImporter importer2 = new NewickImporter("(((1:2.0,2:1.0):1.0,3:1.1):0.1,4:1.0);");

        tree = importer.importTree(null);
        TreeModel baseTreeModel = new BigFastTreeModel(importer2.importTree(null));

        constrainedTreeModel = new ConstrainedTreeModel("testTree",baseTreeModel,tree);
        MutationBranchMap constrainedTreeBranchLengthProvider = new ConstrainedTreeBranchLengthProvider(constrainedTreeModel,tree);
        BranchLengthLikelihoodDelegate branchLengthLikelihoodDelegate = new PoissonBranchLengthLikelihoodDelegate("poissonDelgate",1.0);
        StrictClockBranchRates branchRateModel = new StrictClockBranchRates(new Parameter.Default(1.0));
        thorneyDataLikelihoodDelegate = new ThorneyDataLikelihoodDelegate(constrainedTreeModel,constrainedTreeBranchLengthProvider, branchLengthLikelihoodDelegate);

        treeDataLikelihood = new TreeDataLikelihood(thorneyDataLikelihoodDelegate, constrainedTreeModel, branchRateModel);
        treeDataLikelihood.getLogLikelihood();

        NodeRef rootNode = constrainedTreeModel.getRoot();
        constrainedTreeModel.setNodeHeight(rootNode,3.5);// was 3.1


        double LL = treeDataLikelihood.getLogLikelihood();
        expectedLL= 0;
        double[] expectations = {2d,1d,1d,1.1,1.4,0.5};
        double[] mutations = {1d, 1d, 1d,1.0, 1,0}; // time
        for (int i = 0; i < expectations.length; i++) {
            PoissonDistribution p = new PoissonDistribution(expectations[i]);
            expectedLL += p.logPdf(mutations[i]);
        }

       assertEquals(LL,expectedLL,1E-13);

    }

    private Tree tree;
    private ConstrainedTreeModel constrainedTreeModel;
    private BranchRateModel branchRateModel;
    private ThorneyDataLikelihoodDelegate thorneyDataLikelihoodDelegate;
    private TreeDataLikelihood treeDataLikelihood;
    private double expectedLL;
    private Parameter rateParameter;
}

