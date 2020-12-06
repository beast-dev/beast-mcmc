package test.dr.evomodel.bigfasttree;

import dr.evolution.io.Importer;
import dr.evolution.io.NewickImporter;
import dr.evolution.tree.NodeRef;
import dr.evolution.tree.Tree;

import dr.evolution.tree.TreeUtils;
import dr.evomodel.bigfasttree.*;
import dr.evomodel.bigfasttree.constrainedtree.CladeNodeModel;
import dr.evomodel.bigfasttree.constrainedtree.CladeRef;
import dr.evomodel.bigfasttree.constrainedtree.ConstrainedTreeBranchLengthProvider;
import dr.evomodel.branchratemodel.BranchRateModel;
import dr.evomodel.branchratemodel.StrictClockBranchRates;
import dr.evomodel.tree.TreeModel;
import dr.inference.model.Parameter;
import dr.math.distributions.PoissonDistribution;
import junit.framework.TestCase;

import java.io.IOException;

public class ApproximatePoissonTreeLikelihoodTest extends TestCase {
    public void setUp() throws Exception {
        super.setUp();
        branchRateModel = new StrictClockBranchRates(new Parameter.Default(1.0));

        NewickImporter importer = new NewickImporter("((1:1.0,2:1.0,3:1.0):1.0,4:1.0);");
        NewickImporter importer2 = new NewickImporter("(((1:1.0,2:1.0):0.1,3:1.1):1.0,4:1.0);");

        tree = importer.importTree(null);
        treeModel = new BigFastTreeModel(importer2.importTree(null));

        cladeModel = new CladeNodeModel(tree, treeModel);
        BranchLengthProvider constrainedBranchLengthProvider = new ConstrainedTreeBranchLengthProvider(cladeModel);

        approximatePoissonTreeLikelihood = new ApproximatePoissonTreeLikelihood("approximateTreeLikelihood",
                1,
                treeModel,
                branchRateModel,
                constrainedBranchLengthProvider);

        expectedLL= 0;
        double[] expectations = {1d,1d,1.1,2d,0.1};
        double[] mutations = {1d, 1d, 1.0, 2d, 0}; // time
        for (int i = 0; i < expectations.length; i++) {
            PoissonDistribution p = new PoissonDistribution(expectations[i]);
            expectedLL += p.logPdf(mutations[i]);
        }

        approximatePoissonTreeLikelihood.getLogLikelihood();
    }

    public void testLikelihood() {

        assertEquals(expectedLL,approximatePoissonTreeLikelihood.getLogLikelihood(),1E-13);
    }

    public void testAfterHeightChange(){

        NodeRef insertedNode = treeModel.getParent(treeModel.getNode(0));
        treeModel.setNodeHeight(insertedNode,0.9);
        double ll= 0;
        double[] expectations = {0.9,0.9,1.1,2d,0.2};
        double[] mutations = {1d, 1d, 1d, 2d, 0};
        for (int i = 0; i < expectations.length; i++) {
            PoissonDistribution p = new PoissonDistribution(expectations[i]);
            ll += p.logPdf(mutations[i]);
        }
        assertEquals(ll,approximatePoissonTreeLikelihood.getLogLikelihood(),1E-13);
    }
    public void testAfterTopologyChange(){


        NodeRef selectedNode1 = treeModel.getNode(0); // node 1
        NodeRef selectedNode2 = treeModel.getNode(2); // node 3

        NodeRef parent1 = treeModel.getParent(selectedNode1); // clade root
        NodeRef parent2 = treeModel.getParent(selectedNode2);

        treeModel.beginTreeEdit();
        treeModel.removeChild(parent1, selectedNode1);
        treeModel.removeChild(parent2, selectedNode2);

        treeModel.addChild(parent1, selectedNode2);
        treeModel.addChild(parent2, selectedNode1);
        treeModel.endTreeEdit();


        double LL = approximatePoissonTreeLikelihood.getLogLikelihood();
        approximatePoissonTreeLikelihood.makeDirty();
        double newLL = approximatePoissonTreeLikelihood.getLogLikelihood();

        assertEquals(expectedLL,LL,1E-13);
        assertEquals(LL,newLL);

    }

    public void testAfterPolytomyRootChange() throws TreeUtils.MissingTaxonException {

        NodeRef rootNode = treeModel.getRoot();
        NodeRef polytomyRoot0 = treeModel.getNode(5);
        NodeRef polytomyRoot1 = treeModel.getNode(4);
        NodeRef tip1 = treeModel.getNode(0);

        double polytomyRoot0height = treeModel.getNodeHeight(polytomyRoot0);
        double polytomyRoot1height = treeModel.getNodeHeight(polytomyRoot1);

        treeModel.setNodeHeight(polytomyRoot0, polytomyRoot1height);
        treeModel.setNodeHeight(polytomyRoot1, polytomyRoot0height);


        treeModel.beginTreeEdit();
        treeModel.removeChild(rootNode, polytomyRoot0);
        treeModel.removeChild(polytomyRoot0, polytomyRoot1);
        treeModel.removeChild(polytomyRoot1,tip1);

        treeModel.addChild(rootNode, polytomyRoot1);
        treeModel.addChild(polytomyRoot1,polytomyRoot0);
        treeModel.addChild(polytomyRoot0, tip1);
        treeModel.endTreeEdit();

//        cladeModel.setRootNode(cladeModel.getClade(polytomyRoot1),polytomyRoot1);

        double LL = approximatePoissonTreeLikelihood.getLogLikelihood();
        approximatePoissonTreeLikelihood.makeDirty();
        double newLL = approximatePoissonTreeLikelihood.getLogLikelihood();

        assertEquals(expectedLL,LL,1E-13);
        assertEquals(LL,newLL);

    }


   public void testRootPolytomy() throws IOException, Importer.ImportException, TreeUtils.MissingTaxonException {
       branchRateModel = new StrictClockBranchRates(new Parameter.Default(1.0));

       NewickImporter importer = new NewickImporter("((1:1.0,2:1.0):1.0,3:1.0,4:1.0);");
       NewickImporter importer2 = new NewickImporter("(((1:1.0,2:1.0):1.0,3:1.1):0.1,4:1.0);");

       tree = importer.importTree(null);
       treeModel = new BigFastTreeModel(importer2.importTree(null));

       CladeNodeModel cladeModel = new CladeNodeModel(tree, treeModel);
       BranchLengthProvider constrainedBranchLengthProvider = new ConstrainedTreeBranchLengthProvider(cladeModel);

       approximatePoissonTreeLikelihood = new ApproximatePoissonTreeLikelihood("approximateTreeLikelihood",
               1,
               treeModel,
               branchRateModel,
               constrainedBranchLengthProvider);

       approximatePoissonTreeLikelihood.getLogLikelihood();
       expectedLL= 0;
       double[] expectations = {1d,1d,1d,1.1,1.1};
       double[] mutations = {1d, 1d, 1d,1.0, 1}; // time
       for (int i = 0; i < expectations.length; i++) {
           PoissonDistribution p = new PoissonDistribution(expectations[i]);
           expectedLL += p.logPdf(mutations[i]);
       }

       NodeRef rootNode = treeModel.getRoot();
       NodeRef rootNode1 = treeModel.getNode(5);
       NodeRef tip3 = treeModel.getNode(2);



       CladeRef clade = cladeModel.getClade(rootNode);

       treeModel.beginTreeEdit();

       treeModel.removeChild(rootNode, rootNode1);
       treeModel.setRoot(rootNode1);

       treeModel.removeChild(rootNode1,tip3);
       treeModel.addChild(rootNode,tip3);
       treeModel.addChild(rootNode1,rootNode);

       treeModel.setNodeHeight(rootNode1,treeModel.getNodeHeight(rootNode)+1);
       treeModel.endTreeEdit();


//       cladeModel.setRootNode(clade,rootNode1);

       double LL = approximatePoissonTreeLikelihood.getLogLikelihood();
       approximatePoissonTreeLikelihood.makeDirty();
       double newLL =approximatePoissonTreeLikelihood.getLogLikelihood();

       assertEquals(LL,newLL,1E-13);

   }

    public void testRootPolytomyHeightChange() throws IOException, Importer.ImportException, TreeUtils.MissingTaxonException {
        branchRateModel = new StrictClockBranchRates(new Parameter.Default(1.0));

        NewickImporter importer = new NewickImporter("((1:1.0,2:1.0):1.0,3:1.0,4:1.0,5:1.0);");
        NewickImporter importer2 = new NewickImporter("(((1:1.0,2:1.0):1.0,3:1.1):0.1,(4:1.0,5:1.0):0.1));");

        tree = importer.importTree(null);
        treeModel = new BigFastTreeModel(importer2.importTree(null));

        CladeNodeModel cladeModel = new CladeNodeModel(tree, treeModel);
        BranchLengthProvider constrainedBranchLengthProvider = new ConstrainedTreeBranchLengthProvider(cladeModel);

        approximatePoissonTreeLikelihood = new ApproximatePoissonTreeLikelihood("approximateTreeLikelihood",
                1,
                treeModel,
                branchRateModel,
                constrainedBranchLengthProvider);
        approximatePoissonTreeLikelihood.getLogLikelihood();
        expectedLL= 0;
        double[] expectations = {1d,1d,1d,1.1,1.1};
        double[] mutations = {1d, 1d, 1d,1.0, 1}; // time
        for (int i = 0; i < expectations.length; i++) {
            PoissonDistribution p = new PoissonDistribution(expectations[i]);
            expectedLL += p.logPdf(mutations[i]);
        }

        NodeRef rootNode = treeModel.getRoot();
        treeModel.setNodeHeight(rootNode,2.05);

        double LL = approximatePoissonTreeLikelihood.getLogLikelihood();
        approximatePoissonTreeLikelihood.makeDirty();
        double newLL =approximatePoissonTreeLikelihood.getLogLikelihood();

        assertEquals(LL,newLL);

    }

    // This tests updating when the root and rootchild2  change
    //         +- 6
    //       |
    //+------| gp       +- 4 *
    //|      |          |
    //|      +----------| p       + 9
    //|                 |     +---|
    //|                 +-----|   + 7
    //| r                     |
    //|                       +--- 8
    //|
    //|+- 1
    //+|
    // | +----- 3
    // +-|
    //   +--- 2

    // To
    //+-------------------- 4 *
    //|
    //|       +- 6
    //|p      |
    //|+------| gp                + 9
    //||      |               +---|
    //||      +---------------|   + 7
    //+|                      |
    // | r                    +--- 8
    // |
    // | +- 1
    // +-|
    //   | +---- 3
    //   +-|
    //     +-- 2
    //
    //
    public void testRootAndRootChildUpdates() throws IOException, Importer.ImportException, TreeUtils.MissingTaxonException {
        branchRateModel = new StrictClockBranchRates(new Parameter.Default(1.0));

        NewickImporter importer = new NewickImporter("((1:1.0,2:0.5,3:2.0):0.1,4:0.3,(7:0.2,8:0.1,9:0.3):0.2,6:0.01)");
        NewickImporter importer2 = new NewickImporter("((6:1,(4:1,((9:1,7:1):2,8:2):3.0):6):4,(1:1,(3:3,2:2):1):1);");

        tree = importer.importTree(null);
        treeModel = new BigFastTreeModel(importer2.importTree(null));

        CladeNodeModel cladeModel = new CladeNodeModel(tree, treeModel);
        BranchLengthProvider constrainedBranchLengthProvider = new ConstrainedTreeBranchLengthProvider(cladeModel);

        approximatePoissonTreeLikelihood = new ApproximatePoissonTreeLikelihood("approximateTreeLikelihood",
                1,
                treeModel,
                branchRateModel,
                constrainedBranchLengthProvider);

        approximatePoissonTreeLikelihood.getLogLikelihood();

        NodeRef node = treeModel.getNode(9);
        NodeRef sibling = treeModel.getNode(1);
        NodeRef parent = treeModel.getParent(node);
        NodeRef grandparent = treeModel.getParent(parent);
        NodeRef root = treeModel.getRoot();

        CladeRef clade = cladeModel.getClade(node);

        treeModel.beginTreeEdit();

        treeModel.removeChild(grandparent,parent);
        treeModel.removeChild(parent, sibling);

        treeModel.setNodeHeight(parent,treeModel.getNodeHeight(root)+1);
        treeModel.setRoot(parent);
//        cladeModel.setRootNode(clade, parent);
        treeModel.addChild(parent,root);

        treeModel.addChild(grandparent,sibling);

        treeModel.endTreeEdit();


        double LL = approximatePoissonTreeLikelihood.getLogLikelihood();
        approximatePoissonTreeLikelihood.makeDirty();
        double newLL =approximatePoissonTreeLikelihood.getLogLikelihood();

        assertEquals(LL,newLL);
    }

    public void testOtherChildUpdates() throws IOException, Importer.ImportException, TreeUtils.MissingTaxonException {
        branchRateModel = new StrictClockBranchRates(new Parameter.Default(1.0));

        NewickImporter importer = new NewickImporter("((1:1.0,2:0.5,3:2.0):0.1,4:0.3,(7:0.2,8:0.1,9:0.3):0.2,6:0.01)");
        NewickImporter importer2 = new NewickImporter("((6:1,(4:1,((9:1,7:1):2,8:2):3.0):6):4,(1:1,(3:3,2:2):1):1);");

        tree = importer.importTree(null);
        treeModel = new BigFastTreeModel(importer2.importTree(null));

        CladeNodeModel cladeModel = new CladeNodeModel(tree, treeModel);
        BranchLengthProvider constrainedBranchLengthProvider = new ConstrainedTreeBranchLengthProvider(cladeModel);

        approximatePoissonTreeLikelihood = new ApproximatePoissonTreeLikelihood("approximateTreeLikelihood",
                1,
                treeModel,
                branchRateModel,
                constrainedBranchLengthProvider);

        approximatePoissonTreeLikelihood.getLogLikelihood();

        NodeRef node = treeModel.getNode(1);
        NodeRef sibling = treeModel.getNode(9);
        NodeRef parent = treeModel.getParent(node);
        NodeRef grandparent = treeModel.getParent(parent);
        NodeRef root = treeModel.getRoot();

        CladeRef clade = cladeModel.getClade(parent);

        treeModel.beginTreeEdit();

        treeModel.removeChild(grandparent,parent);
        treeModel.removeChild(parent, sibling);

        treeModel.setNodeHeight(parent,treeModel.getNodeHeight(root)+1);
        treeModel.setRoot(parent);
//        cladeModel.setRootNode(clade, parent);
        treeModel.addChild(parent,root);

        treeModel.addChild(grandparent,sibling);

        treeModel.endTreeEdit();

        double LL = approximatePoissonTreeLikelihood.getLogLikelihood();
        approximatePoissonTreeLikelihood.makeDirty();
        double newLL =approximatePoissonTreeLikelihood.getLogLikelihood();

        assertEquals(LL,newLL);
    }

    private Tree tree;
    private TreeModel treeModel;
    private BranchRateModel branchRateModel;
    private ApproximatePoissonTreeLikelihood approximatePoissonTreeLikelihood;
    private CladeNodeModel cladeModel;
    private double expectedLL;
}

