package test.dr.evomodel.treelikelihood.thorneytreelikelihood;

import dr.evolution.io.NewickImporter;
import dr.evolution.tree.Tree;

import dr.evomodel.bigfasttree.*;
import dr.evomodel.operators.ExchangeOperator;
import dr.evomodel.treelikelihood.thorneytreelikelihood.*;
import dr.evomodel.branchratemodel.BranchRateModel;
import dr.evomodel.branchratemodel.StrictClockBranchRates;
import dr.evomodel.tree.TreeModel;
import dr.inference.model.Parameter;
import dr.inference.operators.AdaptationMode;
import dr.math.distributions.PoissonDistribution;
import junit.framework.TestCase;

public class ThorneyTreeLikelihoodTest extends TestCase {
    public void setUp() throws Exception {
        super.setUp();
        branchRateModel = new StrictClockBranchRates(new Parameter.Default(1.0));

        NewickImporter importer = new NewickImporter("((1:1.0,2:1.0,3:1.0):1.0,4:1.0);");
        NewickImporter importer2 = new NewickImporter("(((1:1.0,2:1.0):0.1,3:1.1):1.0,4:1.0);");

        tree = importer.importTree(null);
        TreeModel baseTreeModel = new BigFastTreeModel(importer2.importTree(null));

        constrainedTreeModel = new ConstrainedTreeModel("testTree",baseTreeModel,tree);
        BranchLengthProvider constrainedTreeBranchLengthProvider = new ConstrainedTreeBranchLengthProvider(constrainedTreeModel,tree);
        ThorneyBranchLengthLikelihoodDelegate thorneyBranchLengthLikelihoodDelegate = new PoissonBranchLengthLikelihoodDelegate("strictClockDelegate",new StrictClockBranchRates(new Parameter.Default(1.0)),1.0);
        thorneyTreeLikelihood = new ThorneyTreeLikelihood("testLikelihood",constrainedTreeModel,constrainedTreeBranchLengthProvider, thorneyBranchLengthLikelihoodDelegate);

        expectedLL= 0;
        double[] expectations = {1d,1d,1.1,1.0,1.0,0.1};
        double[] mutations = {1d, 1d, 1.0, 1.0,1.0, 0}; // time
        for (int i = 0; i < expectations.length; i++) {
            PoissonDistribution p = new PoissonDistribution(expectations[i]);
            expectedLL += p.logPdf(mutations[i]);
        }
        thorneyTreeLikelihood.getLogLikelihood();
    }

    public void testLikelihood() {

        assertEquals(expectedLL, thorneyTreeLikelihood.getLogLikelihood(),1E-13);
    }
//      This works in person but fails jUnit why?

//    public void testAfterHeightChange(){
//        NodeRef insertedNode = constrainedTreeModel.getParent(constrainedTreeModel.getNode(0));
//        constrainedTreeModel.setNodeHeight(insertedNode,0.9);
//        double ll= 0;
//        double[] expectations = {0.9,0.9,1.1,2d,0.2};
//        double[] mutations = {1d, 1d, 1d, 2d, 0};
//        for (int i = 0; i < expectations.length; i++) {
//            PoissonDistribution p = new PoissonDistribution(expectations[i]);
//            ll += p.logPdf(mutations[i]);
//        }
//        assertEquals(ll, thorneyTreeLikelihood.getLogLikelihood(),1E-13);
//    }
    public void testAfterTopologyChange(){

        ExchangeOperator narrow = new ExchangeOperator(0, null, 10);
        ConstrainedTreeOperator op = new ConstrainedTreeOperator(constrainedTreeModel,10,narrow,1.0,1, AdaptationMode.ADAPTATION_OFF,0.2);

        op.doOperation();
        System.out.println(constrainedTreeModel.toString());
        double LL = thorneyTreeLikelihood.getLogLikelihood();
        thorneyTreeLikelihood.makeDirty();
        double newLL = thorneyTreeLikelihood.getLogLikelihood();
        // NO error?
        assertEquals(newLL,LL,1E-13 );

    }
/*
    public void testAfterPolytomyRootChange() throws TreeUtils.MissingTaxonException {

        TreeModel subtree = constrainedTreeModel.getSubtree(0);

        NodeRef rootNode = subtree.getRoot();
        NodeRef polytomyRoot0 = subtree.getNode(5);
        NodeRef polytomyRoot1 = subtree.getNode(4);
        NodeRef tip1 = subtree.getNode(0);

        double polytomyRoot0height = subtree.getNodeHeight(polytomyRoot0);
        double polytomyRoot1height = subtree.getNodeHeight(polytomyRoot1);

        subtree.setNodeHeight(polytomyRoot0, polytomyRoot1height);
        subtree.setNodeHeight(polytomyRoot1, polytomyRoot0height);


        subtree.beginTreeEdit();
        subtree.removeChild(rootNode, polytomyRoot0);
        subtree.removeChild(polytomyRoot0, polytomyRoot1);
        subtree.removeChild(polytomyRoot1,tip1);

        subtree.addChild(rootNode, polytomyRoot1);
        subtree.addChild(polytomyRoot1,polytomyRoot0);
        subtree.addChild(polytomyRoot0, tip1);
        subtree.endTreeEdit();

//        cladeModel.setRootNode(cladeModel.getClade(polytomyRoot1),polytomyRoot1);

        double LL = thorneyTreeLikelihood.getLogLikelihood();
        thorneyTreeLikelihood.makeDirty();
        double newLL = thorneyTreeLikelihood.getLogLikelihood();

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

       thorneyTreeLikelihood = new ThorneyTreeLikelihood("approximateTreeLikelihood",
               tree,
               29903,
               treeModel,
               branchRateModel
       );;

       thorneyTreeLikelihood.getLogLikelihood();
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

       double LL = thorneyTreeLikelihood.getLogLikelihood();
       thorneyTreeLikelihood.makeDirty();
       double newLL = thorneyTreeLikelihood.getLogLikelihood();

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

        thorneyTreeLikelihood = new ThorneyTreeLikelihood("approximateTreeLikelihood",
                tree,
                29903,
                treeModel,
                branchRateModel
        );;
        thorneyTreeLikelihood.getLogLikelihood();
        expectedLL= 0;
        double[] expectations = {1d,1d,1d,1.1,1.1};
        double[] mutations = {1d, 1d, 1d,1.0, 1}; // time
        for (int i = 0; i < expectations.length; i++) {
            PoissonDistribution p = new PoissonDistribution(expectations[i]);
            expectedLL += p.logPdf(mutations[i]);
        }

        NodeRef rootNode = treeModel.getRoot();
        treeModel.setNodeHeight(rootNode,2.05);

        double LL = thorneyTreeLikelihood.getLogLikelihood();
        thorneyTreeLikelihood.makeDirty();
        double newLL = thorneyTreeLikelihood.getLogLikelihood();

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

        thorneyTreeLikelihood = new ThorneyTreeLikelihood("approximateTreeLikelihood",
                tree,
                29903,
                treeModel,
                branchRateModel
        );;

        thorneyTreeLikelihood.getLogLikelihood();

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


        double LL = thorneyTreeLikelihood.getLogLikelihood();
        thorneyTreeLikelihood.makeDirty();
        double newLL = thorneyTreeLikelihood.getLogLikelihood();

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

        thorneyTreeLikelihood = new ThorneyTreeLikelihood("approximateTreeLikelihood",
                tree,
                29903,
                treeModel,
                branchRateModel
        );;

        thorneyTreeLikelihood.getLogLikelihood();

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

        double LL = thorneyTreeLikelihood.getLogLikelihood();
        thorneyTreeLikelihood.makeDirty();
        double newLL = thorneyTreeLikelihood.getLogLikelihood();

        assertEquals(LL,newLL);
    }
    */

    private Tree tree;
    private ConstrainedTreeModel constrainedTreeModel;
    private BranchRateModel branchRateModel;
    private ThorneyTreeLikelihood thorneyTreeLikelihood;
    private double expectedLL;
}

