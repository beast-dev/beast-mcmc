package test.dr.evomodel.tree;

import dr.evolution.io.Importer;
import dr.evolution.io.NewickImporter;
import dr.evolution.tree.NodeRef;
import dr.evolution.tree.Tree;
import dr.evolution.tree.TreeUtils;
import dr.evolution.util.Taxon;
import dr.evomodel.tree.ConstraintsTreeLikelihood;
import dr.evomodel.tree.TreeModel;
import junit.framework.TestCase;

import java.io.IOException;


public class ConstraintsTreeLikelihoodTest extends TestCase {
        public void setUp() throws Exception {
            super.setUp();
            String newickString = "((1:1.0,(7:0.5,8:0.5,9:0.5):0.5,3:1.0,2:1.0):1.0,(4:0.5,5:0.5,6:0.5):0.5);";

            NewickImporter importer = new NewickImporter("(((1:1.0,3:1.0):1.0,2:2.0):1.0,4:3.0);");
            NewickImporter constraintsImporter = new NewickImporter("((1:1.0,3:1.0,2:1.0):1.0,4:1.0);");

            Tree tree = importer.importTree(null);
            Tree constraintsTree = constraintsImporter.importTree(null);

             targetTree = new TreeModel(tree);

            constraintsTreeLikelihood = new ConstraintsTreeLikelihood("MYSTAT",targetTree,constraintsTree);
        }

        public void testShouldPass() {
            assertEquals(0.0,constraintsTreeLikelihood.getLogLikelihood());
        }
        public void testShouldStillPass(){
            constraintsTreeLikelihood.makeDirty();
            assertEquals(0.0,constraintsTreeLikelihood.getLogLikelihood());
        }

        public void testShouldFail()  {


            Taxon selectedTaxon1 = targetTree.getTaxon(targetTree.getTaxonIndex("4"));
            Taxon selectedTaxon2 = targetTree.getTaxon(targetTree.getTaxonIndex("1"));
            NodeRef selectedNode1=null;
            NodeRef selectedNode2=null;
            for (int j = 0; j < targetTree.getExternalNodeCount(); j++) {
                NodeRef tip = targetTree.getExternalNode(j);
                if (targetTree.getNodeTaxon(tip).equals(selectedTaxon1)) {
                    selectedNode1 = tip;
                }else if(targetTree.getNodeTaxon(tip).equals(selectedTaxon2)){
                    selectedNode2=tip;
                }
            }

            NodeRef parent1 = targetTree.getParent(selectedNode1);
            NodeRef parent2 = targetTree.getParent(selectedNode2);

            targetTree.beginTreeEdit();
            targetTree.removeChild(parent1, selectedNode1);
            targetTree.removeChild(parent2,selectedNode2);

            targetTree.addChild(parent1,selectedNode2);
            targetTree.addChild(parent2,selectedNode1);
            targetTree.endTreeEdit();


            assertEquals(Double.NEGATIVE_INFINITY, constraintsTreeLikelihood.getLogLikelihood());
        }

        public void testWithNestedClades() throws IOException, Importer.ImportException, TreeUtils.MissingTaxonException {
            String newickString = "((1:1.0,(7:0.5,8:0.5,9:0.5):0.5,3:1.0,2:1.0):1.0,(4:0.5,5:0.5,6:0.5):0.5);";

            NewickImporter constraintsImporter = new NewickImporter(newickString);

            Tree constraintsTree = constraintsImporter.importTree(null);

            targetTree = new TreeModel(constraintsTree);

            constraintsTreeLikelihood = new ConstraintsTreeLikelihood("MYSTAT",targetTree,constraintsTree);
            assertEquals(0.0,constraintsTreeLikelihood.getLogLikelihood());

        }

        public void testToFailWithNestedClades() throws IOException, Importer.ImportException, TreeUtils.MissingTaxonException {
            String newickString = "((1:1.0,(7:0.5,8:0.5,9:0.5):0.5,3:1.0,2:1.0):1.0,(4:0.5,5:0.5,6:0.5):0.5);";

            NewickImporter constraintsImporter = new NewickImporter(newickString);

            Tree constraintsTree = constraintsImporter.importTree(null);

            targetTree = new TreeModel(constraintsTree);

            constraintsTreeLikelihood = new ConstraintsTreeLikelihood("MYSTAT",targetTree,constraintsTree);

            Taxon selectedTaxon1 = targetTree.getTaxon(targetTree.getTaxonIndex("8"));
            Taxon selectedTaxon2 = targetTree.getTaxon(targetTree.getTaxonIndex("1"));
            NodeRef selectedNode1=null;
            NodeRef selectedNode2=null;
            for (int j = 0; j < targetTree.getExternalNodeCount(); j++) {
                NodeRef tip = targetTree.getExternalNode(j);
                if (targetTree.getNodeTaxon(tip).equals(selectedTaxon1)) {
                    selectedNode1 = tip;
                }else if(targetTree.getNodeTaxon(tip).equals(selectedTaxon2)){
                    selectedNode2=tip;
                }
            }

            NodeRef parent1 = targetTree.getParent(selectedNode1);
            NodeRef parent2 = targetTree.getParent(selectedNode2);

            targetTree.beginTreeEdit();
            targetTree.removeChild(parent1, selectedNode1);
            targetTree.removeChild(parent2,selectedNode2);

            targetTree.addChild(parent1,selectedNode2);
            targetTree.addChild(parent2,selectedNode1);
            targetTree.endTreeEdit();


            assertEquals(Double.NEGATIVE_INFINITY, constraintsTreeLikelihood.getLogLikelihood());


        }




        private ConstraintsTreeLikelihood constraintsTreeLikelihood;
        private TreeModel targetTree;
    }

