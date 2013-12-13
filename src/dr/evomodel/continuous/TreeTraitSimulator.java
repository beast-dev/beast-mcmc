package dr.evomodel.continuous;

import dr.evolution.tree.NodeRef;
import dr.evolution.tree.Tree;
import dr.evolution.util.Taxon;
import dr.evomodel.tree.TreeModel;
import dr.xml.*;

/**
 * Created by IntelliJ IDEA.
 * User: mandevgill
 */
public class TreeTraitSimulator {

    public static final String SIMULATE_TREE_TRAITS = "treeTraitSimulator";
    public static final String TREE_TRAIT_NAME = "traitName";

    public TreeTraitSimulator(TreeTraitNormalDistributionModel treeTraitModel, TreeModel traitTree, String traitName) {
        this.traitTree = traitTree;
        this.treeTraitModel = treeTraitModel;
        this.traitName = traitName;
    }


    public void simulate() {

        double[] simValues = (double[]) treeTraitModel.nextRandom();
        int numTaxa = traitTree.getTaxonCount();
        int dimTrait = treeTraitModel.getDimTrait();
        double temp[] = new double[dimTrait];

//        for(int i = 0; i < numTaxa; i++){
//            for(int j = 0; j < dimTrait; j++){
//                temp[j] = simValues[i*dimTrait + j];
//            }
//            traitTree.setTaxonAttribute(i,traitName,temp);
//        }

        System.err.println("Playing with generated values...");
        Tree treeModel = treeTraitModel.getTree();

        for (int i = 0; i < numTaxa; ++i) {
            NodeRef node = treeModel.getExternalNode(i);
            Taxon taxon = treeModel.getNodeTaxon(node);
            //   System.err.print("taxon: " + taxon.getId());

            StringBuilder sb = new StringBuilder();
            for (int j = 0; j < dimTrait; ++j) {
                if (j > 0) sb.append(" ");
                sb.append(simValues[i * dimTrait + j]);
            }
            taxon.setAttribute(traitName, sb.toString());
            //     System.err.println(" : " + sb.toString());
            //      System.err.println(" : " + taxon.getAttribute(traitName));
        }
        //     System.err.println("Wasn't that fun?");


//        tree.
        //       System.exit(-1);
//		return traitTree;
    }


    // **************************************************************
    // XMLObjectParser
    // **************************************************************

    public static XMLObjectParser PARSER = new AbstractXMLObjectParser() {

        public String getParserName() {
            return SIMULATE_TREE_TRAITS;
        }

        public Object parseXMLObject(XMLObject xo) throws XMLParseException {
//            System.err.println("Getting here!");
//            System.exit(-1);

            String traitName = xo.getStringAttribute(TREE_TRAIT_NAME);
            TreeTraitNormalDistributionModel treeTraitModel = (TreeTraitNormalDistributionModel) xo.getChild(TreeTraitNormalDistributionModel.class);
            TreeModel traitTree = (TreeModel) xo.getChild(TreeModel.class);

            TreeTraitSimulator treeTraitSimulator = new TreeTraitSimulator(treeTraitModel, traitTree, traitName);
//            return
            treeTraitSimulator.simulate();
            return treeTraitSimulator;
        }

        //************************************************************************
        // AbstractXMLObjectParser implementation
        //************************************************************************

        public XMLSyntaxRule[] getSyntaxRules() {
            return rules;
        }

        private XMLSyntaxRule[] rules = new XMLSyntaxRule[]{
                new StringAttributeRule(TREE_TRAIT_NAME, "The name to be given to the trait to be simulated"),
                new ElementRule(TreeTraitNormalDistributionModel.class),
                new ElementRule(TreeModel.class)
        };

        public String getParserDescription() {
            return "Simulates a trait on a tree";
        }

        public Class getReturnType() {
            return TreeTraitSimulator.class;
        }
    };

    TreeModel traitTree = null;
    TreeTraitNormalDistributionModel treeTraitModel = null;
    String traitName = null;
}
