package dr.inference.model;
import dr.evolution.tree.NodeRef;
import dr.evolution.util.Taxon;
import dr.evomodel.tree.TreeModel;
import dr.xml.*;

/**
 * @author Alexander Fisher
 */
public class MaskFromTree {
    // needs to take in a tip from the tree
    // needs to be 'essentially' a mask parameter
    // can currently work by just auto masking root branch
    // will need access to the tree
    public static final String MASK_FROM_TREE = "maskFromTree";
//    public static final String TREE = "tree";

    public MaskFromTree() {

    }


    // **************************************************************
    // XMLObjectParser
    // **************************************************************

    public static XMLObjectParser PARSER = new AbstractXMLObjectParser() {

        public String getParserName() {
            return MASK_FROM_TREE;
        }

        public Object parseXMLObject(XMLObject xo) throws XMLParseException {

            TreeModel tree = (TreeModel) xo.getChild(TreeModel.class);

            Taxon referenceTaxon = (Taxon) xo.getChild(Taxon.class);


            MaskFromTree maskFromTree = new MaskFromTree();
            return maskFromTree;
        }

        //************************************************************************
        // AbstractXMLObjectParser implementation
        //************************************************************************

        public XMLSyntaxRule[] getSyntaxRules() { return rules; }

        private XMLSyntaxRule[] rules = new XMLSyntaxRule[] {
//                new StringAttributeRule(TRAIT_NAME, "The name to be given to the trait to be simulated"),
//                AttributeRule.newBooleanRule(CLONE),
//                AttributeRule.newDoubleRule(INITIAL_VALUE),
                new ElementRule(TreeModel.class),
                new ElementRule(Taxon.class)
        };

        public String getParserDescription() {
            return "Simulates a trait on a tree";
        }

        public Class getReturnType() { return MaskFromTree.class; }
    };

}
