package dr.inference.model;
import dr.evolution.tree.NodeRef;
import dr.evolution.util.Taxon;
import dr.evomodel.tree.TreeModel;
import dr.xml.*;

/**
 * @author Alexander Fisher
 */
public class MaskFromTree {

    public static final String MASK_FROM_TREE = "maskFromTree";
    private TreeModel tree;
    private Taxon taxon;
    private boolean ancestralMaskBranchKnown = false;
    private Parameter mask;

    public MaskFromTree(TreeModel tree, Taxon referenceTaxon) {
        this.tree = tree;
        this.taxon = referenceTaxon;
        int numBranches = tree.getNodeCount() - 1;
        this.mask = new Parameter.Default(numBranches, 1.0);
        updateMask();
    }

    Parameter getAncestralMaskBranch() {
        if (!ancestralMaskBranchKnown) { updateMask(); }
        return (mask);
    }

    void updateMask() {
       int nodeNumber =  tree.getTaxonIndex(taxon.getId());
       NodeRef node = tree.getNode(nodeNumber);
       NodeRef root = tree.getRoot();

       while(tree.getParent(node) != root){
           node = tree.getParent(node);
       }

        int maskIndex = node.getNumber();
        mask.setParameterValue(maskIndex, 0.0);
        ancestralMaskBranchKnown = true;
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


            MaskFromTree maskFromTree = new MaskFromTree(tree, referenceTaxon);
            return maskFromTree;
        }

        //************************************************************************
        // AbstractXMLObjectParser implementation
        //************************************************************************

        public XMLSyntaxRule[] getSyntaxRules() { return rules; }

        private XMLSyntaxRule[] rules = new XMLSyntaxRule[] {
                new ElementRule(TreeModel.class),
                new ElementRule(Taxon.class)
        };

        public String getParserDescription() {
            return "Masks ancestral (off-root) branch to a specific reference taxon on a random tree";
        }

        public Class getReturnType() { return MaskFromTree.class; }
    };

}
