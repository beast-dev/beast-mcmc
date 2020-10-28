package dr.inference.model;

import dr.evolution.tree.NodeRef;
import dr.evolution.util.Taxon;
import dr.evomodel.tree.TreeModel;
import dr.xml.*;

/**
 * @author Alexander Fisher
 */
public class MaskFromTree extends Parameter.Default implements ModelListener {

    public static final String MASK_FROM_TREE = "maskFromTree";
    private TreeModel tree;
    private Taxon taxon;
//    private boolean ancestralMaskBranchKnown = false;

    public MaskFromTree(TreeModel tree, Taxon referenceTaxon, int dim) {
        super(dim);
        this.tree = tree;
        this.taxon = referenceTaxon;
        tree.addModelListener(this);
        tree.addModelRestoreListener(this);
        updateMask();
    }

    void updateMask() {
        // todo: make sure this sets the old 0.0 to 1.0 smarter
        for (int i = 0; i < this.getDimension(); i++) {
            setParameterValueQuietly(i, 1.0);
        }
        int nodeNumber = tree.getTaxonIndex(taxon.getId());
        NodeRef node = tree.getNode(nodeNumber);
        NodeRef root = tree.getRoot();

        while (tree.getParent(node) != root) {
            node = tree.getParent(node);
        }

        int maskIndex = node.getNumber();
        this.setParameterValue(maskIndex, 0.0);
    }

    @Override
    public void setParameterValueNotifyChangedAll(int dim, double value) {
        throw new RuntimeException("Not yet implemented");
    }

    @Override
    public String getParameterName() {
        return null;
    }

    @Override
    public void addBounds(Bounds<Double> bounds) {
        throw new RuntimeException("Not yet implemented");
    }

    @Override
    public Bounds<Double> getBounds() {
        return null;
    }

    @Override
    public void addDimension(int index, double value) {
        throw new RuntimeException("Not yet implemented");
    }

    @Override
    public double removeDimension(int index) {
        return 0;
    }

    @Override
    public void modelChangedEvent(Model model, Object object, int index) {
        super.storeParameterValues();
        updateMask();
    }

    @Override
    public void modelRestored(Model model) {
        super.restoreParameterValues();
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

            int numBranches = tree.getNodeCount() - 1;

            MaskFromTree maskFromTree = new MaskFromTree(tree, referenceTaxon, numBranches);
            return maskFromTree;
        }

        //************************************************************************
        // AbstractXMLObjectParser implementation
        //************************************************************************

        public XMLSyntaxRule[] getSyntaxRules() {
            return rules;
        }

        private XMLSyntaxRule[] rules = new XMLSyntaxRule[]{
                new ElementRule(TreeModel.class),
                new ElementRule(Taxon.class)
        };

        public String getParserDescription() {
            return "Masks ancestral (off-root) branch to a specific reference taxon on a random tree";
        }

        public Class getReturnType() {
            return MaskFromTree.class;
        }
    };
}
