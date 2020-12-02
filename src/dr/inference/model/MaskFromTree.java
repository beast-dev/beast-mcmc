package dr.inference.model;

import dr.evolution.tree.NodeRef;
import dr.evolution.util.Taxon;
import dr.evomodel.tree.TreeModel;
import dr.evomodel.tree.TreeParameterModel;
import dr.xml.*;

import java.util.ArrayList;
import java.util.List;

/**
 * @author Alexander Fisher
 */


/**
 * This class provides a mask for branch-specific parameters based on a user-specified descendant taxon.
 * It works by automagically masking the 'root-branch' that is ancestor to the specified taxon as well
 * as an arbitrary number (maskLength) of daughter branches (default = 0)
 */

public class MaskFromTree extends Parameter.Default implements ModelListener {

    public static final String MASK_FROM_TREE = "maskFromTree";
    public static final String MASK_LENGTH = "maskLength";
    private TreeModel tree;
    private Taxon taxon;
    private TreeParameterModel branchRates;
    private int maskLength; //maskLength determines number of 'daughter' branches to mask as well

    public MaskFromTree(TreeModel tree, Taxon referenceTaxon, int dim, int maskLength) {
        super(dim);
        this.tree = tree;
        this.taxon = referenceTaxon;
        this.branchRates = new TreeParameterModel(tree, this, false); //default includeRoot = false
        tree.addModelListener(this);
        tree.addModelRestoreListener(this);
        this.maskLength = maskLength;

        updateMask();
    }

    void updateMask() {
        // todo: make sure this sets the old 0.0 to 1.0 smarter
        // todo: remove code duplication from maskLength > 1
        for (int i = 0; i < this.getDimension(); i++) {
            setParameterValueQuietly(i, 1.0);
        }
        int nodeNumber = tree.getTaxonIndex(taxon.getId());
        NodeRef node = tree.getNode(nodeNumber);
        NodeRef root = tree.getRoot();

        while (tree.getParent(node) != root) {
            node = tree.getParent(node);
        }
        int maskIndex = branchRates.getParameterIndexFromNodeNumber(node.getNumber());
        this.setParameterValue(maskIndex, 0.0);

        List<NodeRef> nodeList = new ArrayList<>();
        // only happens if maskLength > 0 :
        for (int i = 0; i < maskLength; i++) {

            if (!tree.isExternal(node)) {
                node = tree.getChild(node, 0);
                nodeList.add(node);
            }

//            maskIndex = branchRates.getParameterIndexFromNodeNumber(node.getNumber())
        }
        for (int i = 0; i < nodeList.size(); i++) {
            maskIndex = branchRates.getParameterIndexFromNodeNumber(nodeList.get(i).getNumber());
            this.setParameterValue(maskIndex, 0.0);
        }

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

            int maskLength = xo.getAttribute(MASK_LENGTH, 0);

            int numBranches = tree.getNodeCount() - 1;

            MaskFromTree maskFromTree = new MaskFromTree(tree, referenceTaxon, numBranches, maskLength);
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
                new ElementRule(Taxon.class),
                AttributeRule.newStringRule(MASK_LENGTH, true),
        };

        public String getParserDescription() {
            return "Masks ancestral (off-root) branch to a specific reference taxon on a random tree";
        }

        public Class getReturnType() {
            return MaskFromTree.class;
        }
    };
}
