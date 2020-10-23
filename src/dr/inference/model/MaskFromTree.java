package dr.inference.model;
import dr.evolution.tree.MutableTreeListener;
import dr.evolution.tree.NodeRef;
import dr.evolution.tree.Tree;
import dr.evolution.util.Taxon;
import dr.evomodel.tree.TreeModel;
import dr.xml.*;

/**
 * @author Alexander Fisher
 */
public class MaskFromTree extends Parameter.Abstract implements ModelListener {

    public static final String MASK_FROM_TREE = "maskFromTree";
    private TreeModel tree;
    private Taxon taxon;
    private boolean ancestralMaskBranchKnown = false;
    private Parameter maskParameter;

    public MaskFromTree(TreeModel tree, Taxon referenceTaxon) {
        this.tree = tree;
        this.taxon = referenceTaxon;
        int numBranches = tree.getNodeCount() - 1;
        this.maskParameter = new Parameter.Default(numBranches, 1.0);
        tree.addModelListener(this);
        tree.addModelRestoreListener(this);
        updateMask();
    }

    Parameter getAncestralMaskBranch() {
        if (!ancestralMaskBranchKnown) { updateMask(); }
        return (maskParameter);
    }

    void updateMask() {
       int nodeNumber =  tree.getTaxonIndex(taxon.getId());
       NodeRef node = tree.getNode(nodeNumber);
       NodeRef root = tree.getRoot();

       while(tree.getParent(node) != root){
           node = tree.getParent(node);
       }

        int maskIndex = node.getNumber();
        maskParameter.setParameterValue(maskIndex, 0.0);
        ancestralMaskBranchKnown = true;
    }


    public int getDimension(){
        return maskParameter.getDimension();
    }

    @Override
    public double getParameterValue(int dim) {
        if (!ancestralMaskBranchKnown) { updateMask(); }
        return maskParameter.getParameterValue(dim);
    }

    @Override
    public void setParameterValue(int dim, double value) {
        throw new RuntimeException("Not yet implemented");
    }

    @Override
    public void setParameterValueQuietly(int dim, double value) {
        throw new RuntimeException("Not yet implemented");
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
    protected void storeValues() {
        throw new RuntimeException("Not yet implemented");
    }

    @Override
    protected void restoreValues() {
        throw new RuntimeException("Not yet implemented");
    }

    @Override
    protected void acceptValues() {
        throw new RuntimeException("Not yet implemented");
    }

    @Override
    protected void adoptValues(Parameter source) {
        throw new RuntimeException("Not yet implemented");
    }

    @Override
    public void modelChangedEvent(Model model, Object object, int index) {
        ancestralMaskBranchKnown=false;
    }

    @Override
    public void modelRestored(Model model) {
        ancestralMaskBranchKnown=false;
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
