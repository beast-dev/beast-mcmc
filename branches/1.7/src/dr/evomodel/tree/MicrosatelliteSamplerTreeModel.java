package dr.evomodel.tree;

import dr.inference.model.*;
import dr.evolution.tree.NodeRef;
import dr.evolution.util.Taxon;
import dr.evolution.alignment.Patterns;

import java.util.Map;

/**
 * @author Chieh-Hsi Wu
 *
 * Tree model that stores the probabilites of each branch when the ancestral states are sampled.
 *
 */


public class MicrosatelliteSamplerTreeModel extends AbstractModel {
    protected final TreeModel tree;

    // The tree parameter;
    protected Parameter parameter;

    // the index of the root node.
    protected int rootNodeNumber;
    protected int storedRootNodeNumber;
    private int[] externalValues;
    private double[] logBranchLikelihoods;
    private double[] storedLogBranchLikelihoods;
    private Patterns microsatPat;
    Map<String, Integer> taxaMap;

    /**
     * This class constructs a tree parameter, and will set the dimension of the parameter
     * to match the appropriate number of nodes if necessary.
     *
     * @param id                    model id
     * @param tree                  tree that this parameter corresponds to
     * @param internalValues        values to keep in sync with tree topology moves.
     * @param microsatPat           microsatellite patten
     * @param externalValues        values at the tips of the tree stored in an array.
     * @param taxaMap               maps the taxa to the index of array storing the external values
     * @param internalStateProvided indicates whether internal states has being specified
     *
     */
    public MicrosatelliteSamplerTreeModel(
           String id,
           TreeModel tree,
           Parameter internalValues,
           Patterns microsatPat,
           int[] externalValues,
           Map<String, Integer> taxaMap,
           boolean internalStateProvided) {
        super(id);

        this.tree = tree;
        this.parameter = internalValues;
        this.microsatPat = microsatPat;
        this.externalValues = externalValues;
        this.taxaMap = taxaMap;
        this.logBranchLikelihoods = new double[tree.getNodeCount()];
        this.storedLogBranchLikelihoods = new double[tree.getNodeCount()];

        rootNodeNumber = tree.getRoot().getNumber();

        storedRootNodeNumber = rootNodeNumber;
        if(!internalStateProvided){
            initialiseInternalStates();
        }
        int dim = parameter.getDimension();
        int treeSize = tree.getInternalNodeCount();

        if (dim != treeSize) {
            System.err.println("WARNING: setting dimension of parameter to match tree branch count");
            parameter.setDimension(treeSize);
            Bounds<Double> bounds = parameter.getBounds();
            int upper = (int)(double)bounds.getUpperLimit(0);
            int lower = (int)(double)bounds.getLowerLimit(0);
            parameter.addBounds(new Parameter.DefaultBounds(upper, lower, treeSize));
        }


        addModel(tree);
        addVariable(parameter);
   }


    public void handleModelChangedEvent(Model model, Object object, int index) {}

    protected final void handleVariableChangedEvent(Variable variable, int index, Parameter.ChangeType type) {
        int nodeNumber = getNodeNumberFromParameterIndex(index);

        NodeRef node = tree.getNode(nodeNumber);
        assert (node.getNumber() == nodeNumber) :
                "node.getNumber()=" + node.getNumber() + ", whereas nodeNumber=" + nodeNumber;

        internalNodesChanged = true;
        fireModelChanged(this, nodeNumber);
    }

    boolean internalNodesChanged = false;
    public boolean areInternalNodesChanged(){
        return internalNodesChanged;
    }

    public void setInternalNodesChanged(boolean changed){
         internalNodesChanged = changed;
    }

    protected void storeState() {
        storedRootNodeNumber = rootNodeNumber;
        System.arraycopy(logBranchLikelihoods, 0, storedLogBranchLikelihoods, 0, logBranchLikelihoods.length);
    }

    protected void restoreState() {
        rootNodeNumber = storedRootNodeNumber;
        double[] temp1 = logBranchLikelihoods;
        logBranchLikelihoods = storedLogBranchLikelihoods;
        storedLogBranchLikelihoods = temp1;
    }

    protected void acceptState() {}

    /**
     *
     * @param node  Node reference number
     *
     * @return the state of the node
     *
     */
    public int getNodeValue(NodeRef node) {
        if(tree.isExternal(node)){
            Taxon nodeTaxon = tree.getNodeTaxon(node);
            int externalValIndex = taxaMap.get(nodeTaxon.getId());
            return externalValues[externalValIndex];
        }

        int nodeNumber = node.getNumber();
        int index = getParameterIndexFromNodeNumber(nodeNumber);
        return (int)parameter.getParameterValue(index);

    }

    public String getBranchAttributeLabel() {
        return parameter.getId();
    }

    public String getAttributeForBranch(NodeRef node) {
        return Double.toString(getNodeValue(node));
    }

    public int getNodeNumberFromParameterIndex(int parameterIndex) {
        return parameterIndex + tree.getExternalNodeCount();
    }

    public int getParameterIndexFromNodeNumber(int nodeNumber) {
        return nodeNumber - tree.getExternalNodeCount();
    }

    public TreeModel getTreeModel() {
        return tree;
    }

    public String[] getNodeAttributeLabel() {
        return new String[]{parameter.getId()};
    }

    public String[] getAttributeForNode(NodeRef node) {
        return new String[]{getAttributeForBranch(node)};
    }

    public double getLogBranchLikelihood(NodeRef node){
        return getLogBranchLikelihood(node.getNumber());
    }

    public double getLogBranchLikelihood(int nodeNum){
        return logBranchLikelihoods[nodeNum];
    }

    public void setLogBranchLikelihood(NodeRef node, double value){
        logBranchLikelihoods[node.getNumber()] = value;
    }

    public double getStoredLogBranchLikelihood(int nodeNum){
        return storedLogBranchLikelihoods[nodeNum];
    }

    public Patterns getMicrosatPattern(){
        return microsatPat;
    }


    /**
     * initialise the internal states using th mid point of the states of two children
     */
    public void initialiseInternalStates(){
        int nodeCount = tree.getNodeCount();
        int maxState = microsatPat.getDataType().getStateCount() - 1;
        int internalNodeCount = tree.getInternalNodeCount();
        Parameter.DefaultBounds bounds = new Parameter.DefaultBounds(maxState, 0, internalNodeCount);
        for(int nodeNum = tree.getExternalNodeCount(); nodeNum < nodeCount; nodeNum++){
            NodeRef node = tree.getNode(nodeNum);
            NodeRef leftChild = tree.getChild(node, 0);
            NodeRef rightChild = tree.getChild(node, 1);
            int leftChildState = getNodeValue(leftChild);
            int rightChildState = getNodeValue(rightChild);
            int nodeValue = (leftChildState+rightChildState)/2;
            parameter.setParameterValueQuietly(nodeNum-tree.getExternalNodeCount(), nodeValue);
        }
        parameter.addBounds(bounds);


    }

}



