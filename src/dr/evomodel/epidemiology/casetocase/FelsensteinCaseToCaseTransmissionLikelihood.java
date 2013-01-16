package dr.evomodel.epidemiology.casetocase;

import dr.evolution.tree.FlexibleTree;
import dr.evolution.tree.NodeRef;
import dr.evolution.tree.Tree;
import dr.evolution.util.Date;
import dr.evolution.util.Taxon;
import dr.evolution.util.TaxonList;
import dr.evolution.util.Units;
import dr.evomodel.tree.TreeModel;
import dr.inference.loggers.LogColumn;
import dr.inference.loggers.Loggable;
import dr.inference.model.AbstractModelLikelihood;
import dr.inference.model.Model;
import dr.inference.model.Parameter;
import dr.inference.model.Variable;
import dr.util.Author;
import dr.util.Citable;
import dr.util.Citation;
import dr.xml.*;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;

/**
 * A likelihood function for transmission between identified epidemiological cases. The network is not sampled;
 * instead the overall likelihood of the tree based on current values of epidemiological parameters are calculated
 * using a modification of the Felsenstein pruning algorithm. The maximum a posteriori network can then be reconstructed
 * when the phylogeny is sampled.
 *
 * This version assumes that the TMRCA of two sequences represents the time of the transmission that split
 * the lineages. The plan is to relax this in future.
 *
 * Timescale must be in days at present.
 *
 * Python scripts to write XML for it and analyse the posterior set of networks exist or will; contact MH.
 *
 * @author Matthew Hall
 * @author Andrew Rambaut
 */
public class FelsensteinCaseToCaseTransmissionLikelihood extends AbstractModelLikelihood implements Loggable, Citable {

    /* The phylogenetic tree. */

    private TreeModel virusTree;

    /* Matches cases to external nodes */

    private HashMap<AbstractCase, Integer> tipMap;
    private double estimatedLastSampleTime;
    boolean verbose;

    /**
     * The set of cases
     */
    private AbstractCaseSet cases;
    private boolean likelihoodKnown = false;
    private double logLikelihood;

    // PUBLIC STUFF

    // Name

    public static final String FELSENSTEIN_CASE_TO_CASE_TRANSMISSION_LIKELIHOOD =
            "felsensteinCaseToCaseTransmissionLikelihood";

    // Basic constructor.

    public FelsensteinCaseToCaseTransmissionLikelihood(TreeModel virusTree, AbstractCaseSet caseData)
            throws TaxonList.MissingTaxonException {
        this(FELSENSTEIN_CASE_TO_CASE_TRANSMISSION_LIKELIHOOD, virusTree, caseData);
    }

    // Constructor for an instance with a non-default name

    public FelsensteinCaseToCaseTransmissionLikelihood(String name, TreeModel virusTree, AbstractCaseSet caseData)
            throws TaxonList.MissingTaxonException {

        super(name);
        this.virusTree = virusTree;
        addModel(virusTree);
        Date lastSampleDate = getLatestTaxonDate(virusTree);

        /* We assume samples were taken at the end of the date of the last tip */

        estimatedLastSampleTime = lastSampleDate.getTimeValue();
        cases = caseData;
        addModel(cases);
        verbose = false;

        tipMap = new HashMap<AbstractCase, Integer>();

        for(int i=0; i<virusTree.getExternalNodeCount(); i++){
            TreeModel.Node currentExternalNode = (TreeModel.Node)virusTree.getExternalNode(i);
            Taxon currentTaxon = currentExternalNode.taxon;
            for(AbstractCase thisCase : cases.getCases()){
                for(Taxon caseTaxon: thisCase.getAssociatedTaxa()){
                    if(caseTaxon.equals(currentTaxon)){
                        tipMap.put(thisCase, currentExternalNode.getNumber());
                    }
                }
            }
        }
    }

    /* Get the date of the last tip */

    private static Date getLatestTaxonDate(TreeModel tree){
        Date latestDate = new Date(new java.util.Date(-1000000), Units.Type.DAYS);
        for(int i=0;i<tree.getTaxonCount();i++){
            if(tree.getTaxon(i).getDate().after(latestDate)){
                latestDate = tree.getTaxon(i).getDate();
            }
        }
        return latestDate;
    }

    /* Return the integer day on which the given node occurred */

    private int getNodeDay(NodeRef node){
        Double nodeHeight = getHeight(node);
        Double estimatedNodeTime = estimatedLastSampleTime-nodeHeight;
        Date nodeDate = new Date(estimatedNodeTime, Units.Type.DAYS, false);
        /*Since day t begins at time t-1, to find the day of a given decimal time value, take the ceiling.*/
        return (int)Math.ceil(nodeDate.getTimeValue());
    }

    /* return the height of a given node */

    private double getHeight(NodeRef node){
        return virusTree.getNodeHeight(node);
    }

    // **************************************************************
    // ModelListener IMPLEMENTATION
    // **************************************************************


    protected final void handleModelChangedEvent(Model model, Object object, int index) {
        likelihoodKnown = false;
    }


    // **************************************************************
    // VariableListener IMPLEMENTATION
    // **************************************************************


    protected final void handleVariableChangedEvent(Variable variable, int index, Parameter.ChangeType type) {
        likelihoodKnown = false;
    }

    // **************************************************************
    // Model IMPLEMENTATION
    // **************************************************************

    /**
     * Stores the precalculated state
     */
    protected final void storeState() {
    }

    /**
     * Restores the precalculated state.
     */
    protected final void restoreState() {
        likelihoodKnown = false;
    }

    protected final void acceptState() {
    } // nothing to do

    // **************************************************************
    // Likelihood IMPLEMENTATION
    // **************************************************************

    public final Model getModel() {
        return this;
    }

    public final TreeModel getTree(){
        return virusTree;
    }

    public final double getLogLikelihood() {
        if (!likelihoodKnown) {
            logLikelihood = calculateLogLikelihood();
            likelihoodKnown = true;
        }
        return logLikelihood;
    }

    public final void makeDirty() {
        likelihoodKnown = false;
    }

    /**
     * Calculates the log likelihood of the tree under current values of the epidemiological parameters
     */
    public double calculateLogLikelihood() {
        double totalLikelihood = 0;
        NodeRef root = virusTree.getRoot();
        for(AbstractCase thisCase: cases.getCases()){
            totalLikelihood += L(root,thisCase);
        }
        return Math.log(totalLikelihood);
    }

    /* Check whether the tree will admit a particular painting for a particular node. This is true if a
     * descendant external node is the one corresponding to the case */

    private boolean treeWillAdmitCaseHere (NodeRef node, AbstractCase thisCase){
        if(virusTree.isExternal(node)){
            return tipMap.get(thisCase)==node.getNumber();
        } else {
            boolean check1 = treeWillAdmitCaseHere(virusTree.getChild(node,0), thisCase);
            boolean check2 = treeWillAdmitCaseHere(virusTree.getChild(node,1), thisCase);
            return check1 || check2;
        }
    }

    /* The likelihood of observing a case i at node alpha*/

    private double L(NodeRef alpha, AbstractCase i){
        if(virusTree.isExternal(alpha)){
            return tipMap.get(i)==alpha.getNumber() ? 1:0;
        } else {
            double totalLikelihood = 0;
            NodeRef[] children = getChildren(alpha);
            for(AbstractCase case_x: cases.getCases()){
                for(AbstractCase case_y: cases.getCases()){
                    double P_term = P(i, case_x, case_y, alpha);
                    double L_term = 0;
                    if(P_term!=0){
                        L_term = L(children[0],case_x)*L(children[1],case_y);
                    }
                    totalLikelihood = totalLikelihood + L_term*P_term;
                }
            }
            return totalLikelihood;
        }
    }

    /* The maximum likelihood version of the above for ancestral state reconstruction. */

    private double L_max(NodeRef alpha, AbstractCase i){
        if(virusTree.isExternal(alpha)){
            return tipMap.get(i)==alpha.getNumber() ? 1:0;
        } else {
            double maxLikelihood = 0;
            NodeRef[] children = getChildren(alpha);
            for(AbstractCase case_x: cases.getCases()){
                for(AbstractCase case_y: cases.getCases()){
                    double P_term = P(i, case_x, case_y, alpha);
                    double L_term = 0;
                    if(P_term!=0){
                        L_term = L_max(children[0],case_x)*L_max(children[1],case_y);
                    }
                    if(L_term*P_term > maxLikelihood){
                        maxLikelihood = L_term*P_term;
                    }
                }
            }
            return maxLikelihood;
        }
    }

    /* The probability that you observe an a at alpha given a b at one child and a c at the other*/

    private double P(AbstractCase a, AbstractCase b, AbstractCase c, NodeRef alpha){
        // Internal nodes only
        if(virusTree.isExternal(alpha)){
            throw new RuntimeException("alpha must be external");
        }
        NodeRef[] children = getChildren(alpha);
        int alphaDay = getNodeDay(alpha);
        int[] childDays = getChildDays(alpha);
        //Tree must admit the cases in these positions
        if(!treeWillAdmitCaseHere(alpha, a)
                || !treeWillAdmitCaseHere(children[0],b)
                || !treeWillAdmitCaseHere(children[1],c)) {
            return 0;
        }
        if((a!=b && a!=c) || (b==c)){
            return 0;
        }
        if (a==b){
            double c_term = cases.transmissionBranchLikelihood(a, c, alphaDay, childDays[1]);
            double a_term = cases.noTransmissionBranchLikelihood(a, alphaDay)
                    /cases.noTransmissionBranchLikelihood(a, childDays[0]);
            if(!Double.isNaN(a_term)){
                return c_term*a_term;
            } else {
                return 0;
            }
        } else {
            double b_term = cases.transmissionBranchLikelihood(a, b, alphaDay, childDays[0]);
            double a_term = cases.noTransmissionBranchLikelihood(a, alphaDay)
                    /cases.noTransmissionBranchLikelihood(a, childDays[1]);
            if(!Double.isNaN(a_term)){
                return b_term*a_term;
            } else {
                return 0;
            }
        }
    }

    /* The method to get the maximum likelihood painting and thus network. I'm not yet sure that this procedure
     actually results in a painting, and need to check hard for counterexamples. */

    private HashMap<AbstractCase, AbstractCase> getMaximumLikelihoodNetwork(){
        HashMap<NodeRef, AbstractCase> maxLikelihoodPainting = new HashMap<NodeRef, AbstractCase>();
        NodeRef root = virusTree.getRoot();
        AbstractCase rootPainting = null;
        double maxRootPaintingLikelihood = 0;
        for(AbstractCase thisCase: cases.getCases()){
            double likelihood = L_max(root,thisCase);
            if(likelihood>maxRootPaintingLikelihood){
                maxRootPaintingLikelihood = likelihood;
                rootPainting = thisCase;
            }
        }
        maxLikelihoodPainting.put(root,rootPainting);
        fillUp(root,maxLikelihoodPainting);
        HashMap<AbstractCase, AbstractCase> network = new HashMap<AbstractCase, AbstractCase>();
        for(AbstractCase thisCase: cases.getCases()){
            network.put(thisCase, getInfector(thisCase,maxLikelihoodPainting));
        }
        return network;
    }

    private void fillUp(NodeRef node, HashMap<NodeRef, AbstractCase> painting){
        if(!virusTree.isExternal(node)){
            AbstractCase[] bestChildren = bestChildren(node, painting.get(node));
            for(int i=0; i<bestChildren.length; i++){
                NodeRef child = virusTree.getChild(node,i);
                painting.put(child,bestChildren[i]);
                fillUp(child, painting);
            }
        }
    }

    private AbstractCase[] bestChildren(NodeRef node, AbstractCase nodePainting){
        double maxLikelihood = 0;
        AbstractCase[] bestPair = new AbstractCase[2];
        for(AbstractCase left: cases.getCases()){
            for(AbstractCase right: cases.getCases()){
                double P_term = P(nodePainting, left, right, node);
                double L_term = 0;
                if(P_term>0){
                    L_term = L_max(virusTree.getChild(node,0),left)*L_max(virusTree.getChild(node,1),right);
                }
                if(P_term*L_term>maxLikelihood){
                    bestPair[0] = left;
                    bestPair[1] = right;
                }
            }
        }
        System.out.println("Best pair for node "+node.toString()+" found.");
        return bestPair;
    }

    public AbstractCase getInfector(AbstractCase thisCase, HashMap<NodeRef, AbstractCase> branchMap){
        NodeRef tip = virusTree.getNode(tipMap.get(thisCase));
        return getInfector(tip, branchMap);
    }

    public AbstractCase getInfector(NodeRef node, HashMap<NodeRef, AbstractCase> branchMap){
        if(virusTree.isRoot(node) || node.getNumber()==virusTree.getRoot().getNumber()){
            return null;
        } else {
            AbstractCase nodeCase = branchMap.get(node);
            if(branchMap.get(node)!=nodeCase){
                return branchMap.get(node);
            } else {
                return getInfector(virusTree.getParent(node), branchMap);
            }
        }
    }


    private NodeRef[] getChildren(NodeRef node){
        NodeRef[] children = new NodeRef[virusTree.getChildCount(node)];
        for(int i=0; i<virusTree.getChildCount(node); i++){
            children[i] = virusTree.getChild(node,i);
        }
        return children;
    }

    private int[] getChildDays(NodeRef node){
        int[] childDays = new int[virusTree.getChildCount(node)];
        for(int i=0; i<virusTree.getChildCount(node); i++){
            childDays[i] = getNodeDay(virusTree.getChild(node,i));
        }
        return childDays;
    }


    //************************************************************************
    // AbstractXMLObjectParser implementation
    //************************************************************************

    public static XMLObjectParser PARSER = new AbstractXMLObjectParser() {
        public static final String VIRUS_TREE = "virusTree";

        public String getParserName() {
            return FELSENSTEIN_CASE_TO_CASE_TRANSMISSION_LIKELIHOOD;
        }

        public Object parseXMLObject(XMLObject xo) throws XMLParseException {

            FlexibleTree flexTree = (FlexibleTree) xo.getElementFirstChild(VIRUS_TREE);

            TreeModel virusTree = new TreeModel(flexTree);

            AbstractCaseSet caseSet = (AbstractCaseSet) xo.getChild(AbstractCaseSet.class);

            FelsensteinCaseToCaseTransmissionLikelihood likelihood;

            try {
                likelihood = new FelsensteinCaseToCaseTransmissionLikelihood(virusTree, caseSet);
            } catch (TaxonList.MissingTaxonException e) {
                throw new XMLParseException(e.toString());
            }

            return likelihood;
        }

        public String getParserDescription() {
            return "This element represents a likelihood function for case to case transmission.";
        }

        public Class getReturnType() {
            return FelsensteinCaseToCaseTransmissionLikelihood.class;
        }

        public XMLSyntaxRule[] getSyntaxRules() {
            return rules;
        }

        private final XMLSyntaxRule[] rules = {
                new ElementRule("virusTree", Tree.class, "The (currently fixed) tree"),
                new ElementRule(AbstractCaseSet.class, "The set of cases")
        };
    };

    //************************************************************************
    // Loggable implementation
    //************************************************************************


    public LogColumn[] getColumns(){

        LogColumn[] columns = new LogColumn[cases.size()];

        final HashMap<AbstractCase,AbstractCase> reconstructedNetwork = getMaximumLikelihoodNetwork();

        for(int i=0; i< cases.size(); i++){
            final AbstractCase infected = cases.getCase(i);
            columns[i] = new LogColumn.Abstract(infected.toString()){
                @Override
                protected String getFormattedValue() {
                    if(reconstructedNetwork.get(infected)==null){
                        return "Start";
                    } else {
                        return reconstructedNetwork.get(infected).toString();
                    }
                }
            };
        }

        return columns;

    }

    //************************************************************************
    // Citable implementation
    //************************************************************************


    public List<Citation> getCitations() {
        List<Citation> citations = new ArrayList<Citation>();
        citations.add(new Citation(new Author[]{new Author("M", "Hall"), new Author("A", "Rambaut")},
                Citation.Status.IN_PREPARATION));
        return citations;
    }

}
