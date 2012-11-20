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
import dr.math.MathUtils;
import dr.util.Author;
import dr.util.Citable;
import dr.util.Citation;
import dr.xml.*;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.List;

/**
 * A likelihood function for transmission between identified epidemiological cases
 *
 * Currently works only for fixed trees and estimates the network and epidemiological parameters. Timescale
 * must be in days. Python scripts to write XML for it and analyse the posterior set of networks exist; contact MH.
 *
 * @author Matthew Hall
 * @author Andrew Rambaut
 * @version $Id: $
 */
public class CaseToCaseTransmissionLikelihood extends AbstractModelLikelihood implements Loggable, Citable {

    /* The phylogenetic tree. */

    private TreeModel virusTree;

    /* Mapping of cases to branches on the tree; old version is stored before operators are applied */

    private AbstractCase[] branchMap;
    private AbstractCase[] storedBranchMap;

    /* The log likelihood of the subtree from the parent branch of the referred-to node downwards; old version is stored
    before operators are applied.*/

    private double[] subTreeLogLikelihoods;
    private double[] storedSubTreeLogLikelihoods;

    /* Whether operations have required the recalculation of the log likelihoods of the subtree from the parent branch
    of the referred-to node downwards.
     */

    private boolean[] subTreeRecalculationNeeded;

    /* Matches cases to external nodes */

    private HashMap<AbstractCase, NodeRef> tipMap;
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

    public static final String CASE_TO_CASE_TRANSMISSION_LIKELIHOOD = "caseToCaseTransmissionLikelihood";

    // Basic constructor.

    public CaseToCaseTransmissionLikelihood(TreeModel virusTree, AbstractCaseSet caseData,
                                            String startingNetworkFileName)
            throws TaxonList.MissingTaxonException {
        this(CASE_TO_CASE_TRANSMISSION_LIKELIHOOD, virusTree, caseData, startingNetworkFileName);
    }


    // Constructor for an instance with a non-default name

    public CaseToCaseTransmissionLikelihood(String name, TreeModel virusTree, AbstractCaseSet caseData, String
            startingNetworkFileName)
            throws TaxonList.MissingTaxonException {

        super(name);
        this.virusTree = virusTree;
        addModel(virusTree);
        Date lastSampleDate = getLatestTaxonDate(virusTree);

        /* We assume samples were taken at midday on the date of the last tip */

        estimatedLastSampleTime = lastSampleDate.getTimeValue()+0.5;
        cases = caseData;
        addModel(cases);
        verbose = false;

        tipMap = new HashMap<AbstractCase, NodeRef>();

        subTreeRecalculationNeeded = new boolean[virusTree.getNodeCount()];
        Arrays.fill(subTreeRecalculationNeeded, true);
        subTreeLogLikelihoods = new double[virusTree.getNodeCount()];

        if(startingNetworkFileName==null){
            branchMap = paintRandomNetwork();
        } else {
            branchMap = paintSpecificNetwork(startingNetworkFileName);
        }

        for(int i=0; i<virusTree.getExternalNodeCount(); i++){
            TreeModel.Node currentExternalNode = (TreeModel.Node)virusTree.getExternalNode(i);
            Taxon currentTaxon = currentExternalNode.taxon;
            for(AbstractCase thisCase : cases.getCases()){
                for(Taxon caseTaxon: thisCase.getAssociatedTaxa()){
                    if(caseTaxon.equals(currentTaxon)){
                        tipMap.put(thisCase, currentExternalNode);
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

    // **************************************************************
    // ModelListener IMPLEMENTATION
    // **************************************************************


    protected final void handleModelChangedEvent(Model model, Object object, int index) {
        Arrays.fill(subTreeRecalculationNeeded, true);
        likelihoodKnown = false;

    }


    // **************************************************************
    // VariableListener IMPLEMENTATION
    // **************************************************************


    protected final void handleVariableChangedEvent(Variable variable, int index, Parameter.ChangeType type) {
        Arrays.fill(subTreeRecalculationNeeded, true);
        likelihoodKnown = false;
    }

    // **************************************************************
    // Model IMPLEMENTATION
    // **************************************************************

    /**
     * Stores the precalculated state (in this case the node labels and subtree likelihoods)
     */
    protected final void storeState() {
        storedBranchMap = Arrays.copyOf(branchMap, branchMap.length);
        storedSubTreeLogLikelihoods = Arrays.copyOf(subTreeLogLikelihoods, subTreeLogLikelihoods.length);
    }

    /**
     * Restores the precalculated state.
     */
    protected final void restoreState() {
        branchMap = storedBranchMap;
        subTreeLogLikelihoods = storedSubTreeLogLikelihoods;
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

    public final AbstractCase[] getBranchMap(){
        return branchMap;
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
        Arrays.fill(subTreeRecalculationNeeded, true);
        likelihoodKnown = false;
    }

    /**
     * Calculates the log likelihood of this set of node labels given the tree.
     */
    public double calculateLogLikelihood() {
        return calculateLogLikelihood(branchMap);
    }

    public double calculateLogLikelihood(AbstractCase[] map){
        NodeRef root = virusTree.getRoot();
        return calculateNodeTransmissionLogLikelihood(root, Integer.MIN_VALUE, map);
    }

    /* Return the integer day on which the given node occurred */

    private int getNodeDay(NodeRef node){
        Double nodeHeight = getHeight(node);
        Double estimatedNodeTime = estimatedLastSampleTime-nodeHeight;
        Date nodeDate = new Date(estimatedNodeTime, Units.Type.DAYS, false);
        /*Since day t begins at time t-1, to find the day of a given decimal time value, take the ceiling.*/
        return (int)Math.ceil(nodeDate.getTimeValue());
    }

    /**
     * Given a node, calculates the log likelihood of its parent branch and then goes down the tree and calculates the
     * log likelihoods of lower branches.
     */

    private double calculateNodeTransmissionLogLikelihood(NodeRef node, int parentDay, AbstractCase[]
            currentBranchMap) {

        double logLikelihood=0;
        AbstractCase nodeCase = currentBranchMap[node.getNumber()];
        Integer nodeDay = getNodeDay(node);
        /* Likelihood of the branch leading to the node.*/
        if(subTreeRecalculationNeeded[node.getNumber()]){
            if(virusTree.isRoot(node)){
                /*Likelihood of root node case being infectious by the time of the root.*/
                logLikelihood = logLikelihood + Math.log(cases.rootBranchLikelihood(nodeCase, nodeDay));
            } else {
                /*Likelihood of this case being infected at the time of the parent node and infectious by the time
                of the child node.*/
                NodeRef parent = virusTree.getParent(node);
                logLikelihood = logLikelihood + Math.log(cases.branchLikelihood(currentBranchMap[parent.getNumber()],
                        nodeCase, parentDay, nodeDay));
            }
            /* If this isn't an external node, get the log likelihood of lower branches */
            if(!virusTree.isExternal(node)){
                for(int i=0; i<virusTree.getChildCount(node); i++){
                    logLikelihood = logLikelihood + calculateNodeTransmissionLogLikelihood(virusTree.getChild(node,i),
                            nodeDay, currentBranchMap);
                }
            }
            subTreeLogLikelihoods[node.getNumber()]=logLikelihood;
            subTreeRecalculationNeeded[node.getNumber()]=false;
        }
        else{
            logLikelihood = subTreeLogLikelihoods[node.getNumber()];
        }
        return logLikelihood;
    }

    public boolean[] getRecalculationArray(){
        return subTreeRecalculationNeeded;
    }

    private double getHeight(NodeRef node){
        return virusTree.getNodeHeight(node);
    }

    /* Return the case whose infection resulted from this node (if internal) */

    private AbstractCase getInfected(NodeRef node){
        if(virusTree.isExternal(node)){
            throw new RuntimeException("Node is external");
        } else {
            for(int i=0;i<2;i++){
                NodeRef child = virusTree.getChild(node, i);
                if(branchMap[child.getNumber()]!=branchMap[node.getNumber()]){
                    return branchMap[child.getNumber()];
                }
            }
        }
        throw new RuntimeException("Node appears to have two children with its infector");
    }

    /* Return the case which infected this case */

    public AbstractCase getInfector(AbstractCase thisCase){
        NodeRef tip = tipMap.get(thisCase);
        return getInfector(tip);
    }

    /* Return the case which was the infector in the infection event represented by this node */

    public AbstractCase getInfector(NodeRef node){
        if(virusTree.isRoot(node) || node.getNumber()==virusTree.getRoot().getNumber()){
            return null;
        } else {
            AbstractCase nodeCase = branchMap[node.getNumber()];
            if(branchMap[virusTree.getParent(node).getNumber()]!=nodeCase){
                return branchMap[virusTree.getParent(node).getNumber()];
            } else {
                return getInfector(virusTree.getParent(node));
            }
        }
    }

    /* Return the case (painting) of the parent of this node */

    public AbstractCase getParentCase(NodeRef node){
        return branchMap[virusTree.getParent(node).getNumber()];
    }

    /*Given a new tree with no labels, associates each of the terminal branches with the relevant case and then
    * generates a random painting of the rest of the tree to start off with. If checkNonZero is true in
    * randomlyPaintNode then the network will be checked to prohibit links with zero (or rounded to zero)
    * likelihood first.*/

    private AbstractCase[] paintRandomNetwork(){
        boolean gotOne = false;
        AbstractCase[] map = null;
        int tries = 1;
        System.out.println("Generating a random starting painting of the tree (checking nonzero likelihood for all " +
                "branches and repeating until a start with nonzero likelihood is found)");
        System.out.print("Attempt: ");
        while(!gotOne){
            System.out.print(tries + "...");
            map = prepareExternalNodeMap(new AbstractCase[virusTree.getNodeCount()]);
            TreeModel.Node root = (TreeModel.Node)virusTree.getRoot();
            randomlyPaintNode(root, map, true);
            if(calculateLogLikelihood(map)!=Double.NEGATIVE_INFINITY){
                gotOne = true;
                System.out.print("found.");
            }
            tries++;
        }
        System.out.println();
        return map;
    }

    /* Populates the branch map for external nodes */

    private AbstractCase[] prepareExternalNodeMap(AbstractCase[] map){
        for(int i=0; i<virusTree.getExternalNodeCount(); i++){
            TreeModel.Node currentExternalNode = (TreeModel.Node)virusTree.getExternalNode(i);
            Taxon currentTaxon = currentExternalNode.taxon;
            for(AbstractCase thisCase : cases.getCases()){
                for(Taxon caseTaxon: thisCase.getAssociatedTaxa()){
                    if(caseTaxon.equals(currentTaxon)){
                        map[currentExternalNode.getNumber()]=thisCase;
                    }
                }
            }
        }
        return map;
    }

/*  The CSV file should have a header line, and then a line for each node with its id in the first column and the id
    of its parent in the second. The node with no parent has "Start" in the second column.*/

    private AbstractCase[] paintSpecificNetwork(String networkFileName){
        System.out.println("Using specified starting network.");
        try{
            BufferedReader reader = new BufferedReader (new FileReader(networkFileName));
            HashMap<AbstractCase, AbstractCase> specificParentMap = new HashMap<AbstractCase, AbstractCase>();
            // skip header line
            reader.readLine();
            String currentLine = reader.readLine();
            while(currentLine!=null){
                currentLine = currentLine.replace("\"", "");
                String[] splitLine = currentLine.split("\\,");
                if(!splitLine[1].equals("Start")){
                    specificParentMap.put(cases.getCase(splitLine[0]), cases.getCase(splitLine[1]));
                } else {
                    specificParentMap.put(cases.getCase(splitLine[0]),null);
                }
                currentLine = reader.readLine();
            }
            reader.close();
            return paintSpecificNetwork(specificParentMap);
        } catch(IOException e){
            System.out.println("Problem reading file (IOException)");
            return null;
        }
    }


    /* Takes a HashMap referring each case to its parent, and tries to paint the tree with it */

    private AbstractCase[] paintSpecificNetwork(HashMap<AbstractCase, AbstractCase> map){
        AbstractCase[] branchArray = new AbstractCase[virusTree.getNodeCount()];
        branchArray = prepareExternalNodeMap(branchArray);
        TreeModel.Node root = (TreeModel.Node)virusTree.getRoot();
        specificallyPaintNode(root, branchArray, map);
        return branchArray;
    }

    /* Paints a phylogenetic tree node and its children according to a specified map of child to parent cases */

    private AbstractCase specificallyPaintNode(TreeModel.Node node, AbstractCase[] map,
                                               HashMap<AbstractCase,AbstractCase> parents){
        if(node.isExternal()){
            return map[node.getNumber()];
        } else {
            AbstractCase[] childPaintings = new AbstractCase[2];
            for(int i=0; i<node.getChildCount(); i++){
                childPaintings[i] = specificallyPaintNode(node.getChild(i), map, parents);
            }
            if(parents.get(childPaintings[1])==childPaintings[0]){
                map[node.getNumber()]=childPaintings[0];
            } else if(parents.get(childPaintings[0])==childPaintings[1]){
                map[node.getNumber()]=childPaintings[1];
            } else {
                throw new RuntimeException("This network does not appear to be compatible with the tree");
            }
            return map[node.getNumber()];
        }
    }

    /* Paints a phylogenetic tree with a random compatible painting; if checkNonZero is true, make sure all branch
    likelihoods are nonzero in the process (this sometimes still results in a zero likelihood for the whole tree, but
    is much less likely to).
     */

    private AbstractCase randomlyPaintNode(TreeModel.Node node, AbstractCase[] map, boolean checkNonZero){
        if(node.isExternal()){
            return map[node.getNumber()];
        } else {
            AbstractCase[] choices = new AbstractCase[2];
            for(int i=0; i<node.getChildCount(); i++){
                if((map[node.getChild(i).getNumber()]==null)){
                    choices[i] = randomlyPaintNode(node.getChild(i),map,checkNonZero);
                } else {
                    choices[i] = map[node.getChild(i).getNumber()];
                }
            }
            int randomSelection = MathUtils.nextInt(2);
            AbstractCase decision;
            if(checkNonZero){
                Double[] branchLogLs = new Double[2];
                for(int i=0; i<2; i++){
                    branchLogLs[i]= cases.branchLogLikelihood(choices[i], choices[1-i], getNodeDay(node),
                            getNodeDay(node.getChild(1-i)));
                }
                if(branchLogLs[0]==Double.NEGATIVE_INFINITY && branchLogLs[1]==Double.NEGATIVE_INFINITY){
                    throw new RuntimeException("Both branch possibilities have zero likelihood: "
                            +node.toString()+", cases " + choices[0].getName() + " and " + choices[1].getName() + ".");
                } else if(branchLogLs[0]==Double.NEGATIVE_INFINITY || branchLogLs[1]==Double.NEGATIVE_INFINITY){
                    if(branchLogLs[0]==Double.NEGATIVE_INFINITY){
                        decision = choices[1];
                    } else {
                        decision = choices[0];
                    }
                } else {
                    decision = choices[randomSelection];
                }
            } else {
                decision = choices[randomSelection];
            }
            map[node.getNumber()]=decision;
            return decision;
        }
    }

    //************************************************************************
    // AbstractXMLObjectParser implementation
    //************************************************************************

    public static XMLObjectParser PARSER = new AbstractXMLObjectParser() {
        public static final String VIRUS_TREE = "virusTree";
        public static final String STARTING_NETWORK = "startingNetwork";

        public String getParserName() {
            return CASE_TO_CASE_TRANSMISSION_LIKELIHOOD;
        }

        public Object parseXMLObject(XMLObject xo) throws XMLParseException {

            FlexibleTree flexTree = (FlexibleTree) xo.getElementFirstChild(VIRUS_TREE);

            String startingNetworkFileName=null;

            if(xo.hasChildNamed(STARTING_NETWORK)){
                startingNetworkFileName = (String) xo.getElementFirstChild(STARTING_NETWORK);
            }

            TreeModel virusTree = new TreeModel(flexTree);

            AbstractCaseSet caseSet = (AbstractCaseSet) xo.getChild(AbstractCaseSet.class);

            CaseToCaseTransmissionLikelihood likelihood;

            try {
                likelihood = new CaseToCaseTransmissionLikelihood(virusTree, caseSet, startingNetworkFileName);
            } catch (TaxonList.MissingTaxonException e) {
                throw new XMLParseException(e.toString());
            }

            return likelihood;
        }

        public String getParserDescription() {
            return "This element represents a likelihood function for case to case transmission.";
        }

        public Class getReturnType() {
            return CaseToCaseTransmissionLikelihood.class;
        }

        public XMLSyntaxRule[] getSyntaxRules() {
            return rules;
        }

        private final XMLSyntaxRule[] rules = {
                new ElementRule("virusTree", Tree.class, "The (currently fixed) tree"),
                new ElementRule(AbstractCaseSet.class, "The set of cases"),
                new ElementRule("startingNetwork", String.class, "A CSV file containing a specified starting network",
                        true)
        };
    };


    //************************************************************************
    // Loggable implementation
    //************************************************************************

    public LogColumn[] getColumns(){

        LogColumn[] columns = new LogColumn[cases.size()];

        for(int i=0; i< cases.size(); i++){
            final AbstractCase infected = cases.getCase(i);
            columns[i] = new LogColumn.Abstract(infected.toString()){
                @Override
                protected String getFormattedValue() {
                    if(getInfector(infected)==null){
                        return "Start";
                    } else {
                        return getInfector(infected).toString();
                    }
                }
            };
        }

        return columns;

    }

    public List<Citation> getCitations() {
        List<Citation> citations = new ArrayList<Citation>();
        citations.add(new Citation(new Author[]{new Author("M", "Hall"), new Author("A", "Rambaut")},
                Citation.Status.IN_PREPARATION));
        return citations;
    }
}
