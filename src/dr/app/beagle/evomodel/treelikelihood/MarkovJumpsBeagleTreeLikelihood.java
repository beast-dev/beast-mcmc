package dr.app.beagle.evomodel.treelikelihood;

import dr.app.beagle.evomodel.sitemodel.BranchSiteModel;
import dr.app.beagle.evomodel.sitemodel.SiteRateModel;
import dr.app.beagle.evomodel.substmodel.MarkovJumpsSubstitutionModel;
import dr.app.beagle.evomodel.substmodel.SubstitutionModel;
import dr.evolution.alignment.PatternList;
import dr.evolution.datatype.DataType;
import dr.evolution.tree.NodeRef;
import dr.evolution.tree.Tree;
import dr.evomodel.branchratemodel.BranchRateModel;
import dr.evomodel.tree.TreeModel;
import dr.inference.model.Parameter;
import dr.inference.model.Variable;
import dr.inference.loggers.LogColumn;
import dr.inference.loggers.NumberColumn;

import java.util.List;
import java.util.ArrayList;

/**
 * @author Marc Suchard
 * @author Vladimir Minin
 *         <p/>
 *         A base class for implementing Markov chain-induced counting processes (markovjumps) in BEAST using BEAGLE
 *         This work is supported by NSF grant 0856099
 *         <p/>
 *         Minin VN and Suchard MA (2008) Counting labeled transitions in continous-time Markov models of evolution.
 *         Journal of Mathematical Biology, 56, 391-412.
 */
public class MarkovJumpsBeagleTreeLikelihood extends AncestralStateBeagleTreeLikelihood {

    public MarkovJumpsBeagleTreeLikelihood(PatternList patternList, TreeModel treeModel,
                                           BranchSiteModel branchSiteModel, SiteRateModel siteRateModel,
                                           BranchRateModel branchRateModel, boolean useAmbiguities,
                                           PartialsRescalingScheme scalingScheme, DataType dataType, String stateTag,
                                           SubstitutionModel substModel) {

        super(patternList, treeModel, branchSiteModel, siteRateModel, branchRateModel, useAmbiguities,
                scalingScheme, dataType, stateTag, substModel);

        markovjumps = new ArrayList<MarkovJumpsSubstitutionModel>();
        registerParameter = new ArrayList<Parameter>();
        jumpTag = new ArrayList<String>();
        expectedJumps = new ArrayList<double[][]>();
        
        tmpProbabilities = new double[stateCount * stateCount];
        condJumps = new double[stateCount * stateCount];
    }

    public void addRegister(Parameter addRegisterParameter) {

        if (addRegisterParameter.getDimension() != stateCount * stateCount) {
            throw new RuntimeException("Register parameter of wrong dimension");
        }
        addVariable(addRegisterParameter);
        registerParameter.add(addRegisterParameter);
        markovjumps.add(new MarkovJumpsSubstitutionModel(substitutionModel));
        setupRegistration(numRegisters);
        numRegisters++;
        jumpTag.add(addRegisterParameter.getId());
        expectedJumps.add(new double[treeModel.getNodeCount()][patternCount]);
    }

    public String[] getNodeAttributeLabel() {
        String[] old = super.getNodeAttributeLabel();
        String[] rtn = new String[old.length + numRegisters];
        System.arraycopy(old,0,rtn,0,old.length);
        for(int r = 0; r < numRegisters; r++) {
            rtn[old.length+r] = jumpTag.get(r);
        }
        return rtn;
    }

    public String[] getAttributeForNode(Tree tree, NodeRef node) {
        String[] old = super.getAttributeForNode(tree,node);
        String[] rtn = new String[old.length + numRegisters];
        System.arraycopy(old,0,rtn,0,old.length);
        for(int r = 0; r < numRegisters; r++) {
            double[][] thisExpectedJumps = expectedJumps.get(r);
            rtn[old.length + r] = formattedValue(thisExpectedJumps[node.getNumber()]);
        }
        return rtn;
    }

    private static String formattedValue(double[] values) {
        double total = 0;
        for (double summant : values) {
            total += summant;
        }
        return Double.toString(total); // Currently return the sum across sites
    }

    private void setupRegistration(int whichRegistration) {

        double[] registration = registerParameter.get(whichRegistration).getParameterValues();
        markovjumps.get(whichRegistration).setRegistration(registration);
        areStatesRedrawn = false;
    }

    protected void handleVariableChangedEvent(Variable variable, int index, Parameter.ChangeType type) {
        for(int r = 0; r < numRegisters; r++) {
            if (variable == registerParameter.get(r)) {
                setupRegistration(r);
                return;
            }
        }
        super.handleVariableChangedEvent(variable, index, type);
    }

    protected void hookCalculation(Tree tree, NodeRef parentNode, NodeRef childNode,
                                   int[] parentStates, int[] childStates,
                                   double[] inProbabilities) {

        final int childNum = childNode.getNumber();

        double[] probabilities = inProbabilities;
        if (probabilities == null) { // Leaf will call this hook with a null
            getMatrix(childNum, tmpProbabilities);
            probabilities = tmpProbabilities;
        }
        
        final double branchRate = branchRateModel.getBranchRate(tree, childNode);
        final double branchTime = branchRate * (tree.getNodeHeight(parentNode) - tree.getNodeHeight(childNode));
        double[] categoryRates = this.siteRateModel.getCategoryRates();

        if (categoryRates.length > 1) {
            throw new RuntimeException("MarkovJumps only implemented for one rate category");
        }

        final double substTime = categoryRates[0] * branchTime;

        for(int r = 0; r < markovjumps.size(); r++) {
            // Fill condJumps with conditional mean values for this branch
            markovjumps.get(r).computeCondMeanMarkovJumps(substTime,probabilities,condJumps);

            double[][] thisExpectedJumps = expectedJumps.get(r);

            for(int j=0; j<patternCount; j++) { // Pick out values given parent and child states
                thisExpectedJumps[childNum][j] = condJumps[parentStates[j] * stateCount + childStates[j]];
            }
        }
        
    }

    public LogColumn[] getColumns() {
        LogColumn[] allColumns = new LogColumn[patternCount * numRegisters];
        for(int r=0; r<numRegisters; r++) {
            for(int j=0; j<patternCount; j++) {
                allColumns[r*patternCount + j] = new CountColumn(jumpTag.get(r),r,j);
            }
        }
        return allColumns;
    }

    protected class CountColumn extends NumberColumn {
        private int indexSite;
        private int indexRegistration;

        public CountColumn(String label, int r, int j) {
            super(label+"["+(j+1)+"]");
            indexRegistration = r;
            indexSite = j;
        }

        public double getDoubleValue() {
            double total = 0;
            double[][] thisExpectedJumps = expectedJumps.get(indexRegistration);
            for(int i=0; i<treeModel.getNodeCount(); i++) {
                total += thisExpectedJumps[i][indexSite];
            }
            return total;
        }
    }

    private List<MarkovJumpsSubstitutionModel> markovjumps;
    private List<Parameter> registerParameter;
    private List<String> jumpTag;
    private List<double[][]> expectedJumps;
    private double[] tmpProbabilities;
    private double[] condJumps;
    private int numRegisters;
}
