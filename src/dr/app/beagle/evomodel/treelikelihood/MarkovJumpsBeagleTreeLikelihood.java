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
                                           SubstitutionModel substModel, Parameter registerMatrixParameter,
                                           String jumpTag) {

        super(patternList, treeModel, branchSiteModel, siteRateModel, branchRateModel, useAmbiguities,
                scalingScheme, dataType, stateTag, substModel);

        if (registerMatrixParameter.getDimension() != stateCount * stateCount) {
            throw new RuntimeException("Registration parameter of wrong dimension");
        }

        markovjumps = new MarkovJumpsSubstitutionModel(substModel);        

        this.registerMatrixParameter = registerMatrixParameter;
        addVariable(registerMatrixParameter);
        setupRegistration();

        this.jumpTag = jumpTag;

        expectedJumps = new double[treeModel.getNodeCount()][patternCount];
        tmpProbabilities = new double[stateCount * stateCount];
        condJumps = new double[stateCount * stateCount];
    }

    private String[] copyStringPtrsAndAppend(String[] values, String append) {
        String[] rtn = new String[values.length + 1];
        System.arraycopy(values, 0, rtn, 0, values.length);
        rtn[values.length] = append;
        return rtn;
    }

    public String[] getNodeAttributeLabel() {
        return copyStringPtrsAndAppend(super.getNodeAttributeLabel(), jumpTag);
    }

    public String[] getAttributeForNode(Tree tree, NodeRef node) {
        return copyStringPtrsAndAppend(super.getAttributeForNode(tree, node),
                formattedValue(expectedJumps[node.getNumber()])
        );
    }

    private static String formattedValue(double[] values) {
        double total = 0;
        for (double summant : values) {
            total += summant;
        }
        return Double.toString(total); // Currently return the sum across sites
    }

    private void setupRegistration() {

        double[] registration = registerMatrixParameter.getParameterValues();
        for (int i = 0; i < stateCount; i++) {
            registration[i * stateCount + i] = 0;  // Ensure that the diagonals are zero
        }
        markovjumps.setRegistration(registration);
        areStatesRedrawn = false;
    }

    protected void handleVariableChangedEvent(Variable variable, int index, Parameter.ChangeType type) {
        if (variable == registerMatrixParameter) {
            setupRegistration();
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

        // Fill condJumps with conditional mean values for this branch
        markovjumps.computeCondMeanMarkovJumps(branchTime,probabilities,condJumps);

        for(int j=0; j<patternCount; j++) { // Pick out values given parent and child states
            expectedJumps[childNum][j] = condJumps[parentStates[j] * stateCount + childStates[j]];
        }
        
    }

    private MarkovJumpsSubstitutionModel markovjumps;
    private Parameter registerMatrixParameter;
    private String jumpTag;
    private double[][] expectedJumps;
    private double[] tmpProbabilities;
    private double[] condJumps;
}
