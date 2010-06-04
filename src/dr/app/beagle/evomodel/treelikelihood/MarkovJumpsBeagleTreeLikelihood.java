package dr.app.beagle.evomodel.treelikelihood;

import dr.app.beagle.evomodel.sitemodel.BranchSiteModel;
import dr.app.beagle.evomodel.sitemodel.SiteRateModel;
import dr.app.beagle.evomodel.substmodel.MarkovJumpsSubstitutionModel;
import dr.app.beagle.evomodel.substmodel.SubstitutionModel;
import dr.app.beagle.evomodel.substmodel.UniformizedSubstitutionModel;
import dr.evolution.alignment.PatternList;
import dr.evolution.datatype.DataType;
import dr.evolution.tree.NodeRef;
import dr.evolution.tree.Tree;
import dr.evolution.tree.TreeTrait;
import dr.evomodel.branchratemodel.BranchRateModel;
import dr.evomodel.tree.TreeModel;
import dr.inference.markovjumps.MarkovJumpsRegisterAcceptor;
import dr.inference.model.Parameter;
import dr.inference.model.Variable;
import dr.inference.loggers.LogColumn;
import dr.inference.loggers.NumberColumn;
import dr.inference.markovjumps.MarkovJumpsType;

import java.util.List;
import java.util.ArrayList;
import java.util.Arrays;

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
public class MarkovJumpsBeagleTreeLikelihood extends AncestralStateBeagleTreeLikelihood
        implements MarkovJumpsRegisterAcceptor {

    public MarkovJumpsBeagleTreeLikelihood(PatternList patternList, TreeModel treeModel,
                                           BranchSiteModel branchSiteModel, SiteRateModel siteRateModel,
                                           BranchRateModel branchRateModel, boolean useAmbiguities,
                                           PartialsRescalingScheme scalingScheme, DataType dataType, String stateTag,
                                           SubstitutionModel substModel,
                                           boolean useMAP,
                                           boolean returnMarginalLikelihood,
                                           boolean useUniformization,
                                           boolean reportUnconditionedColumns,
                                           int nSimulants) {

        super(patternList, treeModel, branchSiteModel, siteRateModel, branchRateModel, useAmbiguities,
                scalingScheme, dataType, stateTag, substModel, useMAP, returnMarginalLikelihood);

        this.useUniformization = useUniformization;
        this.reportUnconditionedColumns = reportUnconditionedColumns;
        this.nSimulants = nSimulants;

        markovjumps = new ArrayList<MarkovJumpsSubstitutionModel>();
        registerParameter = new ArrayList<Parameter>();
        jumpTag = new ArrayList<String>();
        expectedJumps = new ArrayList<double[][]>();

        tmpProbabilities = new double[stateCount * stateCount * categoryCount];
        condJumps = new double[categoryCount][stateCount * stateCount];
    }

    public void addRegister(Parameter addRegisterParameter,
                            MarkovJumpsType type,
                            boolean scaleByTime) {

        if ((type == MarkovJumpsType.COUNTS &&
                addRegisterParameter.getDimension() != stateCount * stateCount) ||
                (type == MarkovJumpsType.REWARDS &&
                        addRegisterParameter.getDimension() != stateCount)
                ) {
            throw new RuntimeException("Register parameter of wrong dimension");
        }
        addVariable(addRegisterParameter);
        registerParameter.add(addRegisterParameter);
        if (useUniformization) {
            markovjumps.add(new UniformizedSubstitutionModel(substitutionModel, type, nSimulants));
        } else {
            markovjumps.add(new MarkovJumpsSubstitutionModel(substitutionModel, type));
        }
        setupRegistration(numRegisters);

        final String tag = addRegisterParameter.getId();

        jumpTag.add(tag);
        expectedJumps.add(new double[treeModel.getNodeCount()][patternCount]);

        boolean[] oldScaleByTime = this.scaleByTime;
        int oldScaleByTimeLength = (oldScaleByTime == null ? 0 : oldScaleByTime.length);
        this.scaleByTime = new boolean[oldScaleByTimeLength + 1];
        if (oldScaleByTimeLength > 0) {
            System.arraycopy(oldScaleByTime, 0, this.scaleByTime, 0, oldScaleByTimeLength);
        }
        this.scaleByTime[oldScaleByTimeLength] = scaleByTime;

        treeTraits.addTrait("dwellTimes", new TreeTrait<double[]>() {
            public String getTraitName() {
                return tag;
            }

            public Intent getIntent() {
                return Intent.NODE;
            }

            public Class getTraitClass() {
                return double[].class;
            }

            public int getDimension() {
                return 1;
            }

            public double[][] getTrait(Tree tree, NodeRef node) {
                return new double[][] { getDwellTimesForNodeAndPattern(tree, node, 0) };
            }

            public String[] getTraitString(Tree tree, NodeRef node) {
                return new String[0];
            }
        });

        numRegisters++;
    }

    public double[] getDwellTimesForNodeAndPattern(Tree tree, NodeRef node, int pattern) {
        double[] rtn = new double[numRegisters];
        for (int r = 0; r < numRegisters; r++) {
            double[] mjs = getMarkovJumpsForNodeAndRegister(tree, node, r); 
            rtn[r] = mjs[pattern];
        }
        return rtn;
    }


    public double[] getMarkovJumpsForNodeAndRegister(Tree tree, NodeRef node, int whichRegister) {
        return getMarkovJumpsForRegister(tree, whichRegister)[node.getNumber()];
    }

    public double[][] getMarkovJumpsForNode(Tree tree, NodeRef node) {
        double[][] rtn = new double[numRegisters][];
        for (int r = 0; r < numRegisters; r++) {
            rtn[r] = getMarkovJumpsForNodeAndRegister(tree, node, r);
        }
        return rtn;
    }

    public double[][] getMarkovJumpsForRegister(Tree tree, int whichRegister) {
        if (tree != treeModel) {
            throw new RuntimeException("Must call with internal tree");
        }
        if (!areStatesRedrawn) {
            redrawAncestralStates();
        }
        return expectedJumps.get(whichRegister);
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
        for (int r = 0; r < numRegisters; r++) {
            if (variable == registerParameter.get(r)) {
                setupRegistration(r);
                return;
            }
        }
        super.handleVariableChangedEvent(variable, index, type);
    }

    protected void hookCalculation(Tree tree, NodeRef parentNode, NodeRef childNode,
                                   int[] parentStates, int[] childStates,
                                   double[] inProbabilities, int[] rateCategory) {

        final int childNum = childNode.getNumber();

        double[] probabilities = inProbabilities;
        if (probabilities == null) { // Leaf will call this hook with a null
            getMatrix(childNum, tmpProbabilities);
            probabilities = tmpProbabilities;
        }

        final double branchRate = branchRateModel.getBranchRate(tree, childNode);
        final double substTime = (tree.getNodeHeight(parentNode) - tree.getNodeHeight(childNode));

        for (int r = 0; r < markovjumps.size(); r++) {
            MarkovJumpsSubstitutionModel thisMarkovJumps = markovjumps.get(r);
            if (useUniformization) {
                computeSampledMarkovJumpsForBranch(((UniformizedSubstitutionModel) thisMarkovJumps), substTime,
                        branchRate, childNum, parentStates, childStates, probabilities, scaleByTime[r],
                        expectedJumps.get(r), rateCategory);
            } else {
                computeIntegratedMarkovJumpsForBranch(thisMarkovJumps, substTime, branchRate, childNum, parentStates,
                        childStates, probabilities, condJumps, scaleByTime[r], expectedJumps.get(r), rateCategory);
            }
        }
    }

    private void computeSampledMarkovJumpsForBranch(UniformizedSubstitutionModel thisMarkovJumps,
                                                    double substTime,
                                                    double branchRate,
                                                    int childNum,
                                                    int[] parentStates,
                                                    int[] childStates,
                                                    double[] probabilities,
                                                    boolean scaleByTime,
                                                    double[][] thisExpectedJumps,
                                                    int[] rateCategory) {


        // Fill condJumps with sampled values for this branch for each site
        for (int j = 0; j < patternCount; j++) {
            final int category = rateCategory == null ? 0 : rateCategory[j];
            final double categoryRate = siteRateModel.getRateForCategory(category);
            final int matrixIndex = category * stateCount * stateCount;
            double value = thisMarkovJumps.computeCondStatMarkovJumps(
                    parentStates[j],
                    childStates[j],
                    substTime * branchRate * categoryRate,
                    probabilities[matrixIndex + parentStates[j] * stateCount + childStates[j]]
            );
            if (scaleByTime) {
                value /= branchRate * categoryRate;
            }
            thisExpectedJumps[childNum][j] = value;
        }
    }

    private void computeIntegratedMarkovJumpsForBranch(MarkovJumpsSubstitutionModel thisMarkovJumps,
                                                       double substTime,
                                                       double branchRate,
                                                       int childNum,
                                                       int[] parentStates,
                                                       int[] childStates,
                                                       double[] probabilities,
                                                       double[][] condJumps,
                                                       boolean scaleByTime,
                                                       double[][] thisExpectedJumps,
                                                       int[] rateCategory) {

        // Fill condJumps with conditional mean values for this branch
        for (int i = 0; i < categoryCount; i++) {

            double rate = siteRateModel.getRateForCategory(i);
            if (rate > 0) {

                if (categoryCount == 1) {
                    thisMarkovJumps.computeCondStatMarkovJumps(
                            substTime * branchRate * rate,
                            probabilities,
                            condJumps[i]);
                } else {

                    System.arraycopy(probabilities, i * stateCount * stateCount, tmpProbabilities, 0,
                            stateCount * stateCount);

                    thisMarkovJumps.computeCondStatMarkovJumps(
                            substTime * branchRate * rate,
                            tmpProbabilities,
                            condJumps[i]);
                }

                if (scaleByTime) {
                    double scalar = branchRate * rate;
                    for (int j = 0; j < condJumps[i].length; j++) {
                        condJumps[i][j] /= scalar;
                    }
                }
            } else {
                Arrays.fill(condJumps[i], 0.0);
                if (thisMarkovJumps.getType() == MarkovJumpsType.REWARDS && scaleByTime) {
                    for (int j = 0; j < stateCount; j++) {
                        condJumps[i][j * stateCount + j] = substTime;
                    }
                }
            }
        }

        for (int j = 0; j < patternCount; j++) { // Pick out values given parent and child states
            int category = rateCategory == null ? 0 : rateCategory[j];
            thisExpectedJumps[childNum][j] = condJumps[category][parentStates[j] * stateCount + childStates[j]];
        }
    }

    public LogColumn[] getColumns() {

        int nColumns = patternCount * numRegisters;
        if (reportUnconditionedColumns) {
            if (categoryCount == 1) {
                nColumns += numRegisters;
            } else {
                nColumns *= 2;
            }
        }

        int index = 0;
        LogColumn[] allColumns = new LogColumn[nColumns];
        for (int r = 0; r < numRegisters; r++) {
            for (int j = 0; j < patternCount; j++) {
                allColumns[index++] = new ConditionedCountColumn(jumpTag.get(r), r, j);
                if (reportUnconditionedColumns) {
                    if (categoryCount > 1) {
                        allColumns[index++] = new UnconditionedCountColumn(jumpTag.get(r), r, j, rateCategory);
                    }
                }
            }
            if (reportUnconditionedColumns) {
                if (categoryCount == 1) {
                    allColumns[index++] = new UnconditionedCountColumn(jumpTag.get(r), r);
                }
            }
        }
        return allColumns;
    }

    protected abstract class CountColumn extends NumberColumn {
        protected int indexRegistration;
        protected int indexSite;

        public CountColumn(String label, int r, int j) {
            super(label + (j >= 0 ? "[" + (j + 1) + "]" : ""));
            indexRegistration = r;
            indexSite = j;
        }

        public abstract double getDoubleValue();
    }

    protected class ConditionedCountColumn extends CountColumn {

        public ConditionedCountColumn(String label, int r, int j) {
            super("c_" + label, r, j);
        }

        public double getDoubleValue() {
            double total = 0;
            double[][] values = getMarkovJumpsForRegister(treeModel, indexRegistration);
            for (int i = 0; i < treeModel.getNodeCount(); i++) {
                total += values[i][indexSite];
            }
            return total;
        }
    }

    protected class UnconditionedCountColumn extends CountColumn {
        int[] rateCategory;

        public UnconditionedCountColumn(String label, int r, int j, int[] rateCategory) {
            super("u_" + label, r, j);
            this.rateCategory = rateCategory;
        }

        public UnconditionedCountColumn(String label, int r) {
            this(label, r, -1, null);
        }

        public double getDoubleValue() {
            double value = markovjumps.get(indexRegistration).getMarginalRate() * getExpectedTreeLength();
            if (rateCategory != null) {
                value *= siteRateModel.getRateForCategory(rateCategory[indexSite]);
            }
            return value;
        }

        private double getExpectedTreeLength() {
            double expectedTreeLength = 0;
            for (int i = 0; i < treeModel.getNodeCount(); i++) {
                NodeRef node = treeModel.getNode(i);
                if (!treeModel.isRoot(node)) {
                    expectedTreeLength += branchRateModel.getBranchRate(treeModel, node)
                            * treeModel.getBranchLength(node);
                }
            }
            return expectedTreeLength;
        }
    }

    private List<MarkovJumpsSubstitutionModel> markovjumps;
    private List<Parameter> registerParameter;
    private List<String> jumpTag;
    private List<double[][]> expectedJumps;
    private boolean[] scaleByTime;
    private double[] tmpProbabilities;
    private double[][] condJumps;
    private int numRegisters;
    private final boolean useUniformization;
    private final int nSimulants;
    private final boolean reportUnconditionedColumns;
}
