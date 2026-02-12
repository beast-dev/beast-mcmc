package dr.evomodel.treedatalikelihood.markovjumps;

import dr.evolution.tree.NodeRef;
import dr.evolution.tree.Tree;
import dr.evolution.tree.TreeTrait;
import dr.evolution.tree.TreeTraitProvider;
import dr.evomodel.siteratemodel.SiteRateModel;
import dr.evomodel.substmodel.MarkovJumpsSubstitutionModel;
import dr.evomodel.substmodel.SubstitutionModel;
import dr.evomodel.substmodel.UniformizedSubstitutionModel;
import dr.evomodel.treedatalikelihood.preorder.AbstractRealizedDiscreteTraitDelegate;
import dr.evomodel.treelikelihood.MarkovJumpsTraitProvider;
import dr.inference.markovjumps.MarkovJumpsType;
import dr.inference.markovjumps.StateHistory;

import java.util.ArrayList;
import java.util.List;

public class CompleteHistoryAddOn implements AbstractRealizedDiscreteTraitDelegate.RealizedDiscreteAddOn {

    public CompleteHistoryAddOn(String name,
                                Tree tree,
                                List<SubstitutionModel> branchSubstitutionModels,
                                SiteRateModel siteRateModel,
                                AbstractRealizedDiscreteTraitDelegate ancestralTraitDelegate,
                                MarkovJumpsTraitProvider.ValueScaling valueScaling,
                                boolean logCompleteHistory) {

        this.name = name;
        this.siteRateModel = siteRateModel;
        this.ancestralTraitDelegate = ancestralTraitDelegate;
        this.stateCount = branchSubstitutionModels.get(0).getDataType().getStateCount();

        this.markovJumps = new ArrayList<>();
        for (SubstitutionModel substitutionModel : branchSubstitutionModels) {
            UniformizedSubstitutionModel usm = new UniformizedSubstitutionModel(substitutionModel,
                    MarkovJumpsType.HISTORY, 1);
            usm.setSaveCompleteHistory(true);

            double[] registration = new double[stateCount * stateCount];
            for (int i = 0; i < stateCount; ++i) {
                for (int j = 0; j < stateCount; ++j) {
                    if (i != j) {
                        registration[i * stateCount + j] = 1.0;
                    }
                }
            }
            usm.setRegistration(registration);

            this.markovJumps.add(usm);
        }

        this.valueScaling = valueScaling;

        this.expectedJumps = new double[2 * tree.getNodeCount() - 2][];
        this.histories = new String[2 * tree.getNodeCount() - 2][];
        this.stateHistories = new StateHistory[2 * tree.getNodeCount() - 2][];
    }

    private final SiteRateModel siteRateModel;
    private final AbstractRealizedDiscreteTraitDelegate ancestralTraitDelegate;

    private final List<MarkovJumpsSubstitutionModel> markovJumps;
    private final double[][] expectedJumps;
    private final String[][] histories;
    private final StateHistory[][] stateHistories;

    private final String name;
    private final int stateCount;

    private boolean useCompactHistory = true; // TODO final
    private final MarkovJumpsTraitProvider.ValueScaling valueScaling;

    private final String stateHistoryName = "stateHistory";

    @Override
    public void constructTraits(TreeTraitProvider.Helper treeTraitHelper) {

        TreeTrait<StateHistory[]> stateHistoryTrait = new TreeTrait<StateHistory[]>() {

            @Override
            public String getTraitName() {
                return name + ".raw";
            }

            @Override
            public Intent getIntent() {
                return Intent.BRANCH;
            }

            @Override
            public Class getTraitClass() {
                return StateHistory[].class;
            }

            @Override
            public StateHistory[] getTrait(Tree tree, NodeRef node) {
                ancestralTraitDelegate.getStatesForNode(tree, node);
                return stateHistories[node.getNumber()];
            }

            @Override
            public String getTraitString(Tree tree, NodeRef node) {
                throw new RuntimeException("Not loggable");
            }

            @Override
            public boolean getLoggable() {
                return false;
            }
        };

        treeTraitHelper.addTrait(stateHistoryTrait);

        if (histories != null) {
            TreeTrait<String[]> stateHistoryStringTrait = new TreeTrait.SA() {

                public String getTraitName() {
                    return name + ".history";
                }

                public Intent getIntent() {
                    return Intent.BRANCH;
                }

                public Class getTraitClass() {
                    return String[].class;
                }

                public String[] getTrait(Tree tree, NodeRef node) {
                    ancestralTraitDelegate.getStatesForNode(tree, node);
                    return histories[node.getNumber()];
                }
            };

            treeTraitHelper.addTrait(stateHistoryStringTrait);
        }

        TreeTrait<double[]> totalCountTrait = new TreeTrait.DA() {

            public String getTraitName() {
                return name + ".count";
            }

            public Intent getIntent() {
                return Intent.BRANCH;
            }

            public Class getTraitClass() {
                return double[].class;
            }

            public double[] getTrait(Tree tree, NodeRef node) {
                ancestralTraitDelegate.getStatesForNode(tree, node);
                return expectedJumps[node.getNumber()];
            }
        };

        treeTraitHelper.addTrait(totalCountTrait);
    }

    @Override
    public void hookCalculation(int node, int parent, double[] probabilities) {
        AbstractRealizedDiscreteTraitDelegate.OrderedEvents events = ancestralTraitDelegate.getOrderedInformation(node, parent, probabilities);
        computeSampledMarkovJumpsForBranch((UniformizedSubstitutionModel) markovJumps.get(0), events,
                node, expectedJumps, stateHistories, histories);
    }

    private void computeSampledMarkovJumpsForBranch(UniformizedSubstitutionModel markovJumpsModel,
                                                    AbstractRealizedDiscreteTraitDelegate.OrderedEvents events,
                                                    int destinationIndex,
                                                    double[][] expectedJumps,
                                                    StateHistory[][] stateHistories,
                                                    String[][] histories) {

        final int[] startingStates = events.getStartingStates();
        final int[] endingStates = events.getEndingStates();
        final int[] categories = events.getCategories();
        final double[] probabilities = events.getTransitionMatrix();
        final int patternCount = startingStates.length;

        for (int j = 0; j < patternCount; j++) {

            final int category = categories[j];
            final double categoryRate = siteRateModel.getRateForCategory(category);
            final int matrixIndex = category * stateCount * stateCount;

            double value = markovJumpsModel.computeCondStatMarkovJumps(
                    startingStates[j],
                    endingStates[j],
                    events.getBranchTime() * events.getBranchRate() * categoryRate,
                    probabilities[matrixIndex + startingStates[j] * stateCount + endingStates[j]]);

            if (valueScaling == MarkovJumpsTraitProvider.ValueScaling.BY_TIME) {
                value /= events.getBranchRate() * categoryRate;
            }

            if (expectedJumps[destinationIndex] == null) {
                expectedJumps[destinationIndex] = new double[patternCount];
            }
            expectedJumps[destinationIndex][j] = value;

            StateHistory sh = markovJumpsModel.getStateHistory();
            sh.rescaleTimesOfEvents(events.getStartingTime(), events.getEndingTime());

            if (stateHistories != null) {
                if (stateHistories[destinationIndex] == null) {
                    stateHistories[destinationIndex] = new StateHistory[patternCount];
                }
                stateHistories[destinationIndex][j] = sh;
            }

            if (histories != null) {
                if (histories[destinationIndex] == null) {
                    histories[destinationIndex] = new String[patternCount];
                }

                int site = (useCompactHistory) ? j + 1 : -1;
                String str1 =  sh.toStringChanges(site, markovJumpsModel.getSubstitutionModel().getDataType());
                histories[destinationIndex][j] = str1;
            }

            int m = 0; ++m;

        }
    }

//    @Override
//    public TreeTrait[] getTreeTraits() {
//        return new TreeTrait[0];
//    }
//
//    @Override
//    public TreeTrait getTreeTrait(String key) {
//        return null;
//    }

//    private void computeSampledMarkovJumpsForBranch(UniformizedSubstitutionModel thisMarkovJumps,
//                                                    double substTime,
//                                                    double branchRate,
//                                                    int childNum,
//                                                    int[] parentStates,
//                                                    int[] childStates,
//                                                    double unused1,
//                                                    double childTime,
//                                                    double[] probabilities,
//                                                    boolean scaleByTime,
//                                                    double[][] thisExpectedJumps,
//                                                    int[] rateCategory,
//                                                    boolean saveHistory) {
//
//
//        // Fill condJumps with sampled values for this branch for each site
//        for (int j = 0; j < patternCount; j++) {
//            final int category = rateCategory == null ? 0 : rateCategory[j];
//            final double categoryRate = siteRateModel.getRateForCategory(category);
//            final int matrixIndex = category * stateCount * stateCount;
//            double value = thisMarkovJumps.computeCondStatMarkovJumps(
//                    parentStates[j],
//                    childStates[j],
//                    substTime * branchRate * categoryRate,
//                    probabilities[matrixIndex + parentStates[j] * stateCount + childStates[j]]
//            );
//            if (scaleByTime) {
//                value /= branchRate * categoryRate;
//            }
//            thisExpectedJumps[childNum][j] = value;
//            if (saveHistory) {
//                int site = (useCompactHistory) ? j + 1 : -1;
//                histories[childNum][j] = thisMarkovJumps.getCompleteHistory(site, parentTime, childTime);
//            }
//        }
//    }

}
