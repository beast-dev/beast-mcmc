package dr.evomodel.branchmodel;

import dr.evolution.tree.NodeRef;
import dr.evolution.tree.Tree;
import dr.evomodel.tree.TreeStatistic;
import dr.evomodelxml.branchmodel.EpochBranchBreakingStatisticParser;
import dr.stats.DiscreteStatistics;

import java.util.ArrayList;
import java.util.List;

public class EpochBranchBreakingStatistic extends TreeStatistic {
    private final String mode;
    private final EpochBranchModel epochModel;

//    private List<double[]> proportionList = new ArrayList<>();


    public EpochBranchBreakingStatistic(String name, EpochBranchModel epochModel, String mode){
        super(name);

        this.mode = mode;
        this.epochModel = epochModel;
    }

    @Override
    public void setTree(Tree tree) {

    }

    @Override
    public Tree getTree() {
        return null;
    }

    @Override
    public int getDimension() {
        return 1;
    }

//    private void fillAllBranches() {
//        proportionList.clear();
//
//        Tree tree = epochModel.getTree();
//        int idx = 0;
//        for (int i = 0; i < tree.getNodeCount() - 1; i++) {
//            NodeRef node = tree.getNode(i);
//            if ( !tree.isRoot(node) ) {
//                doBranch(node);
//                idx++;
//            }
//        }
//    }
//
//    private int findEpochInterval(double t, double[] breaks) {
//        for (int i = 0; i < breaks.length - 1; i++) {
//            if (t >= breaks[i] && t < breaks[i+1]) {
//                return i;
//            }
//        }
//        return -1;
//    }
//
//    private void doBranch(NodeRef node) {
//        double[] breaks = epochModel.getEpochTimesDouble();
//        double[] paddedBreaks = new double[breaks.length + 2];
//        System.arraycopy(breaks,0,paddedBreaks,1,breaks.length);
//        paddedBreaks[paddedBreaks.length - 1] = Double.POSITIVE_INFINITY;
//
//        double tYoung = epochModel.getTree().getNodeHeight(node);
//        double tOld = epochModel.getTree().getNodeHeight(epochModel.getTree().getParent(node));
//        double tTot = tOld - tYoung;
//
//        int firstEpochIndex = findEpochInterval(tYoung,paddedBreaks);
//        int lastEpochIndex = findEpochInterval(tOld,paddedBreaks);
//
//        if (firstEpochIndex == lastEpochIndex) {
//            proportionList.add(new double[]{1.0});
//        } else {
//            int nbreaks = lastEpochIndex - firstEpochIndex;
//            double[] props = new double[nbreaks + 1];
//            double tStart = tYoung;
//            for (int i = firstEpochIndex; i < lastEpochIndex; i++) {
//                double tEnd = paddedBreaks[i + 1];
//                props[i - firstEpochIndex] =  (tEnd - tStart)/tTot;
//                tStart = tEnd;
//            }
//            props[nbreaks] = (tOld - paddedBreaks[lastEpochIndex])/tTot;
//
//            double s = 0.0;
//            for (int i = 0; i < props.length; i++) {
//                s += props[i];
//            }
//            if ( Math.abs(s - 1.0) > 1e-6 ) {
//                throw new RuntimeException("Uh oh");
//            }
//
//            proportionList.add(props);
//        }
//    }
//
//    @Override
//    public double getStatisticValue(int dim) {
//        fillAllBranches();
//        if (mode.equals(EpochBranchBreakingStatisticParser.PROP_BROKEN)) {
//            int nbroken = 0;
//            for (int i = 0; i < proportionList.size(); i++) {
//                if (proportionList.get(i).length > 1) {
//                    nbroken++;
//                }
//            }
//            return ((double)nbroken)/((double)proportionList.size());
//        } else if (mode.equals(EpochBranchBreakingStatisticParser.MEAN_MAX)) {
//            double sum = 0.0;
//            for (int i = 0; i < proportionList.size(); i++) {
//                sum += DiscreteStatistics.max(proportionList.get(i));
//            }
//            return sum/proportionList.size();
//        }
//        return Double.NaN;
//    }

    @Override
    public double getStatisticValue(int dim) {
        Tree tree = epochModel.getTree();

        double stat = 0.0;
        if (mode.equals(EpochBranchBreakingStatisticParser.PROP_BROKEN)) {
            for (int i = 0; i < tree.getNodeCount() - 1; i++) {
                NodeRef node = tree.getNode(i);
                if (!tree.isRoot(node)) {
                    BranchModel.Mapping map = epochModel.getBranchModelMapping(node);
                    if (map.getWeights().length > 1) {
                        stat += 1.0;
                    }
                }
            }
            stat /= (tree.getNodeCount() - 1);
        } else if (mode.equals(EpochBranchBreakingStatisticParser.MEAN_MAX)) {
            for (int i = 0; i < tree.getNodeCount() - 1; i++) {
                NodeRef node = tree.getNode(i);
                if (!tree.isRoot(node)) {
                    BranchModel.Mapping map = epochModel.getBranchModelMapping(node);
                    if (map.getWeights().length == 1) {
                        stat += 1.0;
                    } else {
                        double[] weights = map.getWeights();
                        double sum = 0.0;
                        for (double w : weights) {
                            sum += w;
                        }
                        stat += DiscreteStatistics.max(weights)/sum;
                    }
                }
            }
            stat /= (tree.getNodeCount() - 1);
        }
        return stat;
    }

}
