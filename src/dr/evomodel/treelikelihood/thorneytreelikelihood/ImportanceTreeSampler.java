package dr.evomodel.treelikelihood.thorneytreelikelihood;

import dr.evolution.coalescent.Intervals;
import dr.evomodel.coalescent.GMRFMultilocusSkyrideLikelihood;
import dr.evomodel.coalescent.TreeIntervals;
import dr.evomodel.operators.AbstractTreeOperator;
import dr.evomodel.tree.TreeChangedEvent;
import dr.inference.operators.GibbsOperator;
import dr.math.MathUtils;

import java.util.*;
import java.util.concurrent.*;

//An important importance sampler
public class ImportanceTreeSampler extends AbstractTreeOperator  implements GibbsOperator {
    private static final String IMPORTANCE_TREE_SAMPLER = "importanceTreeSampler";

    private final ArrayList<SDismProposal> proposals;
    private ConstrainedTreeModel tree;
    //TODO needs to be treeIntervals I recon
    private Intervals intervalList;
    private ExecutorService threadPool;
    List<Callable<Double>> proposalCallers = new ArrayList<Callable<Double>>();

    public ImportanceTreeSampler(double weight, int samples, int threads, ThorneyTreeLikelihood treeLikelihood,GMRFMultilocusSkyrideLikelihood skygrid) {
        setWeight(weight);
        this.tree = (ConstrainedTreeModel) treeLikelihood.getTreeModel();
        this.intervalList = (Intervals) ((TreeIntervals)skygrid.getIntervalList()).getIntervals(); // this is a hard requirement
        this.proposals = new ArrayList<>();
        this.threadPool = threads > 1 ? Executors.newFixedThreadPool(threads) : Executors.newCachedThreadPool();

        if (samples < threads) {
            throw new IllegalArgumentException("number of samples should be less than threads requested");
        }
        int samplesEach = samples / threads;
        int samplesExtras = samples % threads;
        for (int i = 0; i < threads; i++) {
            int sampleCount;
            if (i < samplesExtras) {
                sampleCount = samplesEach + 1;
            } else {
                sampleCount = samplesEach;
            }
            SDismProposal proposal = new SDismProposal(sampleCount, treeLikelihood, skygrid);
            proposals.add(proposal);
            proposalCallers.add(new ProposalCaller(proposal, i));
        }
    }

class ProposalCaller implements Callable<Double> {

    public ProposalCaller(SDismProposal proposal, int index) {
        this.proposal = proposal;
        this.index = index;
    }

    public Double call() throws Exception {
        return proposal.sampleTree();
    }

    private final SDismProposal proposal;
    private final int index;
}

    @Override
    public String getOperatorName() {
        return IMPORTANCE_TREE_SAMPLER;
    }

    /**
     * Called by operate(), does the actual operation.
     *
     * @return the hastings ratio
     */
    @Override
    public double doOperation()  {

        double[] logWeights = new double[proposals.size()];
        try {
            List<Future<Double>> results = threadPool.invokeAll(proposalCallers);

            int i = 0;

            for (Future<Double> result : results) {
                double logW = result.get();
                logWeights[i] = logW;
                i += 1;
            }
        }catch (InterruptedException | ExecutionException e) {
            e.printStackTrace();
        }

        int choice = MathUtils.randomChoiceLogPDF(logWeights);

        ConstrainedTreeModel chosenTree = proposals.get(choice).getTree();

        tree.copyEdgesAndHeights(chosenTree);
        tree.fireModelChanged(TreeChangedEvent.create());

        Intervals chosenIntervals = (Intervals) proposals.get(choice).getIntervals();
        intervalList.copyIntervals(chosenIntervals);
        return 0.0;
    }


}
