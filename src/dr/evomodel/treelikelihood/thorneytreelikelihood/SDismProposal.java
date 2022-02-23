package dr.evomodel.treelikelihood.thorneytreelikelihood;

import dr.evolution.coalescent.IntervalList;
import dr.evolution.coalescent.Intervals;
import dr.evolution.tree.NodeRef;
import dr.evomodel.branchratemodel.StrictClockBranchRates;
import dr.evomodel.coalescent.GMRFMultilocusSkyrideLikelihood;
import dr.math.Binomial;
import dr.math.LogTricks;
import dr.math.MathUtils;
import dr.util.HeapSort;
import org.apache.commons.math.random.MersenneTwister;
import org.apache.commons.math.random.RandomGenerator;

import java.io.Serializable;
import java.util.*;

public class SDismProposal implements Serializable {
    final ConstrainedTreeModel tree;
    final Intervals intervals;
    final Intervals activeIntervals;
    final BranchLengthProvider branchlengthProvider;
    final PoissonBranchLengthLikelihoodDelegate branchLengthLikelihood;
    final GMRFMultilocusSkyrideLikelihood skygridLikelihood;
    final NodeRef[] sampledTips;
    Map<ConstrainedTreeModel.WrappingSubtree, Clade> clades;
    int activeNodeCount;
    int currentTip;
    int coalescingNodeStat;
    int mutatingNodeCount;
    final int totalCoalescentEvents;
    int simulatedCoalescentEvents;
    int currentGridPoint;
    final double[] gridPoints;
    static final double tol = 0;
    double scaledMutationRate;
    double currentLogWeight;
    double theta;
    final List<NodeRef> availableNodes;
    final RandomGenerator random;
    private int samples;
    public SDismProposal(int samples,ThorneyTreeLikelihood treeLikelihood, GMRFMultilocusSkyrideLikelihood skygridLikelihood) {
        this.tree = new ConstrainedTreeModel( (ConstrainedTreeModel) treeLikelihood.getTreeModel()); // This tree is a copy of the base tree; Note this will not sample the tips
        this.branchlengthProvider = treeLikelihood.getBranchLengthProvider();

        this.branchLengthLikelihood = (PoissonBranchLengthLikelihoodDelegate) treeLikelihood.getThorneyBranchLengthLikelihoodDelegate();

        if (!(this.branchLengthLikelihood.getBranchRateModel() instanceof StrictClockBranchRates)) {
            throw new RuntimeException("Importance tree sampler is only implemented for a stick clock branch rate model");
        }
        this.skygridLikelihood = skygridLikelihood;
        this.samples = samples;

        intervals = new Intervals(tree.getNodeCount(),false);
        activeIntervals = new Intervals(tree.getNodeCount(),false);

        availableNodes = new ArrayList<>();
        activeNodeCount = 0;
        currentTip = 0;
        coalescingNodeStat = 0;
        mutatingNodeCount = 0;

        sampledTips = new NodeRef[tree.getExternalNodeCount()];
        totalCoalescentEvents = tree.getInternalNodeCount();
        simulatedCoalescentEvents = 0;
        currentGridPoint = 0;
        gridPoints = skygridLikelihood.getGridPoints();

        scaledMutationRate = (double) branchLengthLikelihood.getBranchRateModel().getVariable(0).getValue(0) * branchLengthLikelihood.getScale();
        theta = 0;

        // TODO redo if tip heights change
        double[] heights = new double[tree.getExternalNodeCount()];
        for (int i = 0; i < tree.getExternalNodeCount(); i++) {
            heights[i] = tree.getNodeHeight(tree.getExternalNode(i));
        }
        int[] indices = new int[tree.getExternalNodeCount()];
        HeapSort.sort(heights, indices);

        for (int i = 0; i < sampledTips.length; i++) {
            sampledTips[i] = tree.getExternalNode(indices[i]);
        }
        this.random = new MersenneTwister();
        this.random.setSeed(MathUtils.nextInt());
    }



    public double sampleTree() {
        double totalLogWeight = 0;
        for (int i = 0; i < samples; i++) {
            sampleNextTree();
            totalLogWeight = i==0?currentLogWeight: LogTricks.logSum(totalLogWeight, currentLogWeight);
            if (Math.log(random.nextDouble()) < currentLogWeight - totalLogWeight) {
                tree.storeState();
                intervals.copyIntervals(activeIntervals);
            }
        }
        tree.restoreState();
        return totalLogWeight;
    }

    public ConstrainedTreeModel getTree() {
        return tree;
    }
    public IntervalList getIntervals() {
        return intervals;
    }

    private void sampleNextTree() {
        tree.beginTreeEdit();
        setup();

        double currentHeight = 0.0;
        setCurrentHeight(currentHeight);

        while (coalescingNodeStat < 1 && mutatingNodeCount < 1) {
            currentHeight = getMinimumInactiveHeight();
            setCurrentHeight(currentHeight);
        }


        while (!done()) {
            Event event = getNextEvent();
            if (event.time + currentHeight > getNextFixedEventHeight()) { // or if we hit gridpoint
                currentHeight = getNextFixedEventHeight();
                setCurrentHeight(currentHeight);
            } else {
                currentHeight += event.time;
                applyEvent(event.node, currentHeight);
            }
            if (!done()) {
                while (coalescingNodeStat < 1 && mutatingNodeCount < 1) {
                    currentHeight = getNextFixedEventHeight();
                    setCurrentHeight(currentHeight);
                }
            }
        }

        for (Clade clade :
                clades.values()) {
            assert clade.complete;
        }
        assert availableNodes.isEmpty();
        assert activeNodeCount == 0;
        assert currentTip == tree.getExternalNodeCount();
        assert coalescingNodeStat == 0;
        assert mutatingNodeCount == 0;

        tree.endTreeEditQuietly();

    }

    void applyEvent(NodeRef node, double currentHeight) {
        ConstrainedTreeModel.WrappingSubtree subtree = (ConstrainedTreeModel.WrappingSubtree) tree.getSubtree(node, ConstrainedTreeModel.SubtreeContext.IncludeTips);
        NodeRef subtreeNode = subtree.getWrappingNode(node, ConstrainedTreeModel.SubtreeContext.IncludeTips);
        Clade clade = clades.get(subtree);
        StateChange sc = clade.actOnNode(subtreeNode, currentHeight);
        updateState(sc);
    }

    boolean done() {
        return simulatedCoalescentEvents == totalCoalescentEvents;
    }

    void updateState(StateChange change) {
        this.mutatingNodeCount += change.mutDelta;
        this.coalescingNodeStat += change.coalDelta;
        this.currentLogWeight += change.logEventWeight;
    }

    void activateNode(NodeRef nodeRef) {
        assert !tree.isRoot(nodeRef);
        ConstrainedTreeModel.WrappingSubtree subtree = (ConstrainedTreeModel.WrappingSubtree) tree.getSubtree(nodeRef, ConstrainedTreeModel.SubtreeContext.IncludeTips);
        NodeRef subtreeNode = subtree.getWrappingNode(nodeRef, ConstrainedTreeModel.SubtreeContext.IncludeTips);
        assert subtree.isExternal(subtreeNode);
        Clade clade = clades.get(subtree);

        activeNodeCount += 1;
        StateChange change = clade.activateSubtreeNode(subtreeNode);
        updateState(change);

        //update total counts
    }

    /**
     * @return the height of youngest inactive node.
     */
    double getMinimumInactiveHeight() {
        if (currentTip < sampledTips.length) {

            return tree.getNodeHeight(sampledTips[currentTip]);
        } else return Double.POSITIVE_INFINITY;
    }

    double getNextGridPoint() {
        if (currentGridPoint < gridPoints.length) {
            return gridPoints[currentGridPoint];
        } else return Double.POSITIVE_INFINITY;

    }

    double getNextFixedEventHeight() {
        return Math.min(getMinimumInactiveHeight(), getNextGridPoint());
    }

    /**
     * Set the current height.
     */
    void setCurrentHeight(double height) {
        assert height != Double.POSITIVE_INFINITY;
        while (getMinimumInactiveHeight() <= height + tol) {
            NodeRef nodeRef = sampledTips[currentTip];
            activateNode(nodeRef);
            activeIntervals.addSampleEvent(tree.getNodeHeight(nodeRef),nodeRef.getNumber());
            currentTip += 1;
        }
        while (getNextGridPoint() <= height + tol) {
            currentGridPoint += 1;
        }
    }

    /***
     * A private function to clear the topology of the tree.
     * The subtree topology is maintained
     */
    void clearTreeTopology() {
        tree.clearTopology();
    }

    void setup() {


        currentLogWeight = 0;

        clearTreeTopology();
        this.clades = new HashMap<ConstrainedTreeModel.WrappingSubtree, Clade>();
        for (int i = 0; i < tree.getSubtreeCount(); i++) {
            ConstrainedTreeModel.WrappingSubtree subtree = (ConstrainedTreeModel.WrappingSubtree) tree.getSubtree(i);
            Clade clade = new Clade(subtree);
            clades.put(subtree, clade);
        }

        scaledMutationRate = (double) branchLengthLikelihood.getBranchRateModel().getVariable(0).getValue(0) * branchLengthLikelihood.getScale();

        activeNodeCount = 0;
        currentTip = 0;
        coalescingNodeStat = 0;
        mutatingNodeCount = 0;
        simulatedCoalescentEvents = 0;
        currentGridPoint = 0;
        activeIntervals.resetEvents();
    }

    Event getNextEvent() {
        //Time to the next event provided it's an acceptable event
        // mutation in the right place or coalescent in the right place
        // I beleive this matches Donnelly 2000 proprosal distribution which i believe is uniform across lineages.


        double currentNe = Math.exp(skygridLikelihood.getPopSizeParameter().getParameterValue(currentGridPoint));
        double totalMutationRate = scaledMutationRate * (double) activeNodeCount;

        //TODO why no factor of 2?
        theta = currentNe * scaledMutationRate;//? why not *2 here


        double totalCoalescentRate = Binomial.choose2(activeNodeCount) / currentNe;//check

        double totalRate = totalMutationRate + totalCoalescentRate;

        assert mutatingNodeCount + coalescingNodeStat > 0;

        //get time to event;
        double time = nextExponential(totalRate);
        //get type of event;
        EventType type = null;
        // Pick a node at random. The event is then determined by the node.
        NodeRef eventNode = availableNodes.get(random.nextInt(availableNodes.size()));

        return new Event(time, eventNode);
    }
    private double nextExponential(double lambda){
        return -1.0 * Math.log(1 - random.nextDouble()) / lambda;
    }

    private static class StateChange {
        public final int mutDelta;
        public final int coalDelta;
        public final double logEventWeight;

        public StateChange(int mutDelta, int coalDelta, double logEventWeight) {
            this.mutDelta = mutDelta;
            this.coalDelta = coalDelta;
            this.logEventWeight = logEventWeight;
        }
    }

    private enum EventType {
        MUTATION, COALESCENT;
    }

    private static class Event {
        public final double time;
        public final NodeRef node;

        public Event(double time, NodeRef node) {
            this.time = time;
            this.node = node;
        }

    }

    private class Clade {
        private final ConstrainedTreeModel.WrappingSubtree subtree;
        private final List<NodeRef> coalescingNodes;
        private final Map<NodeRef, Integer> mutatingNodes;
        private final Stack<NodeRef> nodeStack;

        private boolean complete;

        private Clade(ConstrainedTreeModel.WrappingSubtree subtree) {
            this.subtree = subtree;
            this.complete = false;

            this.coalescingNodes = new ArrayList<NodeRef>();
            this.mutatingNodes = new HashMap<NodeRef, Integer>();

            this.nodeStack = new Stack<NodeRef>();

            this.nodeStack.add(subtree.getRoot());
            for (int i = 0; i < subtree.getInternalNodeCount(); i++) {
                NodeRef node = subtree.getInternalNode(i);
                if (!subtree.isRoot(node)) {
                    this.nodeStack.add(node);
                }
            }

        }

        public boolean isComplete() {
            return complete;
        }


        public int mutatingNodeCount() {
            return mutatingNodes.size();
        }

        public int coalescingPairsCount() {
            return coalescingNodeStat();

            // This should be the number of coalescing nodes
            // how many pairs are active
//            int liveLin = coalescingNodes.size();
//            return (liveLin - 1) * liveLin / 2;
        }

        public int activePairsCount() {
            int liveLin = coalescingNodes.size() + mutatingNodeCount();
            return (liveLin - 1) * liveLin / 2;
        }

        /**
         * A helper function to activate a node make it available or mutating as
         * indicated by data.
         *
         * @param node
         */
        public StateChange activateSubtreeNode(NodeRef node) {
            assert !coalescingNodes.contains(node);


            // check if needs mutations
            // only needs mutations if this is a tip
            if (branchlengthProvider.getBranchLength(tree, subtree.getNodeInWrappedTree(node)) > 0.5) {
                this.mutatingNodes.put(node, 0);
                availableNodes.add(subtree.getNodeInWrappedTree(node)); //available for events
                return new StateChange(1, 0, 0);
            } else {
                int oldCoalescingCount = coalescingPairsCount();
                if (coalescingNodes.size() > 0) {
                    availableNodes.add(subtree.getNodeInWrappedTree(node)); //available for events
                    if (coalescingNodes.size() == 1) {
                        // this node is now available to coalesce
                        availableNodes.add(subtree.getNodeInWrappedTree(coalescingNodes.get(0)));
                    }
                }

                this.coalescingNodes.add(node);
                int newCoalescingCount = coalescingPairsCount();
                return new StateChange(0, newCoalescingCount - oldCoalescingCount, 0);
            }
        }

        public StateChange actOnNode(NodeRef nodeRef, double time) {
            if (mutatingNodes.containsKey(nodeRef)) {
                return mutateNode(nodeRef);
            } else {
                return coalesceNodes(nodeRef, time);
            }
        }

        /**
         * Add a mutation to the node's lineage and if we are done add this node to
         * those available for coalescence
         */
        public StateChange mutateNode(NodeRef node) {

            int current_muts = mutatingNodes.get(node) + 1;
            int expected_muts = (int) branchlengthProvider.getBranchLength(tree, subtree.getNodeInWrappedTree(node));
            assert current_muts <= expected_muts;
            if (current_muts < expected_muts) {
                mutatingNodes.put(node, current_muts);
                // case 2
                double eventWeight = ((double) coalescingNodeStat + (double)  mutatingNodeCount) / (double) activeNodeCount * (theta / (activeNodeCount - 1.0 + theta));

                return new StateChange(0, 0, Math.log(eventWeight));
            } else {
                double eventWeight;
                if (coalescingNodes.size() > 0) { // there are others of this type now
                    eventWeight = ((double) coalescingNodeStat + (double) mutatingNodeCount) / (double) activeNodeCount * ((coalescingNodes.size() + 1.0) * theta / (activeNodeCount - 1.0 + theta));
                    if (coalescingNodes.size() == 1) {
                        //this node is now available
                        availableNodes.add(subtree.getNodeInWrappedTree(coalescingNodes.get(0)));
                    }
                } else {// same number of total types
                    eventWeight = ((double) coalescingNodeStat + (double) mutatingNodeCount) / (double) activeNodeCount * (theta / (activeNodeCount - 1.0 + theta));
                    // this node now does not have any active options
                    availableNodes.remove(subtree.getNodeInWrappedTree(node));
                }
                int oldCoalescingCount = coalescingPairsCount();
                this.coalescingNodes.add(node);
                mutatingNodes.remove(node);
                int newCoalescingCount = coalescingPairsCount();

                return new StateChange(-1, newCoalescingCount - oldCoalescingCount, Math.log(eventWeight)); //just moved a node from 1 list to another their sum is the same
            }
        }

        public StateChange coalesceNodes(NodeRef left, double height) {
            assert coalescingPairsCount() >= 1;
            int oldCoalescingCount = coalescingPairsCount();
            //pick 2 nodes

            //update the weight
            double eventWeight = ((double) coalescingNodeStat + (double) mutatingNodeCount) / (double) coalescingNodes.size() * (coalescingNodes.size() - 1.0) / (activeNodeCount - 1.0 + theta);

            int node2;
            NodeRef parent;

            coalescingNodes.remove(left);
            node2 = random.nextInt(coalescingNodes.size());
            parent = nodeStack.pop();
            NodeRef right = coalescingNodes.get(node2);


            coalescingNodes.remove(right);

            activeNodeCount -= 2;
            availableNodes.remove(subtree.getNodeInWrappedTree(left));
            availableNodes.remove(subtree.getNodeInWrappedTree(right));


            subtree.addChildQuietly(parent, left);
            subtree.addChildQuietly(parent, right);

            subtree.setNodeHeightQuietly(parent, height);
            simulatedCoalescentEvents += 1;

            // add coalescent event in intervals
            activeIntervals.addCoalescentEvent(height,subtree.getNodeInWrappedTree(parent).getNumber());

            // At this point we need to update the state to reflect the
            // coalescent event
            if (subtree.isRoot(parent)) {
                assert coalescingNodes.size() == 0;
                complete = true;
                NodeRef baseNode = subtree.getNodeInWrappedTree(parent);
                if (!tree.isRoot(baseNode)) {
                    // this will also update the state but that's ok
                    activateNode(baseNode);
                }
            } else {
                if (coalescingNodes.size() == 1) {
                    availableNodes.remove(subtree.getNodeInWrappedTree(coalescingNodes.get(0)));
                }
                StateChange sc = activateSubtreeNode(parent);
                activeNodeCount += 1;
                assert sc.mutDelta == 0;
                // if the parent has mutations then it must be the root and will be activated above
            }
            int newCoalescingCount = coalescingPairsCount();
            return (new StateChange(0, newCoalescingCount - oldCoalescingCount, Math.log(eventWeight)));
        }

        //return the number of nodes current coalescing. If there is only 1 then it can't coalesce. return 0;
        public int coalescingNodeStat() {
            return Math.max(0, coalescingNodes.size() - 1);
        }
    }
}