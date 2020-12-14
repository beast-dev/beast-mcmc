package dr.evomodel.bigfasttree;

import dr.evolution.distance.DistanceMatrix;
import dr.evolution.tree.NodeRef;
import dr.evolution.tree.Tree;
import dr.evomodel.tree.TreeChangedEvent;
import dr.evomodel.tree.TreeModel;
import dr.inference.model.AbstractModel;
import dr.inference.model.Model;
import dr.inference.model.Parameter;
import dr.inference.model.Variable;

import java.util.HashMap;
import java.util.HashSet;
import java.util.Map;
import java.util.Set;

/**
 * A branch length provider in which the branch lengths are calculated from a distance matrix
 * using the methods of Rzhetsky & Nei (1993) MBE.
 *
 * @author Andrew Rambaut
 * @author JT McCrone
 */
public class RzhetskyNeiBranchLengthProvider extends AbstractModel implements BranchLengthProvider {
    public static String RZHETSKY_NEI_BRANCH_LENGTH_PROVIDER = "RzhetskyNeiBranchLengthProvider";

    public RzhetskyNeiBranchLengthProvider( DistanceMatrix distanceMatrix, TreeModel tree) {
        super(RZHETSKY_NEI_BRANCH_LENGTH_PROVIDER);

        this.distanceMatrix = distanceMatrix;
        this.tree = tree;

        branchLengths = new double[tree.getNodeCount()];
        stroredBranchLengths = new double[tree.getNodeCount()];

        distanceSums = new double[tree.getNodeCount()];
        storedDistanceSums = new double[tree.getNodeCount()];

        updateNode = new boolean[tree.getNodeCount()];
        storedUpdatedNodes = new boolean[tree.getNodeCount()];

        allTaxonSet = new HashSet<>(getTaxonSets(this.tree, this.tree.getRoot()));
    }

    private Set<Integer> getTaxonSets(Tree tree, NodeRef node) {

        Set<Integer> taxonSet = new HashSet<>();
        if (tree.isExternal(node)) {
            taxonSet.add(node.getNumber());
        } else {
            assert tree.getChildCount(node) == 2 : "Must be a strictly bifurcating tree";

            for (int i = 0; i < tree.getChildCount(node); i++) {
                taxonSet.addAll(getTaxonSets(tree, tree.getChild(node, i)));
            }
        }

        taxonSetMap.put(node.getNumber(), taxonSet);
        return taxonSet;
    }

    private void calculateBranchLengths(NodeRef node, NodeRef sibling) {

        double length =-1;

        if (tree.isExternal(node)) {

            Set<Integer> taxonSetC = taxonSetMap.get(node.getNumber());
            Set<Integer> taxonSetB = taxonSetMap.get(sibling.getNumber());

            Set<Integer> taxonSetA = new HashSet<>(allTaxonSet);
            taxonSetA.removeAll(taxonSetC);
            taxonSetA.removeAll(taxonSetB);

            double nA = taxonSetA.size();
            double nB = taxonSetB.size();

            double dCA = getSumOfDistances(taxonSetC, taxonSetA);
            double dCB = getSumOfDistances(taxonSetC, taxonSetB);
            double dAB = getSumOfDistances(taxonSetA, taxonSetB);

            // Equation 4 from R&N1993
            length = 0.5 * ((dCA / nA) + (dCB / nB) - (dAB / (nA * nB)));

        } else {

            NodeRef child1 = tree.getChild(node, 0);
            NodeRef child2 = tree.getChild(node, 1);

            calculateBranchLengths(child1, child2);

            calculateBranchLengths(child2, child1);
            if(node!=tree.getRoot()){
                Set<Integer> taxonSetC = taxonSetMap.get(child1.getNumber());
                Set<Integer> taxonSetD = taxonSetMap.get(child2.getNumber());

                Set<Integer> taxonSetB = taxonSetMap.get(sibling.getNumber());

                Set<Integer> taxonSetA = new HashSet<>(allTaxonSet);
                taxonSetA.removeAll(taxonSetC);
                taxonSetA.removeAll(taxonSetD);
                taxonSetA.removeAll(taxonSetB);

                double nA = taxonSetA.size();
                double nB = taxonSetB.size();
                double nC = taxonSetC.size();
                double nD = taxonSetD.size();

                // Equation 3 from R&N1993
                double gamma = (nB * nC + nA * nD) / ((nA + nB) * (nC + nD));

                double dAC = getSumOfDistances(taxonSetA, taxonSetC);
                double dBD = getSumOfDistances(taxonSetB, taxonSetD);
                double dBC = getSumOfDistances(taxonSetB, taxonSetC);
                double dAD = getSumOfDistances(taxonSetA, taxonSetD);
                double dAB = getSumOfDistances(taxonSetA, taxonSetB);
                double dCD = getSumOfDistances(taxonSetC, taxonSetD);

                // Equation 2 from R&N1993
                length = 0.5 * (
                        gamma * ((dAC / nA * nC) + (dBD / nB * nD)) +
                                (1.0 - gamma) * ((dBC / nB * nC) + (dAD / nA * nD)) -
                                (dAB / nA * nB) -
                                (dCD / nC * nD)
                );
            }

        }

        branchLengths[node.getNumber()] = length;
    }


    private double getSumOfDistances(Set<Integer> taxonSet1, Set<Integer> taxonSet2) {
        // This will be slow to do every time.
        // TODO do this in
        double sum = 0.0;
        for (int i : taxonSet1) {
            for (int j : taxonSet2) {
                sum += distanceMatrix.getElement(i, j);
            }
        }

        return sum;
    }

    @Override
    public double getBranchLength(Tree tree, NodeRef node) {
//        if (!branchLengthsKnown) {
            calculateBranchLengths(tree.getRoot(), null);
//        }
        return branchLengths[node.getNumber()];
    }

    /**
     * Set update flag for node and remove it's old contribution to the likelihood.
     * Also handle the root and children so that the 1 branch between children is marked as updated.
     *
     * @param node
     */
    protected void updateNode(NodeRef node) {

        updateNode[node.getNumber()] = true;
        NodeRef parent = tree.getParent(node);
        if (parent != null && !updateNode[parent.getNumber()]) {
            updateNode(parent);
        }
        branchLengthsKnown = false;
    }

    /**
     * Set update flag for a node and its direct children
     */
    protected void updateNodeAndChildren(NodeRef node) {

        updateNode(node);

        for (int i = 0; i < tree.getChildCount(node); i++) {
            NodeRef child = tree.getChild(node, i);
            updateNode(child);
        }

        branchLengthsKnown = false;
    }

    /**
     * Set update flag for all nodes
     */
    protected void updateAllNodes() {
        for (int i = 0; i < tree.getNodeCount(); i++) {
            updateNode[i] = true;
        }
        branchLengthsKnown = false;
    }


    @Override
    protected void handleModelChangedEvent(Model model, Object object, int index) {
        fireModelChanged();

        if (model == tree) {
            if (object instanceof TreeChangedEvent) {

                if (((TreeChangedEvent) object).isNodeChanged()) {
                    // If a node event occurs the node and its two child nodes
                    // are flagged for updating (this will result in everything
                    // above being updated as well. Node events occur when a node
                    // is added to a branch, removed from a branch or its height or
                    // rate changes.
                    NodeRef node = ((TreeChangedEvent) object).getNode();
                    updateNodeAndChildren(node);

                } else if (((TreeChangedEvent) object).isTreeChanged()) {
                    // Full tree events result in a complete updating of the tree likelihood
                    // This event type is now used for EmpiricalTreeDistributions.
//                    System.err.println("Full tree update event - these events currently aren't used\n" +
//                            "so either this is in error or a new feature is using them so remove this message.");
                    updateAllNodes();
                } else {
                    // Other event types are ignored (probably trait changes).
                    //System.err.println("Another tree event has occured (possibly a trait change).");
                }
            }
        }
    }

    /**
     * Additional state information, outside of the sub-model is stored by this call.
     */
    @Override
    protected void storeState() {

    }

    /**
     * After this call the model is guaranteed to have returned its extra state information to
     * the values coinciding with the last storeState call.
     * Sub-models are handled automatically and do not need to be considered in this method.
     */
    @Override
    protected void restoreState() {

    }

    /**
     * This call specifies that the current state is accept. Most models will not need to do anything.
     * Sub-models are handled automatically and do not need to be considered in this method.
     */
    @Override
    protected void acceptState() {

    }

    /**
     * This method is called whenever a parameter is changed.
     * <p/>
     * It is strongly recommended that the model component sets a "dirty" flag and does no
     * further calculations. Recalculation is typically done when the model component is asked for
     * some information that requires them. This mechanism is 'lazy' so that this method
     * can be safely called multiple times with minimal computational cost.
     *
     * @param variable
     * @param index
     * @param type
     */
    @Override
    protected void handleVariableChangedEvent(Variable variable, int index, Parameter.ChangeType type) {

    }


    private final DistanceMatrix distanceMatrix;
    private final Set<Integer> allTaxonSet;
    private boolean branchLengthsKnown;

    private double[] distanceSums;
    private double[] storedDistanceSums;


    private double[] branchLengths;
    private double[] stroredBranchLengths;

    private boolean[] updateNode;
    private boolean[] storedUpdatedNodes;

    private final TreeModel tree;
    private final Map<Integer, Set<Integer>> taxonSetMap = new HashMap<>();
}

