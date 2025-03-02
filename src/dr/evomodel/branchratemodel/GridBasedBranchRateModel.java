package dr.evomodel.branchratemodel;

import dr.evolution.tree.NodeRef;
import dr.evolution.tree.Tree;
import dr.evomodel.bigfasttree.BigFastTreeIntervals;
import dr.evomodel.tree.TreeModel;
import dr.inference.loggers.Loggable;
import dr.inference.model.Model;
import dr.inference.model.Parameter;
import dr.inference.model.Variable;
import dr.xml.Reportable;

import java.util.Arrays;
import java.util.Comparator;

public class GridBasedBranchRateModel extends AbstractBranchRateModel implements Reportable, Loggable {
    private final TreeModel tree;
    private final Parameter gridPoints;
    private final Parameter levelSpecificRates;
    private double[] branchesIntersections; // TODO make this a vector
    private double[] branchRates;
    private boolean ratesKnown;
    private boolean nodesOrderKnown;
    private Integer[] orderedNodesIndexes;

    public GridBasedBranchRateModel(TreeModel tree, Parameter gridPoints, Parameter levelSpecificRates) {
        super("gridBasedBranchRateModel");
        this.tree = tree;
        this.gridPoints = gridPoints;
        this.levelSpecificRates = levelSpecificRates;
        this.ratesKnown = false;
        this.nodesOrderKnown = false;
        this.branchRates = new double[tree.getNodeCount() - 1];
        this.branchesIntersections = new double[(gridPoints.getDimension() + 1) * (tree.getNodeCount() - 1)];
        this.orderedNodesIndexes = new Integer[tree.getNodeCount()];

        addModel(tree);
        addVariable(gridPoints);
        addVariable(levelSpecificRates);

        buildBranchRates();
    }

    private void buildBranchRates() {
        if (!ratesKnown) {
            computeBranchRates();
            ratesKnown = true;
        }
    }

    private void computeBranchRates() {
        computeIntersectionsMatrix();

        for (int nodeIndex = 0; nodeIndex < tree.getNodeCount() - 1; nodeIndex++) {
            double branchRate = 0;
            for (int gridIndex = 0; gridIndex < gridPoints.getDimension() + 1; gridIndex++) {
//                branchRate1 +=  branchesIntersectionsMatrix[gridIndex][nodeIndex];
                branchRate += branchesIntersections[gridIndex + (gridPoints.getDimension() + 1) * nodeIndex] * levelSpecificRates.getParameterValue(gridIndex);
            }
            branchRates[nodeIndex] = branchRate;
        }

    }

    private void orderNodesByHeight() {
        if (!nodesOrderKnown) {
            for (int i = 0; i < tree.getNodeCount(); i++) {
                orderedNodesIndexes[i] = i;
            }
            Arrays.sort(orderedNodesIndexes, Comparator.comparingDouble( i -> tree.getNodeHeight(tree.getNode((i)))));
            nodesOrderKnown = true;
        }
    }

    private void computeIntersectionsMatrix() { // TODO check that grid points and node heights are on the same relative scale
        int minGridIndex = 0;
        orderNodesByHeight();
        for (int i = 0; i < tree.getNodeCount() - 1; i++) {
            int nodeIndex = orderedNodesIndexes[i];
            double tChild = tree.getNodeHeight(tree.getNode(nodeIndex));
            double tParent = tree.getNodeHeight(tree.getParent(tree.getNode(nodeIndex)));
            double tLower = tChild;

            while (minGridIndex < gridPoints.getDimension() && gridPoints.getParameterValue(minGridIndex) < tChild) {
                minGridIndex ++;
            }
            int gridIndex = minGridIndex;

            if (gridIndex < gridPoints.getDimension() && gridPoints.getParameterValue(gridIndex) < tParent) {
                while (gridIndex < gridPoints.getDimension() && gridPoints.getParameterValue(gridIndex) < tParent) {
                    branchesIntersections[gridIndex + (gridPoints.getDimension() + 1) * nodeIndex] = gridPoints.getParameterValue(gridIndex) - tLower;
                    tLower = gridPoints.getParameterValue(gridIndex);
                    gridIndex++;
                }
                branchesIntersections[gridIndex + (gridPoints.getDimension() + 1) * nodeIndex] = tParent - tLower;
            } else {
                branchesIntersections[gridIndex + (gridPoints.getDimension() + 1) * nodeIndex] = tParent - tChild;
            }
        }
    }

    @Override
    public double getBranchRate(Tree tree, NodeRef node) { // TODO I am not really using the "tree" but I should I guess
        buildBranchRates();
        return branchRates[node.getNumber()];
    }

    @Override
    protected void handleModelChangedEvent(Model model, Object object, int index) {
        if (model == tree) {
            ratesKnown = false;
            nodesOrderKnown = false;
        }
    }

    @Override
    protected void storeState() {

    }

    @Override
    protected void restoreState() {

    }

    @Override
    protected void acceptState() {

    }

    @Override
    protected void handleVariableChangedEvent(Variable variable, int index, Parameter.ChangeType type) {
        ratesKnown = false;
    }

    @Override
    public String getReport() {
        return "Branches intersections matrix: " + Arrays.toString(branchesIntersections) + "\n"
                + "Branch rates: " + Arrays.toString(branchRates);
    }
}
