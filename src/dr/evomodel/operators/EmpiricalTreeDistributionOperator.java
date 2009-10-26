package dr.evomodel.operators;

import dr.evolution.tree.NodeRef;
import dr.evolution.tree.Tree;
import dr.evolution.util.Taxon;
import dr.math.MathUtils;
import dr.inference.operators.MCMCOperator;
import dr.inference.operators.OperatorFailedException;
import dr.inference.operators.SimpleMCMCOperator;
import dr.evomodel.tree.TreeModel;
import dr.evomodel.tree.EmpiricalTreeDistributionModel;

import java.util.Iterator;
import java.util.List;

/**
 * @author Andrew Rambaut
 * @version $Id$
 */
public class EmpiricalTreeDistributionOperator extends SimpleMCMCOperator {

    public final static String EMPIRICAL_TREE_DISTRIBUTION_OPERATOR = "empiricalTreeOperator";

    private final EmpiricalTreeDistributionModel treeModel;

    public EmpiricalTreeDistributionOperator(EmpiricalTreeDistributionModel treeModel, double weight) {
        this.treeModel = treeModel;
        setWeight(weight);
    }

    // IMPLEMENTATION: SimpleMCMCOperator
    // =========================================

    public double doOperation() throws OperatorFailedException {
        treeModel.drawTreeIndex();
        return 0;
    }

    public String getPerformanceSuggestion() {
        return "";
    }

    public String getOperatorName() {
        return EMPIRICAL_TREE_DISTRIBUTION_OPERATOR;
    }

}