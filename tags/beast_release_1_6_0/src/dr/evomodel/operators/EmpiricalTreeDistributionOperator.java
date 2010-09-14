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

    public static final String EMPIRICAL_TREE_DISTRIBUTION_OPERATOR = "empiricalTreeDistributionOperator";
    public static final String METROPOLIS_HASTINGS = "metropolisHastings";

    private final EmpiricalTreeDistributionModel treeModel;
    private final boolean metropolisHastings;

    public EmpiricalTreeDistributionOperator(EmpiricalTreeDistributionModel treeModel, boolean metropolisHastings, double weight) {
        this.treeModel = treeModel;
        setWeight(weight);
        this.metropolisHastings = metropolisHastings;
    }

    // IMPLEMENTATION: SimpleMCMCOperator
    // =========================================

    public double doOperation() throws OperatorFailedException {
        treeModel.drawTreeIndex();

        if (metropolisHastings) {           
            return 0.0;
        }

        // if not MetropolisHastings, always accept these moves...
        return Double.POSITIVE_INFINITY;
    }

    public String getPerformanceSuggestion() {
        return "";
    }

    public String getOperatorName() {
        return EMPIRICAL_TREE_DISTRIBUTION_OPERATOR;
    }

}