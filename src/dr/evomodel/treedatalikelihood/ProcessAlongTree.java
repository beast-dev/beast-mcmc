package dr.evomodel.treedatalikelihood;

import dr.evolution.tree.Tree;
import dr.evolution.tree.TreeTrait;
import dr.evolution.tree.TreeTraitProvider;
import dr.evomodel.branchratemodel.BranchRateModel;
import dr.inference.model.ModelListener;

public interface ProcessAlongTree extends TreeTraitProvider {

    Tree getTree();

    BranchRateModel getBranchRateModel();

    void calculatePostOrderStatistics();

    void addModelRestoreListener(ModelListener listener);

    void addModelListener(ModelListener listener);

    double getLogLikelihood();

    void addTrait(TreeTrait treeTrait);
}
