package dr.evomodel.treedatalikelihood;

import dr.evolution.tree.Tree;
import dr.evomodel.branchratemodel.BranchRateModel;
import dr.inference.model.ModelListener;

public interface ProcessAlongTree {

    Tree getTree();

    BranchRateModel getBranchRateModel();

    void calculatePostOrderStatistics();

    void addModelRestoreListener(ModelListener listener);

    void addModelListener(ModelListener listener);

    double getLogLikelihood();
}
