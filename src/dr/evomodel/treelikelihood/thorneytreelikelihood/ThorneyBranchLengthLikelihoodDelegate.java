package dr.evomodel.treelikelihood.thorneytreelikelihood;

import dr.evolution.tree.NodeRef;
import dr.evolution.tree.Tree;


public interface ThorneyBranchLengthLikelihoodDelegate  {
     double getLogLikelihood(double observed, double expected, Tree tree, NodeRef node);

}
