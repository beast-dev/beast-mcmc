package dr.evomodel.treelikelihood.thorneytreelikelihood;

import dr.evolution.tree.NodeRef;
import dr.evolution.tree.Tree;


public interface ThorneyBranchLengthLikelihoodDelegate  {
     double getLogLikelihood(double mutations, Tree tree, NodeRef node);

     double getGradientWrtTime(double mutations, double time);
}