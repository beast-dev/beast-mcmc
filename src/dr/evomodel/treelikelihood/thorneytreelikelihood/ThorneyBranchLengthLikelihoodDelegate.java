package dr.evomodel.treelikelihood.thorneytreelikelihood;

import dr.inference.model.Parameter;

public interface ThorneyBranchLengthLikelihoodDelegate  {
     double getLogLikelihood(double observed,double expected);

}
