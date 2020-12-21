package dr.evomodel.treelikelihood.thorneytreelikelihood;

import dr.inference.model.Model;

public interface ThorneyBranchLengthLikelihoodDelegate extends Model {
     double getLogLikelihood(double mutations,double time);

     double getGradientWrtTime(double mutations, double time);
}
