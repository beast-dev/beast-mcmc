package dr.evomodel.bigfasttree.thorney;


public interface BranchLengthLikelihoodDelegate  {
     double getLogLikelihood(double mutations, double branchLength);

     double getGradientWrtTime(double mutations, double time);
}