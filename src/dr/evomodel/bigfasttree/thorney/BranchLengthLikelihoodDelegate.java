package dr.evomodel.bigfasttree.thorney;


public interface BranchLengthLikelihoodDelegate  {
     double getLogLikelihood(MutationList mutations, double branchLength);
     public double getGradientWrtTime(MutationList mutations, double time, double branchRate);
}