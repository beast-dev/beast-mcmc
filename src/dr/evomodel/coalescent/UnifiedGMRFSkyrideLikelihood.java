package dr.evomodel.coalescent;

import dr.inference.model.Likelihood;
import dr.inference.model.Parameter;
import no.uib.cipr.matrix.SymmTridiagMatrix;

public interface UnifiedGMRFSkyrideLikelihood extends CoalescentIntervalProvider,Likelihood{
    IntervalNodeMapProvider getIntervalNodeMapping();
    int  getCoalescentIntervalDimension();
    Parameter getPopSizeParameter();
    Parameter getPrecisionParameter();
    Parameter getLambdaParameter();

    double[] getSufficientStatistics();
    SymmTridiagMatrix getStoredScaledWeightMatrix(double currentPrecision, double currentLambda);
    SymmTridiagMatrix getScaledWeightMatrix(double proposedPrecision, double proposedLambda);

    //Interval stuff 
    int getIntervalCount();
    boolean isCoalescentInterval(int interval); // because interval types are linked to classes for now
    int getLineageCount(int interval);
}
