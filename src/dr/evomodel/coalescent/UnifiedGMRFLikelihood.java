package dr.evomodel.coalescent;

import java.util.List;

import dr.evolution.coalescent.IntervalEventList;
import dr.evolution.coalescent.IntervalList;
import dr.evolution.coalescent.TreeIntervalList;
import dr.evolution.tree.Tree;
import dr.evomodel.tree.TreeModel;
import dr.inference.model.Likelihood;
import dr.inference.model.MatrixParameter;
import dr.inference.model.Parameter;
import dr.util.Citation;
import no.uib.cipr.matrix.SymmTridiagMatrix;

public interface UnifiedGMRFLikelihood extends CoalescentIntervalProvider, Likelihood{
    interface Skyride extends UnifiedGMRFLikelihood{

    double[] getCoalescentIntervals();

    boolean isCoalescentInterval(int interval);

    void initializationReport();

    double[] getSufficientStatistics();

    String toString();

    IntervalNodeMapProvider getIntervalNodeMapping();

    SymmTridiagMatrix getScaledWeightMatrix(double precision);

    SymmTridiagMatrix getStoredScaledWeightMatrix(double precision);

    SymmTridiagMatrix getScaledWeightMatrix(double precision, double lambda);

    int getCoalescentIntervalDimension();

    double getCoalescentInterval(int i);

    void setupCoalescentIntervals();

    double[] getCoalescentIntervalHeights();

    SymmTridiagMatrix getCopyWeightMatrix();

    SymmTridiagMatrix getStoredScaledWeightMatrix(double precision, double lambda);

    Parameter getPrecisionParameter();

    Parameter getPopSizeParameter();

    Parameter getLambdaParameter();

    SymmTridiagMatrix getWeightMatrix();

    Parameter getBetaParameter();

    MatrixParameter getDesignMatrix();

    double calculateWeightedSSE();

    int getIntervalCount();

    int getLineageCount(int i);

}
    
    interface SkyGrid extends Skyride{

        double getCoalescentEventsStatisticValue(int i);
    
        double getLogLikelihood();
    
        int nLoci();
        
        double[] getNumCoalEvents();
        Tree getTree(int nt);
    
        IntervalEventList getTreeIntervals(int nt);
    
        double getPopulationFactor(int nt);
    
        @Deprecated
        List<Parameter> getBetaListParameter();
    
        List<MatrixParameter> getCovariates();
    
        Parameter getParameter();
    
        int getDimension();
    
        double[] getGradientWrtLogPopulationSize();
    
        double[] getDiagonalHessianWrtLogPopulationSize();
    
        double[] getGradientWrtPrecision();
    
        double[] getDiagonalHessianWrtPrecision();
    
        double[] getGradientWrtRegressionCoefficients();
    
        double[] getDiagonalHessianWrtRegressionCoefficients();
    
        double[] getGridPoints();
    
        
    }
    
}
