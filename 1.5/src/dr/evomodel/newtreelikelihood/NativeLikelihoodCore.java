package dr.evomodel.newtreelikelihood;

import dr.evomodel.sitemodel.SiteModel;
import dr.evomodel.substmodel.SubstitutionModel;
import dr.app.beagle.evomodel.treelikelihood.AbstractTreeLikelihood;

import java.util.logging.Logger;

/*
 * NativeLikelihoodCore.java
 *
 * @author Andrew Rambaut
 *
 */

public class NativeLikelihoodCore implements LikelihoodCore {

    public static final String LIBRARY_NAME = "NativeLikelihoodCore";
    
    public NativeLikelihoodCore() { }

    public NativeLikelihoodCore(int stateCount) {
        StringBuffer sb = new StringBuffer();
        sb.append("Constructing native likelihood core\n");
        Logger.getLogger("dr.evomodel.treelikelihood").info(sb.toString());
    }

    public boolean canHandleTipPartials() {
        return true;
    }

    public boolean canHandleTipStates() {
        return true;
    }

    public boolean canHandleDynamicRescaling() {
    	return true;
    }
    
    public native void initialize(int nodeCount, int stateTipCount, int patternCount, int matrixCount);

    public void finalize() throws Throwable {
        super.finalize();
        freeNativeMemory();
    }

    private native void freeNativeMemory();

    public native void setTipPartials(int tipIndex, double[] partials);

    public native void setTipStates(int tipIndex, int[] states);

    public void updateSubstitutionModel(SubstitutionModel substitutionModel) {
        updateRootFrequencies(substitutionModel.getFrequencyModel().getFrequencies());
        updateEigenDecomposition(
                substitutionModel.getEigenVectors(),
                substitutionModel.getInverseEigenVectors(),
                substitutionModel.getEigenValues());
    }

    protected native void updateRootFrequencies(double[] frequencies);

    protected native void updateEigenDecomposition(double[][] eigenVectors,
                                                 double[][] inverseEigenValues,
                                                 double[] eigenValues);

    public void updateSiteModel(SiteModel siteModel) {
        if (rates == null) {
            rates = new double[siteModel.getCategoryCount()];
        }
        for (int i = 0; i < rates.length; i++) {
            rates[i] = siteModel.getRateForCategory(i);
        }
        updateCategoryRates(rates);
        updateCategoryProportions(siteModel.getCategoryProportions());
    }

    /**
     * A utility array to transfer category rates
     */
    private double[] rates = null;
    
    public void updatePartials(int[] operations, int[] dependencies, int operationCount, boolean rescale) {
    	updatePartials(operations, dependencies, operationCount);
    }
 
    protected native void updateCategoryRates(double[] rates);

    protected native void updateCategoryProportions(double[] proportions);

    public native void updateMatrices(int[] branchUpdateIndices, double[] branchLengths, int branchUpdateCount);

    public native void updatePartials(int[] operations, int[] dependencies, int operationCount);

    public native void calculateLogLikelihoods(int rootNodeIndex, double[] outLogLikelihoods);

    public native void storeState();

    public native void restoreState();

    /* Library loading routines */

    public static class LikelihoodCoreLoader implements LikelihoodCoreFactory.LikelihoodCoreLoader {

        public String getLibraryName() { return LIBRARY_NAME; }

        public LikelihoodCore createLikelihoodCore(int[] configuration, AbstractTreeLikelihood treeLikelihood) {
            int stateCount = configuration[0];
            try {
                System.loadLibrary(getLibraryName() + "-" + stateCount);
            } catch (UnsatisfiedLinkError e) {
                return null;
            }
            return new NativeLikelihoodCore(stateCount);
        }
    }

}