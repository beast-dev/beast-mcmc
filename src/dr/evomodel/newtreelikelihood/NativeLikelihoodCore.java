package dr.evomodel.newtreelikelihood;

import dr.evomodel.sitemodel.SiteModel;
import dr.evomodel.substmodel.SubstitutionModel;

/*
 * NativeLikelihoodCore.java
 *
 * @author Andrew Rambaut
 *
 */

public class NativeLikelihoodCore implements LikelihoodCore {

    public static final String LIBRARY_NAME = "NativeLikelihoodCore";

    public NativeLikelihoodCore() {
        // don't need to do any thing;
    }

    public native void initialize(int nodeCount, int patternCount, int matrixCount);

    public void finalize() throws Throwable {
        super.finalize();
        freeNativeMemory();
    }

    private native void freeNativeMemory();

    public native void setTipPartials(int tipIndex, double[] partials);

    public native void updateSubstitutionModel(SubstitutionModel substitutionModel);

    public native void updateSiteModel(SiteModel siteModel);

    public native void updateMatrices(int[] branchUpdateIndices, double[] branchLengths, int branchUpdateCount);

    public native void updatePartials(int[] operations, int[] dependencies, int operationCount);

    public native void calculateLogLikelihoods(int rootNodeIndex, double[] outLogLikelihoods);

    public native void storeState();

    public native void restoreState();

    /* Library loading routines */

    public static class LikelihoodCoreLoader implements LikelihoodCoreFactory.LikelihoodCoreLoader {

        public String getLibraryName() { return "NativeLikelihoodCore"; }

        public LikelihoodCore createLikelihoodCore(int[] configuration, AbstractTreeLikelihood treeLikelihood) {
            int stateCount = configuration[0];
            try {
                System.loadLibrary(getLibraryName()+"-"+stateCount);
            } catch (UnsatisfiedLinkError e) {
                return null;
            }
            return new NativeLikelihoodCore();
        }
    }

}