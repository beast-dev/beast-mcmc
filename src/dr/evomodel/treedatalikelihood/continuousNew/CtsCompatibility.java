package dr.evomodel.treedatalikelihood.continuousNew;

import dr.evomodel.treedatalikelihood.DataLikelihoodDelegate;
import dr.evomodel.treedatalikelihood.continuous.ContinuousRateTransformation;
import dr.evomodel.treedatalikelihood.continuous.ContinuousTraitPartialsProvider;
import dr.evomodel.treedatalikelihood.continuous.DiffusionProcessDelegate;
import dr.evomodel.treedatalikelihood.continuous.RootProcessDelegate;
import dr.evomodel.treedatalikelihood.continuous.cdi.ContinuousDiffusionIntegrator;
import dr.evomodel.treedatalikelihood.continuous.cdi.PrecisionType;
import dr.evomodel.treedatalikelihood.continuousNew.canonical.adapter.CanonicalOUGradientAdapter;

public interface CtsCompatibility extends DataLikelihoodDelegate {
    //WishartStatisticsWrapper is the only one that actually required this to be an abstractModel
    void enableCanonicalLikelihood();
    void enableCdiBMLikelihood();
    void enableMomentBMLikelihood();
    void enableCanonicalOULikelihood();
    PrecisionType getPrecisionType();
    void setExtensionHelper();
    int getTraitDim();
    int getTraitCount();
    RootProcessDelegate getRootProcessDelegate();
    int getActiveNodeIndex(int nodeIndex);
    int getActiveMatrixIndex(int nodeIndex);
    int getActivePrecisionIndex(int precisionIndex);
    DiffusionProcessDelegate getDiffusionProcessDelegate();
    ContinuousDiffusionIntegrator getIntegrator();
    dr.evomodel.continuous.MultivariateDiffusionModel getDiffusionModel();
    int getPartialBufferCount();
    ContinuousTraitPartialsProvider getDataModel();
    void getPostOrderPartial(int nodeIndex, double[] buffer);
    void addBranchConditionalDensityTrait(String traitName);
    boolean usesCanonicalOULikelihood();
    boolean usesCdiBMLikelihood();
    boolean usesMomentBMLikelihood();
    boolean hasCanonicalPostOrderMessageProvider();
    boolean hasCanonicalPreOrderMessageProvider();
    CanonicalOUGradientAdapter getCanonicalOUGradientAdapter();
    ContinuousRateTransformation getRateTransformation();
}
