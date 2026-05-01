package dr.evomodel.treedatalikelihood.continuous.canonical.branch;

import dr.evomodel.treedatalikelihood.continuous.CanonicalDebugOptions;
import dr.evomodel.treedatalikelihood.continuous.CanonicalGradientFallbackPolicy;
import dr.evomodel.treedatalikelihood.continuous.OUGaussianBranchTransitionProvider;

import dr.evomodel.treedatalikelihood.continuous.canonical.message.CanonicalNumerics;
import org.ejml.data.DenseMatrix64F;

final class OUCanonicalBranchSpdDebugContext implements CanonicalNumerics.DenseSpdFailureDump {

    private final DenseMatrix64F actualizationMatrix;
    private final DenseMatrix64F displacementVector;
    private final DenseMatrix64F branchPrecisionMatrix;
    private final DenseMatrix64F branchCovarianceMatrix;
    private final DenseMatrix64F aboveChildPrecisionMatrix;
    private final DenseMatrix64F aboveChildCovarianceMatrix;
    private final DenseMatrix64F aboveParentPrecisionMatrix;
    private final DenseMatrix64F aboveParentCovarianceMatrix;

    OUCanonicalBranchSpdDebugContext(final DenseMatrix64F actualizationMatrix,
                                     final DenseMatrix64F displacementVector,
                                     final DenseMatrix64F branchPrecisionMatrix,
                                     final DenseMatrix64F branchCovarianceMatrix,
                                     final DenseMatrix64F aboveChildPrecisionMatrix,
                                     final DenseMatrix64F aboveChildCovarianceMatrix,
                                     final DenseMatrix64F aboveParentPrecisionMatrix,
                                     final DenseMatrix64F aboveParentCovarianceMatrix) {
        this.actualizationMatrix = actualizationMatrix;
        this.displacementVector = displacementVector;
        this.branchPrecisionMatrix = branchPrecisionMatrix;
        this.branchCovarianceMatrix = branchCovarianceMatrix;
        this.aboveChildPrecisionMatrix = aboveChildPrecisionMatrix;
        this.aboveChildCovarianceMatrix = aboveChildCovarianceMatrix;
        this.aboveParentPrecisionMatrix = aboveParentPrecisionMatrix;
        this.aboveParentCovarianceMatrix = aboveParentCovarianceMatrix;
    }

    @Override
    public void appendJson(final StringBuilder sb,
                           final String context,
                           final DenseMatrix64F originalSource,
                           final DenseMatrix64F symmetrizedSource,
                           final double jitterBase) {
        sb.append(",\n");
        sb.append("\"actualization\":").append(CanonicalNumerics.jsonMatrix(actualizationMatrix)).append(",\n");
        sb.append("\"displacement\":").append(CanonicalNumerics.jsonVector(displacementVector)).append(",\n");
        sb.append("\"branchPrecisionBuffer\":").append(CanonicalNumerics.jsonMatrix(branchPrecisionMatrix)).append(",\n");
        sb.append("\"branchCovarianceBuffer\":").append(CanonicalNumerics.jsonMatrix(branchCovarianceMatrix)).append(",\n");
        sb.append("\"aboveChildPrecisionBuffer\":").append(CanonicalNumerics.jsonMatrix(aboveChildPrecisionMatrix)).append(",\n");
        sb.append("\"aboveChildCovarianceBuffer\":").append(CanonicalNumerics.jsonMatrix(aboveChildCovarianceMatrix)).append(",\n");
        sb.append("\"aboveParentPrecisionBuffer\":").append(CanonicalNumerics.jsonMatrix(aboveParentPrecisionMatrix)).append(",\n");
        sb.append("\"aboveParentCovarianceBuffer\":").append(CanonicalNumerics.jsonMatrix(aboveParentCovarianceMatrix)).append("\n");
    }
}
