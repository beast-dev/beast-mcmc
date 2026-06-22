package dr.evomodel.continuous.ou.diagonal;

import dr.evomodel.continuous.ou.kernel.OUPreparedBranch;

public final class DiagonalOUPreparedBranch implements OUPreparedBranch {

    private final double[] actualizationDiagonal;
    private final double[] displacement;
    private final double[] precision;
    private double effectiveBranchLength;
    private double logNormalizer;
    private boolean valid;

    DiagonalOUPreparedBranch(final int dimension) {
        this.actualizationDiagonal = new double[dimension];
        this.displacement = new double[dimension];
        this.precision = new double[dimension * dimension];
        this.effectiveBranchLength = Double.NaN;
        this.logNormalizer = Double.NaN;
        this.valid = false;
    }

    public double[] actualizationDiagonal() {
        return actualizationDiagonal;
    }

    public double[] displacement() {
        return displacement;
    }

    public double[] precision() {
        return precision;
    }

    public double logNormalizer() {
        return logNormalizer;
    }

    public boolean isCurrent(final double branchLength) {
        return valid
                && Double.doubleToLongBits(effectiveBranchLength)
                == Double.doubleToLongBits(branchLength);
    }

    @Override
    public void markInvalid() {
        valid = false;
    }

    void markCurrent(final double branchLength, final double logNormalizer) {
        this.effectiveBranchLength = branchLength;
        this.logNormalizer = logNormalizer;
        this.valid = true;
    }
}
