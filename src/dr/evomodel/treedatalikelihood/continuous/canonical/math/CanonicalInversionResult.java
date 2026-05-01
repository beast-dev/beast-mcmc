package dr.evomodel.treedatalikelihood.continuous.canonical.math;

/**
 * Result metadata for canonical flat-matrix inversions that may include
 * collapsed or missing dimensions.
 */
public final class CanonicalInversionResult {

    public enum Code {
        FULLY_OBSERVED,
        NOT_OBSERVED,
        PARTIALLY_OBSERVED
    }

    private final Code returnCode;
    private final int effectiveDimension;
    private final double logDeterminant;

    public CanonicalInversionResult(final Code returnCode,
                                    final int effectiveDimension,
                                    final double logDeterminant) {
        this.returnCode = returnCode;
        this.effectiveDimension = effectiveDimension;
        this.logDeterminant = logDeterminant;
    }

    public Code getReturnCode() {
        return returnCode;
    }

    public int getEffectiveDimension() {
        return effectiveDimension;
    }

    public double getLogDeterminant() {
        return logDeterminant;
    }

    public double getDeterminant() {
        return Math.exp(logDeterminant);
    }

    public static Code getCode(final int fullDimension, final int effectiveDimension) {
        if (effectiveDimension == fullDimension) {
            return Code.FULLY_OBSERVED;
        } else if (effectiveDimension == 0) {
            return Code.NOT_OBSERVED;
        }
        return Code.PARTIALLY_OBSERVED;
    }
}
