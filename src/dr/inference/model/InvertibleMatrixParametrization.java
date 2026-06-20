package dr.inference.model;

import java.util.List;

/**
 * Base class for square matrix parameters that are nonsingular by construction.
 */
public abstract class InvertibleMatrixParametrization extends AbstractComputedCompoundMatrix {

    private final Parameter[] nativeParameters;
    private final int[] nativeParameterOffsets;
    private final CompoundParameter nativeParameter;
    private final int nativeDimension;
    private final double[] nativeGradientScratch;

    protected InvertibleMatrixParametrization(final String name,
                                             final int dimension,
                                             final List<Parameter> parameters) {
        super(name, checkDimension(dimension), parameters);
        this.nativeParameters = parameters.toArray(new Parameter[parameters.size()]);
        this.nativeParameterOffsets = new int[nativeParameters.length];
        this.nativeParameter = new CompoundParameter(name + ".native");

        int offset = 0;
        for (int i = 0; i < nativeParameters.length; i++) {
            nativeParameterOffsets[i] = offset;
            offset += nativeParameters[i].getDimension();
            nativeParameter.addParameter(nativeParameters[i]);
        }
        this.nativeDimension = offset;
        this.nativeGradientScratch = new double[nativeDimension];
    }

    public final boolean isInvertible() {
        return true;
    }

    public abstract double getLogAbsDeterminant();

    public abstract double getDeterminantSign();

    public double getDeterminant() {
        return getDeterminantSign() * Math.exp(getLogAbsDeterminant());
    }

    public abstract void fillInverse(double[] inverseRowMajor);

    public final CompoundParameter getNativeParameter() {
        return nativeParameter;
    }

    public final int getNativeDimension() {
        return nativeDimension;
    }

    public final int getNativeParameterCount() {
        return nativeParameters.length;
    }

    public final Parameter getNativeSubParameter(final int index) {
        return nativeParameters[index];
    }

    public final boolean supportsNativeParameter(final Parameter parameter) {
        return parameter == nativeParameter
                || findNativeSubParameterIndex(parameter) >= 0
                || matchesNativeCompound(parameter);
    }

    public final int getPullBackGradientDimension(final Parameter parameter) {
        if (parameter == nativeParameter || matchesNativeCompound(parameter)) {
            return nativeDimension;
        }

        final int nativeIndex = findNativeSubParameterIndex(parameter);
        if (nativeIndex >= 0) {
            return nativeParameters[nativeIndex].getDimension();
        }

        throw unsupportedNativeParameter(parameter);
    }

    public final double[] pullBackGradientForParameter(final Parameter parameter,
                                                       final double[] gradientWrtMatrixRowMajor,
                                                       final int dimension) {
        final double[] gradient = new double[getPullBackGradientDimension(parameter)];
        fillPullBackGradientForParameter(parameter, gradientWrtMatrixRowMajor, dimension, gradient, 0);
        return gradient;
    }

    public final synchronized void fillPullBackGradientForParameter(final Parameter parameter,
                                                                    final double[] gradientWrtMatrixRowMajor,
                                                                    final int dimension,
                                                                    final double[] out,
                                                                    final int offset) {
        final int outDimension = getPullBackGradientDimension(parameter);
        if (out == null || out.length < offset + outDimension) {
            throw new IllegalArgumentException(
                    "out must have room for " + outDimension + " entries at offset " + offset);
        }
        if (parameter == nativeParameter || matchesNativeCompound(parameter)) {
            fillPullBackGradientFlat(gradientWrtMatrixRowMajor, dimension, out, offset);
            return;
        }

        final int nativeIndex = findNativeSubParameterIndex(parameter);
        if (nativeIndex >= 0) {
            final int length = nativeParameters[nativeIndex].getDimension();
            fillPullBackGradientFlat(gradientWrtMatrixRowMajor, dimension, nativeGradientScratch, 0);
            System.arraycopy(nativeGradientScratch, nativeParameterOffsets[nativeIndex], out, offset, length);
            return;
        }

        throw unsupportedNativeParameter(parameter);
    }

    public final void fillPullBackGradientFlat(final double[] gradientWrtMatrixRowMajor,
                                               final int dimension,
                                               final double[] out) {
        fillPullBackGradientFlat(gradientWrtMatrixRowMajor, dimension, out, 0);
    }

    public abstract void fillPullBackGradientFlat(double[] gradientWrtMatrixRowMajor,
                                                  int dimension,
                                                  double[] out,
                                                  int offset);

    private int findNativeSubParameterIndex(final Parameter parameter) {
        for (int i = 0; i < nativeParameters.length; i++) {
            if (parameter == nativeParameters[i]) {
                return i;
            }
        }
        return -1;
    }

    private boolean matchesNativeCompound(final Parameter parameter) {
        if (!(parameter instanceof CompoundParameter)) {
            return false;
        }
        final CompoundParameter compound = (CompoundParameter) parameter;
        if (compound.getParameterCount() != nativeParameters.length) {
            return false;
        }
        for (int i = 0; i < nativeParameters.length; i++) {
            if (compound.getParameter(i) != nativeParameters[i]) {
                return false;
            }
        }
        return true;
    }

    private IllegalArgumentException unsupportedNativeParameter(final Parameter parameter) {
        return new IllegalArgumentException("Unsupported native parameter: "
                + (parameter == null ? "null" : parameter.getId()));
    }

    private static int checkDimension(final int dimension) {
        if (dimension < 1) {
            throw new IllegalArgumentException("dimension must be positive");
        }
        return dimension;
    }
}
