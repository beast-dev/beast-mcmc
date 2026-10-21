package dr.inference.model;

/**
 * Matrix-valued parameters that expose a native optimization coordinate system
 * can pull dense matrix-entry gradients back to those native coordinates.
 */
public interface NativeMatrixGradientProvider {

    Parameter getNativeParameter();

    boolean supportsNativeParameter(Parameter parameter);

    int getPullBackGradientDimension(Parameter parameter);

    double[] pullBackGradientForParameter(Parameter parameter,
                                          double[] gradientWrtMatrixRowMajor,
                                          int dimension);

    void fillPullBackGradientForParameter(Parameter parameter,
                                          double[] gradientWrtMatrixRowMajor,
                                          int dimension,
                                          double[] out,
                                          int offset);
}
