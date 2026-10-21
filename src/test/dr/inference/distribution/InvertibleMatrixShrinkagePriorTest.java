package test.dr.inference.distribution;

import dr.inference.distribution.InvertibleMatrixShrinkagePrior;
import dr.inference.model.GivensRotationMatrixParameter;
import dr.inference.model.LDUInvertibleMatrixParametrization;
import dr.inference.model.MatrixParameter;
import dr.inference.model.NormalizedMatrixParameter;
import dr.inference.model.Parameter;
import dr.inference.model.QRUnitUpperMatrixParametrization;
import test.dr.math.MathTestCase;

public class InvertibleMatrixShrinkagePriorTest extends MathTestCase {

    public void testLogLikelihoodMatchesClosedFormForLdu() {
        final LDUInvertibleMatrixParametrization matrix =
                new LDUInvertibleMatrixParametrization(
                        "ldu",
                        3,
                        new Parameter.Default("lower", new double[]{0.2, -0.3, 0.4}),
                        new Parameter.Default("logDiagonal", new double[]{
                                Math.log(2.0), Math.log(0.5), Math.log(1.5)}),
                        new Parameter.Default("upper", new double[]{-0.1, 0.6, 0.25}));
        final Parameter lambda = new Parameter.Default("lambda", new double[]{0.7});
        final InvertibleMatrixShrinkagePrior prior =
                new InvertibleMatrixShrinkagePrior("prior", matrix, lambda);

        double frobeniusSquared = 0.0;
        for (int row = 0; row < 3; row++) {
            for (int col = 0; col < 3; col++) {
                final double value = matrix.getParameterValue(row, col);
                frobeniusSquared += value * value;
            }
        }
        final double expected = lambda.getParameterValue(0) * matrix.getLogAbsDeterminant()
                - 0.5 * lambda.getParameterValue(0) * frobeniusSquared;

        assertEquals(expected, prior.getLogLikelihood(), 1e-12);
    }

    public void testDimensionScaledLogLikelihoodMatchesClosedFormForLdu() {
        final LDUInvertibleMatrixParametrization matrix =
                new LDUInvertibleMatrixParametrization(
                        "ldu",
                        3,
                        new Parameter.Default("lower", new double[]{0.2, -0.3, 0.4}),
                        new Parameter.Default("logDiagonal", new double[]{
                                Math.log(2.0), Math.log(0.5), Math.log(1.5)}),
                        new Parameter.Default("upper", new double[]{-0.1, 0.6, 0.25}));
        final Parameter lambda = new Parameter.Default("lambda", new double[]{0.7});
        final InvertibleMatrixShrinkagePrior prior =
                new InvertibleMatrixShrinkagePrior("prior", matrix, lambda, true);

        double frobeniusSquared = 0.0;
        for (int row = 0; row < 3; row++) {
            for (int col = 0; col < 3; col++) {
                final double value = matrix.getParameterValue(row, col);
                frobeniusSquared += value * value;
            }
        }
        final double expected = lambda.getParameterValue(0) * matrix.getLogAbsDeterminant()
                - 0.5 * lambda.getParameterValue(0) * 3.0 * frobeniusSquared;

        assertEquals(expected, prior.getLogLikelihood(), 1e-12);
    }

    public void testGradientMatchesFiniteDifferenceForLdu() {
        final LDUInvertibleMatrixParametrization matrix =
                new LDUInvertibleMatrixParametrization(
                        "ldu",
                        3,
                        new Parameter.Default("lower", new double[]{0.2, -0.3, 0.4}),
                        new Parameter.Default("logDiagonal", new double[]{
                                Math.log(2.0), Math.log(0.5), Math.log(1.5)}),
                        new Parameter.Default("upper", new double[]{-0.1, 0.6, 0.25}));
        final Parameter lambda = new Parameter.Default("lambda", new double[]{0.7});
        final InvertibleMatrixShrinkagePrior prior =
                new InvertibleMatrixShrinkagePrior("prior", matrix, lambda);

        checkGradientAgainstFiniteDifference(prior);
    }

    public void testDimensionScaledGradientMatchesFiniteDifferenceForLdu() {
        final LDUInvertibleMatrixParametrization matrix =
                new LDUInvertibleMatrixParametrization(
                        "ldu",
                        3,
                        new Parameter.Default("lower", new double[]{0.2, -0.3, 0.4}),
                        new Parameter.Default("logDiagonal", new double[]{
                                Math.log(2.0), Math.log(0.5), Math.log(1.5)}),
                        new Parameter.Default("upper", new double[]{-0.1, 0.6, 0.25}));
        final Parameter lambda = new Parameter.Default("lambda", new double[]{0.7});
        final InvertibleMatrixShrinkagePrior prior =
                new InvertibleMatrixShrinkagePrior("prior", matrix, lambda, true);

        checkGradientAgainstFiniteDifference(prior);
    }

    public void testGradientMatchesFiniteDifferenceForQr() {
        final Parameter angles = new Parameter.Default("angles", new double[]{0.2, -0.4, 0.7});
        final Parameter upper = new Parameter.Default("upper", new double[]{-0.1, 0.6, 0.25});
        final QRUnitUpperMatrixParametrization matrix =
                new QRUnitUpperMatrixParametrization(
                        "qr",
                        new GivensRotationMatrixParameter("orthogonal", 3, angles),
                        upper);
        final Parameter lambda = new Parameter.Default("lambda", new double[]{1.3});
        final InvertibleMatrixShrinkagePrior prior =
                new InvertibleMatrixShrinkagePrior("prior", matrix, lambda);

        checkGradientAgainstFiniteDifference(prior);
    }

    public void testGradientMatchesFiniteDifferenceForPlainDenseMatrix() {
        final MatrixParameter matrix = denseMatrix("dense");
        final Parameter lambda = new Parameter.Default("lambda", new double[]{0.9});
        final InvertibleMatrixShrinkagePrior prior =
                new InvertibleMatrixShrinkagePrior("prior", matrix, lambda);

        assertSame(matrix, prior.getParameter());
        checkGradientAgainstFiniteDifference(prior);
    }

    public void testGradientMatchesFiniteDifferenceForGlobalNormalizedDenseMatrix() {
        final MatrixParameter raw = denseMatrix("raw");
        final NormalizedMatrixParameter matrix =
                new NormalizedMatrixParameter("normalized", raw, NormalizedMatrixParameter.Normalization.GLOBAL);
        final Parameter lambda = new Parameter.Default("lambda", new double[]{1.1});
        final InvertibleMatrixShrinkagePrior prior =
                new InvertibleMatrixShrinkagePrior("prior", matrix, lambda);

        assertSame(raw, prior.getParameter());
        checkGradientAgainstFiniteDifference(prior);
    }

    private static MatrixParameter denseMatrix(final String id) {
        return new MatrixParameter(id, new Parameter[]{
                new Parameter.Default(id + ".col0", new double[]{2.0, 0.1, -0.2}),
                new Parameter.Default(id + ".col1", new double[]{0.3, 1.5, 0.4}),
                new Parameter.Default(id + ".col2", new double[]{-0.4, 0.2, 1.2})
        });
    }

    private static void checkGradientAgainstFiniteDifference(final InvertibleMatrixShrinkagePrior prior) {
        final Parameter nativeParameter = prior.getParameter();
        final double[] analytic = prior.getGradientLogDensity();
        assertEquals(nativeParameter.getDimension(), analytic.length);

        final double h = 1e-6;
        for (int i = 0; i < nativeParameter.getDimension(); i++) {
            final double original = nativeParameter.getParameterValue(i);

            nativeParameter.setParameterValue(i, original + h);
            final double plus = prior.getLogLikelihood();

            nativeParameter.setParameterValue(i, original - h);
            final double minus = prior.getLogLikelihood();

            nativeParameter.setParameterValue(i, original);

            final double finiteDifference = (plus - minus) / (2.0 * h);
            assertEquals("gradient mismatch at native index " + i, finiteDifference, analytic[i], 1e-5);
        }
    }
}
