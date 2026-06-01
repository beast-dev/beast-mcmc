package test.dr.inference.timeseries;

import dr.evomodel.continuous.ou.OUProcessModel;
import dr.evomodel.treedatalikelihood.continuous.canonical.math.MatrixOps;
import dr.inference.timeseries.core.LatentProcessModel;
import dr.inference.timeseries.core.TimeGrid;
import dr.inference.timeseries.model.gaussian.EulerOUProcessModel;
import dr.inference.timeseries.model.gaussian.OUTimeSeriesProcessAdapter;
import dr.inference.timeseries.representation.GaussianTransitionRepresentation;
import dr.inference.timeseries.representation.RepresentableProcess;

public final class OUTimeSeriesTestSupport {

    private OUTimeSeriesTestSupport() {
        // no instances
    }

    public static RepresentableProcess representable(final OUProcessModel process) {
        return new OUTimeSeriesProcessAdapter(process);
    }

    public static LatentProcessModel latent(final OUProcessModel process) {
        return new OUTimeSeriesProcessAdapter(process);
    }

    public static RepresentableProcess representable(final EulerOUProcessModel process) {
        return process;
    }

    public static <T> T representation(final OUProcessModel process, final Class<T> representationClass) {
        return new OUTimeSeriesProcessAdapter(process).getRepresentation(representationClass);
    }

    public static <T> T representation(final EulerOUProcessModel process, final Class<T> representationClass) {
        return process.getRepresentation(representationClass);
    }

    public static <T> T representation(final RepresentableProcess process, final Class<T> representationClass) {
        return process.getRepresentation(representationClass);
    }

    public static boolean supportsRepresentation(final OUProcessModel process, final Class<?> representationClass) {
        return new OUTimeSeriesProcessAdapter(process).supportsRepresentation(representationClass);
    }

    public static GaussianTransitionRepresentation transitionRepresentation(final OUProcessModel process) {
        return representation(process, GaussianTransitionRepresentation.class);
    }

    public static void getTransitionMatrix(final OUProcessModel process,
                                    final int fromIndex,
                                    final int toIndex,
                                    final TimeGrid timeGrid,
                                    final double[][] out) {
        final int dim = out.length;
        final double[] flat = new double[dim * dim];
        transitionRepresentation(process).getTransitionMatrixFlat(fromIndex, toIndex, timeGrid, flat);
        MatrixOps.fromFlat(flat, out, dim);
    }

    public static void getTransitionMatrix(final EulerOUProcessModel process,
                                           final int fromIndex,
                                           final int toIndex,
                                           final TimeGrid timeGrid,
                                           final double[][] out) {
        final int dim = out.length;
        final double[] flat = new double[dim * dim];
        process.getTransitionMatrixFlat(fromIndex, toIndex, timeGrid, flat);
        MatrixOps.fromFlat(flat, out, dim);
    }

    public static void getTransitionOffset(final OUProcessModel process,
                                    final int fromIndex,
                                    final int toIndex,
                                    final TimeGrid timeGrid,
                                    final double[] out) {
        transitionRepresentation(process).getTransitionOffset(fromIndex, toIndex, timeGrid, out);
    }

    public static void getTransitionOffset(final EulerOUProcessModel process,
                                           final int fromIndex,
                                           final int toIndex,
                                           final TimeGrid timeGrid,
                                           final double[] out) {
        process.getTransitionOffset(fromIndex, toIndex, timeGrid, out);
    }

    public static void getTransitionCovariance(final OUProcessModel process,
                                        final int fromIndex,
                                        final int toIndex,
                                        final TimeGrid timeGrid,
                                        final double[][] out) {
        final int dim = out.length;
        final double[] flat = new double[dim * dim];
        transitionRepresentation(process).getTransitionCovarianceFlat(fromIndex, toIndex, timeGrid, flat);
        MatrixOps.fromFlat(flat, out, dim);
    }

    public static void getTransitionCovariance(final EulerOUProcessModel process,
                                               final int fromIndex,
                                               final int toIndex,
                                               final TimeGrid timeGrid,
                                               final double[][] out) {
        final int dim = out.length;
        final double[] flat = new double[dim * dim];
        process.getTransitionCovarianceFlat(fromIndex, toIndex, timeGrid, flat);
        MatrixOps.fromFlat(flat, out, dim);
    }
}
