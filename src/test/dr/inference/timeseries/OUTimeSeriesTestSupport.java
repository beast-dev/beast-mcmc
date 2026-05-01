package test.dr.inference.timeseries;

import dr.evomodel.continuous.ou.OUProcessModel;
import dr.inference.timeseries.core.LatentProcessModel;
import dr.inference.timeseries.core.TimeGrid;
import dr.inference.timeseries.gaussian.EulerOUProcessModel;
import dr.inference.timeseries.gaussian.OUTimeSeriesProcessAdapter;
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
        transitionRepresentation(process).getTransitionMatrix(fromIndex, toIndex, timeGrid, out);
    }

    public static void getTransitionMatrix(final EulerOUProcessModel process,
                                           final int fromIndex,
                                           final int toIndex,
                                           final TimeGrid timeGrid,
                                           final double[][] out) {
        process.getTransitionMatrix(fromIndex, toIndex, timeGrid, out);
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
        transitionRepresentation(process).getTransitionCovariance(fromIndex, toIndex, timeGrid, out);
    }

    public static void getTransitionCovariance(final EulerOUProcessModel process,
                                               final int fromIndex,
                                               final int toIndex,
                                               final TimeGrid timeGrid,
                                               final double[][] out) {
        process.getTransitionCovariance(fromIndex, toIndex, timeGrid, out);
    }
}
