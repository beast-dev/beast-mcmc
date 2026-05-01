package dr.evomodel.treedatalikelihood.continuous.canonical;

import java.util.logging.Logger;

/**
 * Logging entry point for canonical OU diagnostics.
 */
public final class CanonicalDiagnosticsLog {

    private static final Logger LOGGER = Logger.getLogger("dr.evomodel.treedatalikelihood.continuous.canonical");

    private CanonicalDiagnosticsLog() { }

    public static void warning(final String message) {
        LOGGER.warning(message);
    }

    public static void info(final String message) {
        LOGGER.info(message);
    }
}
