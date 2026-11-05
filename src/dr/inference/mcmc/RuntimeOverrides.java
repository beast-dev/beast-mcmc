/*
 * RuntimeOverrides.java
 *
 * Copyright © 2002-2026 the BEAST Development Team
 * http://beast.community/about
 *
 * This file is part of BEAST.
 * See the NOTICE file distributed with this work for additional
 * information regarding copyright ownership and licensing.
 *
 * BEAST is free software; you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as
 * published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.
 *
 *  BEAST is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with BEAST; if not, write to the
 * Free Software Foundation, Inc., 51 Franklin St, Fifth Floor,
 * Boston, MA  02110-1301  USA
 */

package dr.inference.mcmc;

import java.util.HashSet;
import java.util.Set;
import java.util.logging.Logger;

/**
 * The length of the chain and the frequency of logging, as given in the XML file, unless they have been
 * overridden for this run. The command line options -chain_length and -log_every set the properties read
 * here, so that a file can be run short, or sampled more or less often, without being edited.
 *
 * @author Andrew Rambaut
 */
public class RuntimeOverrides {

    public static final String CHAIN_LENGTH_PROPERTY = "mcmc.chain_length";
    public static final String LOG_EVERY_PROPERTY = "mcmc.log_every";

    /**
     * @param chainLength the length given in the file
     * @return the length to use for this run
     */
    public static long getChainLength(long chainLength) {
        Long override = getPositiveProperty(CHAIN_LENGTH_PROPERTY);

        if (override == null) {
            return chainLength;
        }

        if (!reportedChainLength) {
            Logger.getLogger("dr.inference").info("Overriding chain length: running for " + override +
                    " states rather than the " + chainLength + " given in the file");
            reportedChainLength = true;
        }
        return override;
    }

    /**
     * @param logEvery the frequency given in the file
     * @return the frequency to use for this run
     */
    public static int getLogEvery(int logEvery) {
        Long override = getPositiveProperty(LOG_EVERY_PROPERTY);

        // a logger that does not log at a fixed frequency is left as it is
        if (override == null || logEvery <= 0) {
            return logEvery;
        }

        if (!reportedLogEvery) {
            Logger.getLogger("dr.inference").info("Overriding logging frequency: logging every " + override +
                    " states");
            reportedLogEvery = true;
        }
        return (int) Math.min(override, Integer.MAX_VALUE);
    }

    /**
     * @return the value of the property, or null if it is not set; a value that is not a positive whole
     * number is reported and ignored, rather than stopping a run that would otherwise be fine
     */
    private static Long getPositiveProperty(String property) {
        String value = System.getProperty(property);

        if (value == null) {
            return null;
        }

        try {
            long parsed = Long.parseLong(value.trim());
            if (parsed > 0) {
                return parsed;
            }
            if (!reported.add(property)) {
                return null;
            }
            // a warning through the logger is only shown when running verbosely, and this should be seen
            System.err.println("WARNING: " + property + " must be greater than zero, but is " + value +
                    "; the value in the file will be used");
        } catch (NumberFormatException nfe) {
            if (reported.add(property)) {
                System.err.println("WARNING: " + property + " is not a whole number: " + value +
                        "; the value in the file will be used");
            }
        }
        return null;
    }

    /** the properties already complained about, so that each is reported once however many loggers there are */
    private static final Set<String> reported = new HashSet<String>();

    private static boolean reportedChainLength = false;
    private static boolean reportedLogEvery = false;
}
