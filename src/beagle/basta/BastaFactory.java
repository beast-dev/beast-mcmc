/*
 * BastaFactory.java
 *
 * Copyright © 2002-2024 the BEAST Development Team
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
 *
 */

package beagle.basta;

import beagle.*;

import java.util.logging.Logger;

public class BastaFactory extends BeagleFactory {

    public static BeagleBasta loadBastaInstance(
            int tipCount,
            int coalescentBufferCount,
            int maxCoalescentIntervalCount,
            int partialsBufferCount,
            int compactBufferCount,
            int stateCount,
            int patternCount,
            int eigenBufferCount,
            int matrixBufferCount,
            int categoryCount,
            int scaleBufferCount,
            int[] resourceList,
            long preferenceFlags,
            long requirementFlags) {

        BeagleJNIWrapper beagleWrapper = getBeagleJNIWrapper();
        if (beagleWrapper != null) {

            BastaJNIWrapper bastaWrapper = getBastaJNIWrapper();
            if (bastaWrapper != null) {

                try {
                    BeagleBasta beagle = new BastaJNIImpl(
                            tipCount,
                            coalescentBufferCount,
                            maxCoalescentIntervalCount,
                            partialsBufferCount,
                            compactBufferCount,
                            stateCount,
                            patternCount,
                            eigenBufferCount,
                            matrixBufferCount,
                            categoryCount,
                            scaleBufferCount,
                            resourceList,
                            preferenceFlags,
                            requirementFlags
                    );

                    // In order to know that it was a CPU instance created, we have to let BEAGLE
                    // to make the instance and then override it...

                    InstanceDetails details = beagle.getDetails();

                    if (details != null) // If resourceList/requirements not met, details == null here
                        return beagle;

                } catch (BeagleException beagleException) {
                    Logger.getLogger("beagle").info("  " + beagleException.getMessage());
                }
            } else {
                throw new RuntimeException("No acceptable BEAGLE-BASTA library plugin found. " +
                        "Make sure that BEAGLE-BASTA is properly installed or try changing resource requirements.");
            }
        }

        throw new RuntimeException("No acceptable BEAGLE library plugins found. " +
                "Make sure that BEAGLE is properly installed or try changing resource requirements.");
    }

    private static BastaJNIWrapper getBastaJNIWrapper() {
        if (BastaJNIWrapper.INSTANCE == null) {
            try {
                BastaJNIWrapper.loadBastaLibrary();
            } catch (UnsatisfiedLinkError ule) {
                System.err.println("Failed to load BEAGLE-BASTA library: " + ule.getMessage());
            }
        }

        return BastaJNIWrapper.INSTANCE;
    }
}
