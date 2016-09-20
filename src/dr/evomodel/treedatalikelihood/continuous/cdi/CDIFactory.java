/*
 * CDIFactory.java
 *
 * Copyright (c) 2002-2016 Alexei Drummond, Andrew Rambaut and Marc Suchard
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

package dr.evomodel.treedatalikelihood.continuous.cdi;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.logging.Logger;

/**
 * @author Marc A. Suchard
 */
public class CDIFactory {

    private static Map<Integer, ResourceDetails> resourceDetailsMap = new HashMap<Integer, ResourceDetails>();

    public static String getVersionInformation() {
        getCDIJNIWrapper();

        if (CDIJNIWrapper.INSTANCE != null) {

            try {
                return CDIJNIWrapper.INSTANCE.getCitation();
            } catch (UnsatisfiedLinkError ule) {
                return "Using CDI library v0.1 for accelerated, parallel likelihood evaluation";

            }
        }

        return "CDI not installed/found";
    }

    public static String getVersion() {
        getCDIJNIWrapper();

        if (CDIJNIWrapper.INSTANCE != null) {

            try {
                return CDIJNIWrapper.INSTANCE.getVersion();
            } catch (UnsatisfiedLinkError ule) {
                return "0.1";
            }
        }

        return "CDI not installed/found";
    }

    public static List<ResourceDetails> getResourceDetails() {
        getCDIJNIWrapper();

        return new ArrayList<ResourceDetails>(resourceDetailsMap.values());
    }

    public static ResourceDetails getResourceDetails(int resourceNumber) {
        getCDIJNIWrapper();
        return resourceDetailsMap.get(resourceNumber);
    }

    public static ContinuousDiffusionIntegrator loadCDIInstance(
            final PrecisionType precisionType,
            final int numTraits,
            final int dimTrait,
            final int bufferCount,
            final int diffusionCount
    ) {
        boolean forceJava = Boolean.valueOf(System.getProperty("java.only"));

        if (!forceJava) {
            getCDIJNIWrapper();
            if (CDIJNIWrapper.INSTANCE != null) {
                try {
                    // Do nothing yet
                    ContinuousDiffusionIntegrator cdi = null;

//                Beagle beagle = new BeagleJNIImpl(
//                        tipCount,
//                );

                    InstanceDetails details = cdi.getDetails();

                    if (details != null) // If resourceList/requirements not met, details == null here
                        return cdi;

                } catch (CDIException cdiException) {
                    Logger.getLogger("cdi").info("  " + cdiException.getMessage());
                }

            }
        }

        return new ContinuousDiffusionIntegrator.Basic(
                precisionType, numTraits, dimTrait, bufferCount, diffusionCount
        );
    }

    private static CDIJNIWrapper getCDIJNIWrapper() {
        if (CDIJNIWrapper.INSTANCE == null) {
            try {
                CDIJNIWrapper.loadCDILibrary();

            } catch (UnsatisfiedLinkError ule) {
                System.err.println("Failed to load CDI library: " + ule.getMessage());
            }

            if (CDIJNIWrapper.INSTANCE != null) {
                for (ResourceDetails details : CDIJNIWrapper.INSTANCE.getResourceList()) {
                    resourceDetailsMap.put(details.getNumber(), details);
                }
            }

        }

        return CDIJNIWrapper.INSTANCE;
    }
}
