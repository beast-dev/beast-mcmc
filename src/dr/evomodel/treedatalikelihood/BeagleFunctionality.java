/*
 * BeagleFunctionality.java
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

package dr.evomodel.treedatalikelihood;

import beagle.BeagleInfo;

import java.util.ArrayList;
import java.util.List;

/**
 * @author Marc A. Suchard
 */
public class BeagleFunctionality {

    private static boolean checkGTEVersion(int[] versionNumbers){
        int[] beagleVersionNumbers = BeagleInfo.getVersionNumbers();
        if (versionNumbers.length == 0 || beagleVersionNumbers.length == 0)
            return false;
        for (int i = 0; i < versionNumbers.length && i < beagleVersionNumbers.length; i++){
            if (beagleVersionNumbers[i] > versionNumbers[i])
                return true;
            if (beagleVersionNumbers[i] < versionNumbers[i])
                return false;
        }
        return true;
    }

    public static boolean IS_THREAD_COUNT_COMPATIBLE() {
        return checkGTEVersion(new int[]{3,1});
    }

    public static boolean IS_ODD_STATE_SSE_FIXED() {
        // SSE for odd state counts fixed in BEAGLE 3.1.3
        return checkGTEVersion(new int[]{3,1,3});
    }

    public static boolean IS_MATRIX_VECTOR_SUPPORTED() {
        // the spectral representation, top pre-order partials (updatePrePartials_v5), adjoint cross-products, degree-2
        // operations (a second child of Beagle.NONE) and buffer growth (ensureBufferCounts) are in native BEAGLE from
        // 4.2.0; a 4.1.x version may or may not have them (development builds of matrix_vector were 4.1.0 to 4.1.2),
        // and the Java implementation (-java) has none of them
        return !Boolean.getBoolean("java.only") && checkGTEVersion(new int[]{4,2,0});
    }

    public static void checkMatrixVectorSupported(String feature) {
        if (!IS_MATRIX_VECTOR_SUPPORTED()) {
            throw new IllegalArgumentException(feature + " needs native BEAGLE 4.2.0 or later" +
                    (Boolean.getBoolean("java.only") ? ", not the Java implementation (-java)" :
                            "; found " + BeagleInfo.getVersion()));
        }
    }

    public static boolean IS_BASTA_GRADIENT_SUPPORTED() {
        // BASTA gradients through native BEAGLE-BIT (allocateCoalescentGradBuffers, the 13-argument
        // accumulateBastaPartialsGrad, getBastaSlabConstants) are in BEAGLE from 4.2.0, while the BASTA likelihood
        // already works on 4.1.0; version only, not java.only, because BastaFactory does not honour -java
        return checkGTEVersion(new int[]{4,2,0});
    }

    static boolean IS_PRE_ORDER_SUPPORTED() {
        int[] versionNumbers = BeagleInfo.getVersionNumbers();
        return checkGTEVersion(new int[]{3,2});
    }

    static boolean IS_MULTI_PARTITION_COMPATIBLE() {
        int[] versionNumbers = BeagleInfo.getVersionNumbers();
        return checkGTEVersion(new int[]{3});
    }

    public static List<Integer> parseSystemPropertyIntegerArray(String propertyName) {
        List<Integer> order = new ArrayList<>();
        String r = System.getProperty(propertyName);
        if (r != null) {
            String[] parts = r.split(",");
            for (String part : parts) {
                try {
                    int n = Integer.parseInt(part.trim());
                    order.add(n);
                } catch (NumberFormatException nfe) {
                    System.err.println("Invalid entry '" + part + "' in " + propertyName);
                }
            }
        }
        return order;
    }

    public static List<Long> parseSystemPropertyLongArray(String propertyName) {
        List<Long> order = new ArrayList<>();
        String r = System.getProperty(propertyName);
        if (r != null) {
            String[] parts = r.split(",");
            for (String part : parts) {
                try {
                    long n = Long.parseLong(part.trim());
                    order.add(n);
                } catch (NumberFormatException nfe) {
                    System.err.println("Invalid entry '" + part + "' in " + propertyName);
                }
            }
        }
        return order;
    }

    public static List<String> parseSystemPropertyStringArray(String propertyName) {

        List<String> order = new ArrayList<>();

        String r = System.getProperty(propertyName);
        if (r != null) {
            String[] parts = r.split(",");
            for (String part : parts) {
                try {
                    String s = part.trim();
                    order.add(s);
                } catch (NumberFormatException nfe) {
                    System.err.println("Invalid entry '" + part + "' in " + propertyName);
                }
            }
        }
        return order;
    }
}
