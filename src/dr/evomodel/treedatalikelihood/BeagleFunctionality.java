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
