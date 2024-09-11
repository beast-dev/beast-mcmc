/*
 * BastaJNIWrapper.java
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

public class BastaJNIWrapper {

    private static final String LIBRARY_NAME = getPlatformSpecificLibraryName();

    private BastaJNIWrapper() { }

    public native int allocateCoalescentBuffers(int instance,
                                                int bufferCount,
                                                int maxCoalescentIntervalCount); // TODO buffers have different sizes

    public native int getBastaBuffer(int instance,
                                     int index,
                                     double[] buffer);

    public native int updateBastaPartials(int instance,
                                          final int[] operations,
                                          int operationCount,
                                          final int[] intervals,
                                          int intervalCount,
                                          int populationSizeIndex,
                                          int coalescentProbabilityIndex);

    public native int accumulateBastaPartials(int instance,
                                              final int[] operations,
                                              int operationCount,
                                              final int[] intervals,
                                              int intervalCount,
                                              final double[] intervalLengths,
                                              int populationSizesIndex,
                                              int coalescentProbabilityIndex,
                                              double[] result);

    private static String getPlatformSpecificLibraryName() {
        String osName = System.getProperty("os.name").toLowerCase();
        String osArch = System.getProperty("os.arch").toLowerCase();
        if (osName.startsWith("windows")) {
            if (osArch.equals("x86") || osArch.equals("i386")) return "hmsbeagle-basta32";
            if (osArch.startsWith("amd64") || osArch.startsWith("x86_64")) return "hmsbeagle-basta64";
        }
        return "hmsbeagle-jni-basta";
    }

    public static void loadBastaLibrary() throws UnsatisfiedLinkError {
        String path = "";
        if (System.getProperty("beagle.library.path") != null) {
            path = System.getProperty("beagle.library.path");
            if (path.length() > 0 && !path.endsWith("/")) {
                path += "/";
            }
        }

        System.loadLibrary(path + LIBRARY_NAME);
        INSTANCE = new BastaJNIWrapper();
    }

    public static BastaJNIWrapper INSTANCE;

    public native int updateBastaPartialsGrad(int instance, int[] operations, int operationCount, int[] intervals, int intervalCount, int populationSizeIndex, int coalescentProbabilityIndex);

    public native int updateTransitionMatricesGrad(int instance, int[] transitionMatrixIndices, double[] branchLengths, int count);

    public native int accumulateBastaPartialsGrad(int instance, int[] operations, int operationCount, int[] intervalStarts, int intervalCount, double[] intervalLengths, int populationSizeIndex, int coalescentProbabilityIndex, double[] result);
}
