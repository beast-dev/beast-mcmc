/*
 * RegressionJNIWrapper.java
 *
 * Copyright (c) 2002-2013 Alexei Drummond, Andrew Rambaut and Marc Suchard
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

package dr.inference.regression;

/**
 * @author Marc A. Suchard
 */
public class RegressionJNIWrapper {

    public static final String LIBRARY_NAME = "bsccs_jni";
    public static final String LIBRARY_PATH_LABEL = "bsccs.library.path";
    public static final String LIBRARY_PLATFORM_NAME = getPlatformSpecificLibraryName();

    public static final int NO_PRIOR = 0;
    public static final int LAPLACE_PRIOR = 1;
    public static final int NORMAL_PRIOR = 2;

    /**
     * private constructor to enforce singleton instance
     */
    private RegressionJNIWrapper() {
    }

    public native int loadData(String fileName);

    public native double getLogLikelihood(int instance);

    public native double getLogPrior(int instance);

    public native double getBeta(int instance, int index);

    public native int getBetaSize(int instance);

    public native void setBeta(int instance, int index, double value);

    public native void setBeta(int instance, double[] values);

    public native double getHyperprior(int instance);

    public native void setHyperprior(int instance, double value);

    public native void findMode(int instance);

    public native int getUpdateCount(int instance);

    public native int getLikelihoodCount(int instance);

    public native void setPriorType(int instance, int type);

    public native void makeDirty(int instance);

//    public native int createInstance(
//            int tipCount,
//            int partialsBufferCount,
//            int compactBufferCount,
//            int stateCount,
//            int patternCount,
//            int eigenBufferCount,
//            int matrixBufferCount,
//            int categoryCount,
//            int scaleBufferCount,
//          less t  final int[] resourceList,
//            int resourceCount,
//            long preferenceFlags,
//            long requirementFlags,
//            InstanceDetails returnInfo);
//
//    public native int finalize(int instance);
//
//    public native int setPatternWeights(int instance,
//                                        final double[] patternWeights);
//
//    public native int setTipStates(int instance, int tipIndex, final int[] inStates);
//
//    public native int getTipStates(int instance, int tipIndex, final int[] inStates);
//
//    public native int setTipPartials(int instance, int tipIndex, final double[] inPartials);
//
//    public native int setPartials(int instance, int bufferIndex, final double[] inPartials);
//
//    public native int getPartials(int instance, int bufferIndex, int scaleIndex,
//                                  final double[] outPartials);
//
//
//    public native int setEigenDecomposition(int instance,
//                                            int eigenIndex,
//                                            final double[] eigenVectors,
//                                            final double[] inverseEigenValues,
//                                            final double[] eigenValues);
//
//    public native int setStateFrequencies(int instance,
//                                          int stateFrequenciesIndex,
//                                          final double[] stateFrequencies);
//
//    public native int setCategoryWeights(int instance,
//                                         int categoryWeightsIndex,
//                                         final double[] categoryWeights);
//
//    public native int setCategoryRates(int instance,
//                                       final double[] inCategoryRates);
//
//    public native int setTransitionMatrix(int instance, int matrixIndex, final double[] inMatrix, double paddedValue);
//
//    public native int getTransitionMatrix(int instance, int matrixIndex, final double[] outMatrix);
//
//	public native int convolveTransitionMatrices(int instance,
//			                                     final int[] firstIndices,
//			                                     final int[] secondIndices,
//			                                     final int[] resultIndices,
//			                                     int matrixCount);
//
//    public native int updateTransitionMatrices(int instance, int eigenIndex,
//                                               final int[] probabilityIndices,
//                                               final int[] firstDerivativeIndices,
//                                               final int[] secondDervativeIndices,
//                                               final double[] edgeLengths,
//                                               int count);
//
//    public native int updatePartials(final int instance,
//                                     final int[] operations,
//                                     int operationCount,
//                                     int cumulativeScalingIndex);
//
//    public native int waitForPartials(final int instance,
//                                      final int[] destinationPartials,
//                                      int destinationPartialsCount);
//
//    public native int accumulateScaleFactors(final int instance,
//                                             final int[] scaleIndices,
//                                             final int count,
//                                             final int cumulativeScalingIndex);
//
//    public native int removeScaleFactors(final int instance,
//                                         final int[] scaleIndices,
//                                         final int count,
//                                         final int cumulativeScalingIndex);
//
//    public native int resetScaleFactors(final int instance,
//                                        final int cumulativeScalingIndex);
//
//    public native int copyScaleFactors(final int instance,
//                                       final int destScalingIndex,
//                                       final int srcScalingIndex);
//
//    public native int calculateRootLogLikelihoods(int instance,
//                                                  final int[] bufferIndices,
//                                                  final int[] categoryWeightsIndices,
//                                                  final int[] stateFrequenciesIndices,
//                                                  final int[] cumulativeScaleIndices,
//                                                  int count,
//                                                  final double[] outSumLogLikelihood);
//
//    public native int calculateEdgeLogLikelihoods(int instance,
//                                                  final int[] parentBufferIndices,
//                                                  final int[] childBufferIndices,
//                                                  final int[] probabilityIndices,
//                                                  final int[] firstDerivativeIndices,
//                                                  final int[] secondDerivativeIndices,
//                                                  final int[] categoryWeightsIndices,
//                                                  final int[] stateFrequenciesIndices,
//                                                  final int[] scalingFactorsIndices,
//                                                  int count,
//                                                  final double[] outSumLogLikelihood,
//                                                  final double[] outSumFirstDerivative,
//                                                  final double[] outSumSecondDerivative);
//
//    public native int getSiteLogLikelihoods(final int instance,
//                                            final double[] outLogLikelihoods);

    /* Library loading routines */

    private static String getPlatformSpecificLibraryName() {
        String osName = System.getProperty("os.name").toLowerCase();
        String osArch = System.getProperty("os.arch").toLowerCase();
        if (osName.startsWith("windows")) {
            if(osArch.equals("i386")) return LIBRARY_NAME + "32";
            if(osArch.startsWith("amd64")||osArch.startsWith("x86_64")) return LIBRARY_NAME + "64";
        }
        return LIBRARY_NAME;
    }

    public static RegressionJNIWrapper loadLibrary() throws UnsatisfiedLinkError {

        if (INSTANCE == null) {
            System.err.println("Trying to load BSCCS library...");
            String path = "";
            if (System.getProperty(LIBRARY_PATH_LABEL) != null) {
                path = System.getProperty(LIBRARY_PATH_LABEL);
                if (path.length() > 0 && !path.endsWith("/")) {
                    path += "/";
                }
            }

            System.loadLibrary(path + LIBRARY_PLATFORM_NAME);
            INSTANCE = new RegressionJNIWrapper();
            System.err.println("BSCCS library loaded.");
        }

        return INSTANCE;
    }

    private static RegressionJNIWrapper INSTANCE = null;
}
