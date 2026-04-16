/*
 * EigenJNIWrapper.java
 *
 * Copyright (c) 2002-2024 Alexei Drummond, Andrew Rambaut and Marc Suchard
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
package dr.evomodel.substmodel.eigen;

/*
 * EigenJNIWrapper
 *
 * @author Xiang Ji
 *
 */

public class EigenJNIWrapper implements AutoCloseable {
    private long nativeHandle;

    static {
        System.loadLibrary("eigen-jni");
    }

    public EigenJNIWrapper() {
    }

    public void createInstance(int matrixCount, int stateCount) {
        nativeHandle = nativeCreateInstance(matrixCount, stateCount);
        if (nativeHandle == 0) {
            throw new RuntimeException("Failed to create native EigenImpl");
        }
    }

    public int setMatrix(int matrix, int[] indices, double[] values,
                         int nonZeroCount, int stateCount) {
        return nativeSetMatrix(nativeHandle, matrix, indices, values,
                               nonZeroCount, stateCount);
    }

    public int getEigenDecomposition(int matrix, double[] eigenValues,
                                     double[] eigenVectors,
                                     double[] inverseEigenVectors) {
        return nativeGetEigenDecomposition(nativeHandle, matrix,
                                           eigenValues, eigenVectors,
                                           inverseEigenVectors);
    }

    @Override
    public void close() {
        if (nativeHandle != 0) {
            nativeDestroyInstance(nativeHandle);
            nativeHandle = 0;
        }
    }

    public native String getVersion();
    private static native long nativeCreateInstance(int matrixCount, int stateCount);
    private static native void nativeDestroyInstance(long handle);
    private static native int nativeSetMatrix(long handle, int matrix,
                                              int[] indices, double[] values,
                                              int nonZeroCount, int stateCount);
    private static native int nativeGetEigenDecomposition(long handle, int matrix,
                                                          double[] eigenValues,
                                                          double[] eigenVectors,
                                                          double[] inverseEigenVectors);
}
