/*
 * SpectraJNIWrapper.java
 *
 * Copyright (c) 2002-2017 Alexei Drummond, Andrew Rambaut and Marc Suchard
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


package dr.evomodel.substmodel.spectra;

/*
 * SpectraJNIWrapper
 *
 * @author Xiang Ji
 *
 */
public class SpectraJNIWrapper {

    private static SpectraJNIWrapper INSTANCE = null;

//    static {
//        System.load("/usr/local/lib/libspectra-jni.jnilib");
//    }

    private static final String LIBRARY_NAME = "spectra-jni";
    private static final String LIBRARY_PATH_LABEL = "java.library.path";
    private static final String LIBRARY_PLATFORM_NAME = getPlatformSpecificLibraryName();
    private static final String LIBRARY_PLATFORM_EXTENSION = getPlatformSpecificLibraryExtension();
    private static final String LIBRARY_PLATFORM_PREFIX = getPlatformSpecificLibraryPrefix();

    public SpectraJNIWrapper() {
    }

    public native int createInstance(int matrixCount,
                                     int stateCount);

    public native int setMatrix(int matrix,
                                int[] indices,
                                double[] values,
                                int nonZeroCount);

    public native int getEigenVectors(int matrix,
                                      int numEigenValues,
                                      double sigma,
                                      double[] eigenValues,
                                      double[] eigenVectors);

    public native String getVersion();

    public static SpectraJNIWrapper loadLibrary() throws UnsatisfiedLinkError {

        if (INSTANCE == null) {
            System.err.println("Trying to load Spectra library...");
            String path = "";
            if (System.getProperty(LIBRARY_PATH_LABEL) != null) {
                path = System.getProperty(LIBRARY_PATH_LABEL);
                if (path.length() > 0 && !path.endsWith("/")) {
                    path += "/";
                }
                System.load(path + LIBRARY_PLATFORM_PREFIX + LIBRARY_NAME + LIBRARY_PLATFORM_EXTENSION);
            } else {
                System.loadLibrary(LIBRARY_PLATFORM_NAME);
            }

            INSTANCE = new SpectraJNIWrapper();
            System.err.println("Spectra library loaded.");
        }

        return INSTANCE;
    }

    private static String getPlatformSpecificLibraryName() {
        String osName = System.getProperty("os.name").toLowerCase();
        String osArch = System.getProperty("os.arch").toLowerCase();
        if (osName.startsWith("windows")) {
            if (osArch.equals("x86") || osArch.equals("i386")) return LIBRARY_NAME + "32";
            if (osArch.startsWith("amd64") || osArch.startsWith("x86_64")) return LIBRARY_NAME + "64";
        }
        return LIBRARY_NAME;
    }
    private static String getPlatformSpecificLibraryExtension() {
        String osName = System.getProperty("os.name").toLowerCase();
        if (osName.startsWith("windows")) {
            return ".dll";
        } else if (osName.startsWith("mac")) {
            return ".jnilib";
        } else {
            return ".so";
        }
    }
    private static String getPlatformSpecificLibraryPrefix() {
        String osName = System.getProperty("os.name").toLowerCase();
        if (osName.startsWith("windows")) {
            return "";
        } else {
            return "lib";
        }
    }

}