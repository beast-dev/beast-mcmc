/*
 * NativeMDSSingleton.java
 *
 * Copyright (c) 2002-2014 Alexei Drummond, Andrew Rambaut and Marc Suchard
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

package dr.app.beagle.multidimensionalscaling;

/**
 * NativeMDSSingleton
 *
 * @author Marc Suchard
 * @author Andrew Rambaut
 * @version $Id$
 *          <p/>
 *          $HeadURL$
 *          <p/>
 *          $LastChangedBy$
 *          $LastChangedDate$
 *          $LastChangedRevision$
 */
public class NativeMDSSingleton {

    public static final String LIBRARY_NAME = "mds_jni";
    public static final String LIBRARY_PATH_LABEL = "mds.library.path";
    public static final String LIBRARY_PLATFORM_NAME = getPlatformSpecificLibraryName();


    private NativeMDSSingleton() {
    } // ensure singleton

    private static String getPlatformSpecificLibraryName() {
        String osName = System.getProperty("os.name").toLowerCase();
        String osArch = System.getProperty("os.arch").toLowerCase();
        if (osName.startsWith("windows")) {
            if (osArch.equals("i386")) return LIBRARY_NAME + "32";
            if (osArch.startsWith("amd64") || osArch.startsWith("x86_64")) return LIBRARY_NAME + "64";
        }
        return LIBRARY_NAME;
    }

    public static NativeMDSSingleton loadLibrary() throws UnsatisfiedLinkError {

        if (INSTANCE == null) {
            System.err.println("Trying to load MDS library...");
            String path = "";
            if (System.getProperty(LIBRARY_PATH_LABEL) != null) {
                path = System.getProperty(LIBRARY_PATH_LABEL);
                if (path.length() > 0 && !path.endsWith("/")) {
                    path += "/";
                }
            }

            System.loadLibrary(path + LIBRARY_PLATFORM_NAME);
            INSTANCE = new NativeMDSSingleton();
            System.err.println("MDS library loaded.");
        }

        return INSTANCE;
    }

    private static NativeMDSSingleton INSTANCE = null;

    public native int initialize(int dimensionCount, int locationCount, long flags);

    public native void updateLocations(int instance, int updateCount, double[] locations);

    public native double calculateLogLikelihood(int instance);

    public native void storeState(int instance);

    public native void restoreState(int instance);

    public native void makeDirty(int instance);

    public native void setPairwiseData(int instance, double[] observations);

    public native void setParameters(int instance, double[] parameters);

}
