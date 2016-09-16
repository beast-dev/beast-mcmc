/*
 * CDIJNIWrapper.java
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

/**
 * Created by msuchard on 9/14/16.
 */
public class CDIJNIWrapper {

    public static final String LIBRARY_NAME = getPlatformSpecificLibraryName();

    /**
     * private constructor to enforce singleton instance
     */
    private CDIJNIWrapper() { }

    public native String getVersion();

    public native String getCitation();

    public native ResourceDetails[] getResourceList();

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
//            final int[] resourceList,
//            int resourceCount,
//            long preferenceFlags,
//            long requirementFlags,
//            InstanceDetails returnInfo);
//
//    public native int finalize(int instance);
//

    /* Library loading routines */

    private static String getPlatformSpecificLibraryName() {
        String osName = System.getProperty("os.name").toLowerCase();
        String osArch = System.getProperty("os.arch").toLowerCase();
        if (osName.startsWith("windows")) {
            if(osArch.equals("x86")||osArch.equals("i386")) return "cdi32";
            if(osArch.startsWith("amd64")||osArch.startsWith("x86_64")) return "cdi64";
        }
        return "cdi-jni";
    }

    public static void loadCDILibrary() throws UnsatisfiedLinkError {
        String path = "";
        if (System.getProperty("cdi.library.path") != null) {
            path = System.getProperty("cdi.library.path");
            if (path.length() > 0 && !path.endsWith("/")) {
                path += "/";
            }
        }

        System.loadLibrary(path + LIBRARY_NAME);
        INSTANCE = new CDIJNIWrapper();
    }

    public static CDIJNIWrapper INSTANCE;
}
