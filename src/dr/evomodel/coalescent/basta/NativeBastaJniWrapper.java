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

package dr.evomodel.coalescent.basta;

import dr.evomodel.treedatalikelihood.continuous.cdi.ResourceDetails;

/**
 * @author Marc A. Suchard
 */
public class NativeBastaJniWrapper {

    /**
     * private constructor to enforce singleton instance
     */
    private NativeBastaJniWrapper() { }

    public native void functionOne();

    public native void functionTwo();

    public native String getVersion();

    public native String getCitation();

    public native ResourceDetails[] getResourceList();

    public static NativeBastaJniWrapper getBastaJniWrapper() {
        if (NativeBastaJniWrapper.INSTANCE == null) {
            try {
                NativeBastaJniWrapper.loadBastaLibrary();

            } catch (UnsatisfiedLinkError ule) {
                System.err.println("Failed to load native BASTA library: " + ule.getMessage());
            }
        }

        return NativeBastaJniWrapper.INSTANCE;
    }

    /* Library loading routines */

    private static String getPlatformSpecificLibraryName() {
        String osName = System.getProperty("os.name").toLowerCase();
        String osArch = System.getProperty("os.arch").toLowerCase();
        if (osName.startsWith("windows")) {
            if(osArch.equals("x86")||osArch.equals("i386")) return "native_basta32";
            if(osArch.startsWith("amd64")||osArch.startsWith("x86_64")) return "native_basta64";
        }
        return "native_basta-jni";
    }

    private static void loadBastaLibrary() throws UnsatisfiedLinkError {
        String path = "";
        if (System.getProperty("basta.library.path") != null) {
            path = System.getProperty("basta.library.path");
            if (path.length() > 0 && !path.endsWith("/")) {
                path += "/";
            }
        }

        System.loadLibrary(path + LIBRARY_NAME);
        INSTANCE = new NativeBastaJniWrapper();
    }

    private static final String LIBRARY_NAME = getPlatformSpecificLibraryName();
    private static NativeBastaJniWrapper INSTANCE;
}
