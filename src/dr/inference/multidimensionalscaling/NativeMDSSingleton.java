/*
 * NativeMDSSingleton.java
 *
 * Copyright (c) 2002-2015 Alexei Drummond, Andrew Rambaut and Marc Suchard
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

package dr.inference.multidimensionalscaling;

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

    private static final String LIBRARY_NAME = "mds_jni";
    private static final String LIBRARY_PATH_LABEL = "mds.library.path";
    private static final String LIBRARY_PLATFORM_NAME = getPlatformSpecificLibraryName();
    private static final String LIBRARY_PLATFORM_EXTENSION = getPlatformSpecificLibraryExtension();
    private static final String LIBRARY_PLATFORM_PREFIX = getPlatformSpecificLibraryPrefix();

    static final String THREADS = "mds.threads";
    static final String MDS_RESOURCE = "mds.resource";

//    private int getThreads() {
//        String r = System.getProperty(THREADS);
//        int i = 1;
//        if (r != null) {
//            i = Integer.parseInt(r.trim());
//        }
//        return i;
//    }

    private NativeMDSSingleton() {
    } // ensure singleton

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
            return ".dylib";
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

    public static NativeMDSSingleton loadLibrary() throws UnsatisfiedLinkError {

        if (INSTANCE == null) {
            System.err.println("Trying to load MDS library...");
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

            INSTANCE = new NativeMDSSingleton();
            System.err.println("MDS library loaded.");
        }

        return INSTANCE;
    }

    private static NativeMDSSingleton INSTANCE = null;

    public int initialize(int dimensionCount, int locationCount, MultiDimensionalScalingCore.CoreInformation information) {
        return initialize(dimensionCount, locationCount,
                information.flags, information.deviceNumber, information.numThreads);
    }

    private native int initialize(int dimensionCount, int locationCount, long flags, int deviceNumber, int threads);

    public native void updateLocations(int instance, int updateCount, double[] locations);

    public native double getSumOfIncrements(int instance);

    public native void storeState(int instance);

    public native void restoreState(int instance);

    public native void acceptState(int instance);

    public native void makeDirty(int instance);

    public native void setPairwiseData(int instance, double[] observations);

    public native void setParameters(int instance, double[] parameters);

    public native double[] getPairwiseData(int instance);

    public native void getLocationGradient(int instance, double[] gradient);

    public native int getInternalDimension(int instance);

//jsize size = env->GetArrayLength( arr );
//std::vector<double> input( size );
//env->GetDoubleArrayRegion( arr, 0, size, &input[0] );
//
////  ...
//
//jdoubleArray output = env->NewDoubleArray( results.size() );
//env->SetDoubleArrayRegion( output, 0, results.size(), &results[0] );

    public class ResourceDetails {
        private final int number;
        private String name;
        private String description;
        private long flags;

        public ResourceDetails(int number) {
            this.number = number;
        }

        public int getNumber() {
            return this.number;
        }

        public String getName() {
            return this.name;
        }

        public void setName(String name) {
            this.name = name;
        }

        public String getDescription() {
            return this.description;
        }

        public void setDescription(String description) {
            this.description = description;
        }

        public long getFlags() {
            return this.flags;
        }

        public void setFlags(long flags) {
            this.flags = flags;
        }

        public String toString() {
            StringBuilder sb = new StringBuilder();
            sb.append("").append(getNumber()).append(" : ").append(getName()).append("\n");
            if (this.getDescription() != null) {
                sb.append("    ").append(getDescription());
            }

            sb.append("    Flags:");
            sb.append(getFlags());
            sb.append("\n");
            return sb.toString();
        }
    }
}
