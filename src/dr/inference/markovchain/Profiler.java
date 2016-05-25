/*
 * Profiler.java
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

package dr.inference.markovchain;

import java.util.HashMap;
import java.util.Iterator;
import java.util.Map;

/**
 * @author Andrew Rambaut
 *
 * @version $Id: Profiler.java,v 1.3 2004/12/16 10:25:01 alexei Exp $
 */
public class Profiler {

    private static boolean profilerAvailable = true;

    private static Map<String, Profile> profiles = new HashMap<String, Profile>();
    //private static String currentProfile = null;
    private static long startTime;

    public static boolean startProfile() { //String name) {
        if (profilerAvailable) {
            //@todo unfortunately I can't get my Ant Build to work with asserts yet. Help!
            //assert currentProfile == null : "Profile " + currentProfile + " is already recording";
            //currentProfile = name;
            startTime = getCurrentThreadCpuTime();
        }
        return true;
    }

    public static boolean stopProfile(String name) {
        if (profilerAvailable) {
            long stopTime = getCurrentThreadCpuTime();
            //@todo unfortunately I can't get my Ant Build to work with asserts yet. Help!
            //assert name.equals(currentProfile) : "Profile " + name + " is not recording";

            long count = 1;
            long totalTime = stopTime - startTime;
            Profile profile = profiles.get(name);
            if (profile != null) {
                totalTime += profile.time;
                count += profile.count;
            }

            profiles.put(name, new Profile(count, totalTime));
            //currentProfile = null;
        }
        return true;
    }

    public static void report() {
        if (profilerAvailable) {
            Iterator<String> iter = profiles.keySet().iterator();
            while (iter.hasNext()) {
                String name = iter.next();
                Profile profile = profiles.get(name);
                long average = profile.time / profile.count;
                System.err.println("PROFILE: " + name + " [" + profile.time + " ms, " + profile.count + " calls, " + average + " ms / call]");
            }
        }
    }

    private static class Profile {
        long time;
        long count;

        Profile(long count, long time) {
            this.count = count;
            this.time = time;
        }
    }

    public static native long getCurrentThreadCpuTime();

    static {
        try {
            System.loadLibrary("mcmcprof");
        } catch (UnsatisfiedLinkError ule) {
            profilerAvailable = false;
        }
    }
}
