/*
 * Utils.java
 *
 * Copyright (C) 2002-2006 Alexei Drummond and Andrew Rambaut
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

package org.virion.jam.mac;

import java.lang.reflect.Method;

public class Utils {

    protected static boolean MAC_OS_X;
    protected static String MAC_OS_X_VERSION;

    public static boolean isMacOSX() {
        return MAC_OS_X;
    }

    public static void macOSXRegistration(org.virion.jam.framework.Application application) {
        if (MAC_OS_X) {

            Class osxAdapter = null;

            try {
                osxAdapter = Class.forName("org.virion.jam.maconly.OSXAdapter");
            } catch (Exception e) {
                System.err.println("This version of Mac OS X does not support the Apple EAWT.");
            }

            try {
                if (osxAdapter != null) {

                    Class[] defArgs = {org.virion.jam.framework.Application.class};
                    Method registerMethod = osxAdapter.getDeclaredMethod("registerMacOSXApplication", defArgs);

                    if (registerMethod != null) {
                        Object[] args = {application};
                        registerMethod.invoke(osxAdapter, args);
                    }

                    // This is slightly gross.  to reflectively access methods with boolean args,
                    // use "boolean.class", then pass a Boolean object in as the arg, which apparently
                    // gets converted for you by the reflection system.
                    defArgs[0] = boolean.class;
                    Method prefsEnableMethod = osxAdapter.getDeclaredMethod("enablePrefs", defArgs);
                    if (prefsEnableMethod != null) {
                        Object args[] = {Boolean.TRUE};
                        prefsEnableMethod.invoke(osxAdapter, args);
                    }
                }

            } catch (Exception e) {
                System.err.println("Exception while loading the OSXAdapter:");
                e.printStackTrace();
            }
        }
    }

    static {
        MAC_OS_X_VERSION = System.getProperty("mrj.version");
        MAC_OS_X = MAC_OS_X_VERSION != null;
    }
}