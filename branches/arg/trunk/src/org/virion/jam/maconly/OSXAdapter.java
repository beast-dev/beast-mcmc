/*
 * OSXAdapter.java
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

package org.virion.jam.maconly;

import com.apple.eawt.ApplicationAdapter;
import com.apple.eawt.ApplicationEvent;

public class OSXAdapter extends ApplicationAdapter {

    // pseudo-singleton model; no point in making multiple instances
    // of the EAWT application or our adapter
    private static OSXAdapter theAdapter;
    private static com.apple.eawt.Application theApplication;

    // reference to the app where the existing quit, about, prefs code is
    private org.virion.jam.framework.Application application;

    private OSXAdapter(org.virion.jam.framework.Application application) {
        this.application = application;
    }

    // implemented handler methods.  These are basically hooks into existing
    // functionality from the main app, as if it came over from another platform.
    public void handleAbout(ApplicationEvent ae) {
        if (application != null) {
            ae.setHandled(true);
            application.doAbout();
        } else {
            throw new IllegalStateException("handleAbout: MyApp instance detached from listener");
        }
    }

    public void handlePreferences(ApplicationEvent ae) {
        if (application != null) {
            application.doPreferences();
            ae.setHandled(true);
        } else {
            throw new IllegalStateException("handlePreferences: MyApp instance detached from listener");
        }
    }

    public void handleQuit(ApplicationEvent ae) {
        if (application != null) {
            /*
            /	You MUST setHandled(false) if you want to delay or cancel the quit.
            /	This is important for cross-platform development -- have a universal quit
            /	routine that chooses whether or not to quit, so the functionality is identical
            /	on all platforms.  This example simply cancels the AppleEvent-based quit and
            /	defers to that universal method.
            */
            ae.setHandled(false);
            application.doQuit();
        } else {
            throw new IllegalStateException("handleQuit: MyApp instance detached from listener");
        }
    }


    // The main entry-point for this functionality.  This is the only method
    // that needs to be called at runtime, and it can easily be done using
    // reflection.
    public static void registerMacOSXApplication(org.virion.jam.framework.Application application) {
        if (theApplication == null) {
            theApplication = new com.apple.eawt.Application();
        }

        if (theAdapter == null) {
            theAdapter = new OSXAdapter(application);
        }
        theApplication.addApplicationListener(theAdapter);
    }

    // Another static entry point for EAWT functionality.  Enables the
    // "Preferences..." menu item in the application menu.
    public static void enablePrefs(boolean enabled) {
        if (theApplication == null) {
            theApplication = new com.apple.eawt.Application();
        }
        theApplication.setEnabledPreferencesMenu(enabled);
    }
}