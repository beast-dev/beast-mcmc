/*
 * MacDesktopIntegration.java
 *
 * Copyright © 2002-2024 the BEAST Development Team
 * http://beast.community/about
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
 *
 */

package dr.app.util;

import jam.framework.Application;

import javax.swing.*;
import java.awt.*;
import java.io.File;
import java.lang.reflect.InvocationHandler;
import java.lang.reflect.Method;
import java.lang.reflect.Proxy;
import java.util.List;

/**
 * Hooks a jam {@link Application} into the macOS application menu and Dock: the About and Preferences
 * items, Cmd-Q, files opened from the Finder or dropped on the Dock, the menu bar shown when no window is
 * open, and the Dock icon. It uses java.awt.Desktop and java.awt.Taskbar (Java 9), through reflection so
 * that the code still compiles for Java 8; jam's own adapter needs com.apple.eawt, which modern JDKs no
 * longer have, so without this the generic Java handlers are used.
 *
 * @author Andrew Rambaut
 */
public final class MacDesktopIntegration {

    private MacDesktopIntegration() {
    }

    /**
     * Registers the application's handlers; does nothing (apart from a note on stderr) where the Java 9
     * desktop API is not available.
     *
     * @param app      the application whose doAbout, doPreferences, doQuit and doOpenFile are to be called
     * @param dockIcon the image for the Dock (an ImageIcon), or null to leave the Dock icon alone
     */
    public static void register(final Application app, Icon dockIcon) {
        try {
            if (!Desktop.isDesktopSupported()) {
                return;
            }
            final Desktop desktop = Desktop.getDesktop();

            setHandler(desktop, "setAboutHandler", "java.awt.desktop.AboutHandler", new InvocationHandler() {
                public Object invoke(Object proxy, Method method, Object[] args) {
                    app.doAbout();
                    return null;
                }
            });
            setHandler(desktop, "setPreferencesHandler", "java.awt.desktop.PreferencesHandler", new InvocationHandler() {
                public Object invoke(Object proxy, Method method, Object[] args) {
                    app.doPreferences();
                    return null;
                }
            });
            setHandler(desktop, "setQuitHandler", "java.awt.desktop.QuitHandler", new InvocationHandler() {
                public Object invoke(Object proxy, Method method, Object[] args) throws Exception {
                    // doQuit exits the application unless the user cancels, in which case the quit is cancelled
                    app.doQuit();
                    Object response = args[1];
                    response.getClass().getMethod("cancelQuit").invoke(response);
                    return null;
                }
            });
            setHandler(desktop, "setOpenFileHandler", "java.awt.desktop.OpenFilesHandler", new InvocationHandler() {
                @SuppressWarnings("unchecked")
                public Object invoke(Object proxy, Method method, Object[] args) throws Exception {
                    Object event = args[0];
                    List<File> files = (List<File>) event.getClass().getMethod("getFiles").invoke(event);
                    for (File file : files) {
                        app.doOpenFile(file);
                    }
                    return null;
                }
            });

            // the menu bar shown when no window is open
            JMenuBar menuBar = new JMenuBar();
            Application.getMenuBarFactory().populateMenuBar(menuBar, null);
            Desktop.class.getMethod("setDefaultMenuBar", JMenuBar.class).invoke(desktop, menuBar);

            setDockIcon(dockIcon);
        } catch (Exception e) {
            // (ReflectiveOperationException on Java 8; UnsupportedOperationException where a feature is missing)
            System.err.println("The macOS application menu could not be set up: " + e);
        }
    }

    /**
     * Registers only an About handler and the Dock icon: for a program that shows a dialog without a jam
     * Application (such as BEAST's options dialog).
     *
     * @param about    shows the about box
     * @param dockIcon the image for the Dock (an ImageIcon), or null to leave the Dock icon alone
     */
    public static void registerAbout(final Runnable about, Icon dockIcon) {
        try {
            if (!Desktop.isDesktopSupported()) {
                return;
            }
            setHandler(Desktop.getDesktop(), "setAboutHandler", "java.awt.desktop.AboutHandler", new InvocationHandler() {
                public Object invoke(Object proxy, Method method, Object[] args) {
                    about.run();
                    return null;
                }
            });
            setDockIcon(dockIcon);
        } catch (Exception e) {
            System.err.println("The macOS application menu could not be set up: " + e);
        }
    }

    private static void setDockIcon(Icon dockIcon) throws Exception {
        if (dockIcon instanceof ImageIcon) {
            Class<?> taskbarClass = Class.forName("java.awt.Taskbar");
            Object taskbar = taskbarClass.getMethod("getTaskbar").invoke(null);
            taskbarClass.getMethod("setIconImage", Image.class).invoke(taskbar, ((ImageIcon) dockIcon).getImage());
        }
    }

    /**
     * Calls desktop.<setter>(handler) with a proxy of the handler interface that forwards the interface's
     * calls (but not Object's equals, hashCode and toString) to the invocation handler.
     */
    private static void setHandler(Desktop desktop, String setter, String handlerInterface, final InvocationHandler handler)
            throws Exception {
        Class<?> type = Class.forName(handlerInterface);
        Object proxy = Proxy.newProxyInstance(type.getClassLoader(), new Class<?>[]{type}, new InvocationHandler() {
            public Object invoke(Object proxy, Method method, Object[] args) throws Throwable {
                if (method.getDeclaringClass() == Object.class) {
                    if (method.getName().equals("equals")) {
                        return proxy == args[0];
                    } else if (method.getName().equals("hashCode")) {
                        return System.identityHashCode(proxy);
                    } else {
                        return handler.toString();
                    }
                }
                return handler.invoke(proxy, method, args);
            }
        });
        Desktop.class.getMethod(setter, type).invoke(desktop, proxy);
    }
}
