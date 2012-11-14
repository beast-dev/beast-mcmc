/*
 * InstantiableTracerApp.java
 *
 * Copyright (C) 2002-2012 Alexei Drummond and Andrew Rambaut
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
 * BEAST is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with BEAST; if not, write to the
 * Free Software Foundation, Inc., 51 Franklin St, Fifth Floor,
 * Boston, MA  02110-1301  USA
 */

package dr.app.tracer.application;

import dr.app.util.OSType;
import dr.inference.trace.LogFileTraces;
import jam.framework.DocumentFrame;
import jam.framework.SingleDocApplication;

import javax.swing.*;
import java.awt.*;
import java.io.File;
import java.security.Permission;
import java.util.Locale;

/**
 *
 * A class for calling an instance of Tracer within a BEAST (or any other Java) analysis.
 * Allows a threaded call of Tracer (that won't terminate the JVM). A lot of the code is
 * copied from TracerApp.java (but is a SingleDocApplication instead of a
 * MultiDocApplication)
 *
 *
 * @author Wai Lok Sibon Li
 * @version $Id: InstantiableTracerApp.java,v 1.0 2012/09/17 15:23:33 sibon.li Exp $
 */
public class InstantiableTracerApp extends SingleDocApplication {
    private boolean exiting;
    public InstantiableTracerApp (String nameString, String aboutString, Icon icon,
                         String websiteURLString, String helpURLString, boolean exiting) {
            super(new TracerMenuBarFactory(), nameString, aboutString, icon,
                    websiteURLString, helpURLString);

            addPreferencesSection(new GeneralPreferencesSection());
        this.exiting = exiting;
    }


    @Override
    public void doQuit() {
        if(exiting) {   // Implemented this way because documentFrame is private in SingleDocApplication
            super.doQuit();
        }
        else {
            if (documentFrame == null) {
                return;
            }
            if (documentFrame.requestClose()) {

                documentFrame.setVisible(false);
                documentFrame.dispose();
    //            try {
    //                System.exit(0);
    //            }catch (ExitException e) {
    //                System.setSecurityManager(null);
    //            }
            }
        }
    }

//    @Override
//    public void doQuit() {
//        //super.doQuit();
//        if (documentFrame == null) {
//            return;
//        }
//        if (documentFrame.requestClose()) {
//
//            documentFrame.setVisible(false);
//            documentFrame.dispose();
////            try {
////                System.exit(0);
////            }catch (ExitException e) {
////                System.setSecurityManager(null);
////            }
//        }
//    }

    private DocumentFrame documentFrame = null;


    private static void loadTracerInstance(String nameString, String logFileName, long bi, boolean exiting) {
        final String name = nameString;
        final String fileName = logFileName;
        final long burnin = bi;
        final boolean exit = exiting;

        Thread thread = new Thread() {
            public void run() {
                try {

                    // There is a major issue with languages that use the comma as a decimal separator.
                    // To ensure compatibility between programs in the package, enforce the US locale.
                    Locale.setDefault(Locale.US);

                    boolean lafLoaded = false;

                    if (OSType.isMac()) {
                        System.setProperty("apple.awt.graphics.UseQuartz", "true");
                        System.setProperty("apple.awt.antialiasing","true");
                        System.setProperty("apple.awt.rendering","VALUE_RENDER_QUALITY");

                        System.setProperty("apple.laf.useScreenMenuBar","true");
                        System.setProperty("apple.awt.draggableWindowBackground","true");
                        System.setProperty("apple.awt.showGrowBox","true");

                        // set the Quaqua Look and Feel in the UIManager
                        try {
                            UIManager.setLookAndFeel(
                                    "ch.randelshofer.quaqua.QuaquaLookAndFeel"
                            );
                            lafLoaded = true;


                        } catch (Exception e) {
                            //
                        }

                        UIManager.put("SystemFont", new Font("Lucida Grande", Font.PLAIN, 13));
                        UIManager.put("SmallSystemFont", new Font("Lucida Grande", Font.PLAIN, 11));
                    }

                    if (!lafLoaded) {
                        try {
                            UIManager.setLookAndFeel(UIManager.getSystemLookAndFeelClassName());
                        } catch (Exception e1) {
                            e1.printStackTrace();
                        }
                    }

                    java.net.URL url = TracerApp.class.getResource("images/Mapper.png");
                    Icon icon = null;

                    if (url != null) {
                        icon = new ImageIcon(url);
                    }

                    InstantiableTracerApp app = new InstantiableTracerApp(name,
                            "Tracer tool running through "+ name + ". Authors Wai Lok Sibon Li & Andrew Rambaut",
                            icon, "http://beast.bio.ed.ac.uk/", "http://beast.bio.ed.ac.uk/Tracer", exit);

                    TracerFrame frame = new TracerFrame(name);
                    app.setDocumentFrame(frame);
                    app.initialize();

                    File file = new File(fileName);
                    LogFileTraces[] traces = { new LogFileTraces(fileName, file) };
                    traces[0].setBurnIn((int) burnin);
                    frame.processTraces(traces);
                }
                catch(Exception e) {
                    e.printStackTrace();
                    throw new RuntimeException("Exception in InstantiableTracerApp");
                }
            }
        };
//        System.setSecurityManager(new NoExitSecurityManager());
        thread.start();

    }

    /*
     *Instantiates a Tracer window which will terminate the JRE (similar to Tracer but burnin can be set)
     */
    public static void loadExitingTracerInstance(String nameString, String logFileName, long bi) {
        loadTracerInstance(nameString, logFileName, bi, true);
    }

    /*
     * Loads an instance of tracer that does not terminate the JRE. Also does not call System.exit() when the
     * window is closed by the user.
     */
    public static void loadNonExitingTracerInstance(String nameString, String logFileName, long bi) {
        loadTracerInstance(nameString, logFileName, bi, false);
    }

}

@Deprecated
class NoExitSecurityManager extends SecurityManager {
    @Override
    public void checkPermission(Permission perm) {
//            super.checkPermission(perm);
    }
    @Override
    public void checkPermission(Permission perm, Object context) {
//            super.checkPermission(perm, context);
    }
    @Override
    public void checkExit(int status) {
        super.checkExit(status);
        throw new ExitException(status);//ExitException(status);

    }
}

@Deprecated
class ExitException extends SecurityException {
    public final int status;
    public ExitException(int status) {
        super();
        this.status = status;
    }
}