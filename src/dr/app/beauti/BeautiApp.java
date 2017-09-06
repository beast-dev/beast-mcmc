/*
 * BeautiApp.java
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

package dr.app.beauti;

import dr.app.beast.BeastVersion;
import dr.app.util.Arguments;
import dr.app.util.OSType;
import dr.util.Version;
import jam.framework.*;

import javax.swing.*;
import java.awt.*;
import java.lang.reflect.Method;
import java.security.AccessControlException;
import java.util.HashSet;
import java.util.Locale;
import java.util.Set;


/**
 * @author Andrew Rambaut
 * @author Alexei Drummond
 * @version $Id: BeautiApp.java,v 1.18 2006/09/09 16:07:05 rambaut Exp $
 */
public class BeautiApp extends MultiDocApplication {
    public final static Version VERSION = new BeastVersion();

    public BeautiApp(String nameString, String aboutString, Icon icon,
                     String websiteURLString, String helpURLString) {
        super(new BeautiMenuBarFactory(), nameString, aboutString, icon, websiteURLString, helpURLString);
    }

    public static void centreLine(String line, int pageWidth) {
        int n = pageWidth - line.length();
        int n1 = n / 2;
        for (int i = 0; i < n1; i++) {
            System.out.print(" ");
        }
        System.out.println(line);
    }

    public static void printTitle() {
        System.out.println();
        centreLine("BEAUti " + VERSION.getVersionString() + ", " + VERSION.getDateString(), 60);
        centreLine("Bayesian Evolutionary Analysis Utility", 60);
        for (String creditLine : VERSION.getCredits()) {
            centreLine(creditLine, 60);
        }
        System.out.println();
    }

    public static void printUsage(Arguments arguments) {

        arguments.printUsage("beauti", "[<input-file-name> ...]");
        System.out.println();
        System.out.println("  Example: beauti test.nex");
        System.out.println("  Example: beauti -help");
        System.out.println();
    }


//    /**
//     * In a departure from the standard UI, there is no "Open" command for this application
//     * Instead, the user can create a New window, Import a NEXUS file and Apply a Template file.
//     * None of these operations result in a file being associated with the DocumentFrame. All
//     * these actions are located in the BeautiFrame class. This overriden method should never
//     * be called and throw a RuntimeException if it is.
//     *
//     * @return the action
//     */
//    public Action getOpenAction() {
//        throw new UnsupportedOperationException("getOpenAction is not supported");
//    }

    private static boolean lafLoaded = false;

    // Main entry point
    static public void main(String[] args) {

        // There is a major issue with languages that use the comma as a decimal separator.
        // To ensure compatibility between programs in the package, enforce the US locale.
        Locale.setDefault(Locale.US);

        Arguments arguments = new Arguments(
                new Arguments.Option[]{
                        new Arguments.Option("advanced", "Enable advanced & developer features"),
                        new Arguments.Option("version", "Print the version and credits and stop"),
                        new Arguments.Option("help", "Print this information and stop"),
                });

        int argumentCount = 0;

        try {
            argumentCount = arguments.parseArguments(args);
        } catch (Arguments.ArgumentException ae) {
            printTitle();
            System.out.println();
            System.out.println(ae.getMessage());
            System.out.println();
            printUsage(arguments);
            System.exit(1);
        }

        if (arguments.hasOption("help")) {
            printTitle();
            printUsage(arguments);
            System.exit(0);
        }

        if (arguments.hasOption("version")) {
            printTitle();
        }

        advanced = arguments.hasOption("advanced");

        String[] args2 = arguments.getLeftoverArguments();

        if (args2.length > 1) {
            System.err.println("Unknown option: " + args2[1]);
            System.err.println();
            printUsage(arguments);
            return;
        }

//        if (args2.length > 1) {
//
//            if (args.length != 3) {
//                printTitle();
//                printUsage(arguments);
//                System.exit(1);
//            }
//
//            String inputFileName = args[0];
//            String templateFileName = args[1];
//            String outputFileName = args[2];
//
//            new CommandLineBeauti(inputFileName, templateFileName, outputFileName);
//
//        } else {
        String inputFileName = null;
        if (args2.length == 1) {
            inputFileName = args2[0];
        }

        if (OSType.isMac()) {
            // Explicitly turn on font antialiasing.
            System.setProperty("swing.aatext", "true");

            System.setProperty("apple.awt.graphics.UseQuartz", "true");
            System.setProperty("apple.awt.antialiasing","true");
            System.setProperty("apple.awt.rendering","VALUE_RENDER_QUALITY");

            System.setProperty("apple.laf.useScreenMenuBar","true");
            System.setProperty("apple.awt.draggableWindowBackground","true");
            System.setProperty("apple.awt.showGrowBox","true");

            try {
                // set the Quaqua Look and Feel in the UIManager
                javax.swing.SwingUtilities.invokeAndWait(new Runnable() {
                    public void run() {
                        try {
//                                try {
//                                    // We need to do this using dynamic class loading to avoid other platforms
//                                    // having to link to this class. If the Quaqua library is not on the classpath
//                                    // it simply won't be used.
//                                    Class<?> qm = Class.forName("ch.randelshofer.quaqua.QuaquaManager");
//                                    Method method = qm.getMethod("setExcludedUIs", Set.class);
//
//                                    Set<String> excludes = new HashSet<String>();
//                                    excludes.add("Button");
//                                    excludes.add("ToolBar");
//                                    method.invoke(null, excludes);
//
//                                }
//                                catch (Throwable e) {
//                                }

                            //set the Quaqua Look and Feel in the UIManager
                            UIManager.setLookAndFeel("ch.randelshofer.quaqua.QuaquaLookAndFeel");
                            lafLoaded = true;
                        } catch (Exception e) {
                        }
                    }
                });
            } catch (Exception e) {
            }

            UIManager.put("SystemFont", new Font("Lucida Grande", Font.PLAIN, 13));
            UIManager.put("SmallSystemFont", new Font("Lucida Grande", Font.PLAIN, 11));
        }

        try {

            if (!lafLoaded) {
                try {
                    // set the System Look and Feel in the UIManager
                    javax.swing.SwingUtilities.invokeAndWait(new Runnable() {
                        public void run() {
                            try {
                                UIManager.setLookAndFeel(UIManager.getSystemLookAndFeelClassName());
                            } catch (Exception e) {
                                e.printStackTrace();
                            }
                        }
                    });
                } catch (Exception e) {
                    e.printStackTrace();
                }
            }

            java.net.URL url = BeautiApp.class.getResource("images/beauti.png");
            Icon icon = null;

            if (url != null) {
                icon = new ImageIcon(url);
            }

            final String nameString = "BEAUti";
            final String versionString = VERSION.getVersionString();

            String aboutString = "<html>" +
                    "<div style=\"font-family:HelveticaNeue-Light, 'Helvetica Neue Light', Helvetica, Arial, 'Lucida Grande',sans-serif; font-weight: 100\">" +
                    "<center>" +
                    "<div style=\"font-size:13\"><p>Bayesian Evolutionary Analysis Utility<br>" +
                    "Version " + versionString + ", " + VERSION.getDateString() + "</p>" +
                    "<p>by Alexei J. Drummond, Andrew Rambaut, Marc A. Suchard and Walter Xie</p></div>" +
                    "<hr><div style=\"font-size:11;\">Part of the BEAST package:" +
                    VERSION.getHTMLCredits() +
                    "</div></center></div></html>";

            String websiteURLString = "http://beast.community/BEAUti";
            String helpURLString = "http://beast.community/BEAUti";

            System.setProperty("BEAST & BEAUTi Version", VERSION.getVersion());

            BeautiApp app = new BeautiApp(nameString, aboutString, icon,
                    websiteURLString, helpURLString);
            app.setDocumentFrameFactory(new DocumentFrameFactory() {
                public DocumentFrame createDocumentFrame(Application app, MenuBarFactory menuBarFactory) {
                    return new BeautiFrame(nameString);
                }
            });
            app.initialize();

            if (inputFileName != null) {
                app.doOpen(inputFileName);
            } else {
                app.doNew();
            }
        } catch (Exception e) {
            JOptionPane.showMessageDialog(new JFrame(), "Fatal exception: " + e,
                    "Please report this to the authors",
                    JOptionPane.ERROR_MESSAGE);
            e.printStackTrace();
        }
//        }
    }

    public static boolean advanced = false;
}
