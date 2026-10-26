/*
 * TreeLimeApp.java
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

package dr.app.treelime;

import dr.app.beast.BeastMain;
import dr.app.beast.BeastVersion;
import dr.app.util.Arguments;
import dr.app.util.OSType;
import dr.util.Version;
import jam.framework.*;
import jam.util.IconUtils;

import javax.swing.*;
import java.awt.*;
import java.util.Locale;


/**
 * @author Andrew Rambaut
 * @author Alexei Drummond
 */
public class TreeLimeApp extends MultiDocApplication {
    public final static Version VERSION = new BeastVersion();

    public final static String backgroundColor =  "#35484F";
    public final static String foregroundColor = "#CBB944";

    public TreeLimeApp(String nameString, String aboutString, Icon icon,
                     String websiteURLString, String helpURLString) {
        super(new TreeLimeMenuBarFactory(), nameString, nameString, aboutString, Color.decode(backgroundColor), Color.decode(foregroundColor), icon, websiteURLString, helpURLString);
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
        centreLine("TreeLime " + VERSION.getVersionString() + ", " + VERSION.getDateString(), 60);
        centreLine("Maximum likelihood and MAP estimation for BEAST models", 60);
        for (String creditLine : VERSION.getCredits()) {
            centreLine(creditLine, 60);
        }
        System.out.println();
    }

    public static void printUsage(Arguments arguments) {

        arguments.printUsage("treelime", "[<input-file-name> ...]");
        System.out.println();
        System.out.println("  Example: treelime test.nex");
        System.out.println("  Example: treelime -help");
        System.out.println();
    }


//    /**
//     * In a departure from the standard UI, there is no "Open" command for this application
//     * Instead, the user can create a New window, Import a NEXUS file and Apply a Template file.
//     * None of these operations result in a file being associated with the DocumentFrame. All
//     * these actions are located in the TreeLimeFrame class. This overriden method should never
//     * be called and throw a RuntimeException if it is.
//     *
//     * @return the action
//     */
//    public Action getOpenAction() {
//        throw new UnsupportedOperationException("getOpenAction is not supported");
//    }

    private static boolean lafLoaded = false;

    /**
     * Sets the Dock (or taskbar) icon where the platform has one. java.awt.Taskbar is Java 9 and later, and the build
     * targets Java 8, so it is reached by reflection; without it the icon stays the default.
     */
    private static void setDockIcon(Image image) {
        if (image == null) {
            return;
        }
        try {
            Class<?> taskbarClass = Class.forName("java.awt.Taskbar");
            if (!(Boolean) taskbarClass.getMethod("isTaskbarSupported").invoke(null)) {
                return;
            }
            Object taskbar = taskbarClass.getMethod("getTaskbar").invoke(null);
            Class<?> featureClass = Class.forName("java.awt.Taskbar$Feature");
            Object iconFeature = featureClass.getField("ICON_IMAGE").get(null);
            if ((Boolean) taskbarClass.getMethod("isSupported", featureClass).invoke(taskbar, iconFeature)) {
                taskbarClass.getMethod("setIconImage", Image.class).invoke(taskbar, image);
            }
        } catch (Exception e) {
            // Java 8, or a platform without a taskbar icon
        }
    }

    // Main entry point
    static public void main(String[] args) {

        // There is a major issue with languages that use the comma as a decimal separator.
        // To ensure compatibility between programs in the package, enforce the US locale.
        Locale.setDefault(Locale.US);

        Arguments arguments = new Arguments(
                new Arguments.Option[]{
                        new Arguments.Option("advanced", "a", "Enable advanced & developer features"),
                        new Arguments.Option("version", "v", "Print the version and credits and stop"),
                        new Arguments.Option("help", "h", "Print this information and stop"),
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
            System.setProperty("apple.awt.antialiasing", "true");
            System.setProperty("apple.awt.rendering", "VALUE_RENDER_QUALITY");

            System.setProperty("apple.laf.useScreenMenuBar", "true");
            System.setProperty("apple.awt.draggableWindowBackground", "true");
            System.setProperty("apple.awt.showGrowBox", "true");

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

            final Image iconImage = IconUtils.getImage(TreeLimeApp.class, "images/treelime.png");
            javax.swing.Icon icon = IconUtils.resize(IconUtils.getIcon(TreeLimeApp.class, "images/treelime.png"), 128, 128);
            setDockIcon(iconImage);

            final String nameString = "TreeLime";
            final String versionString = VERSION.getVersionString();

            String aboutString = "<html>" +
                    "<div style=\"font-family:HelveticaNeue-Light, 'Helvetica Neue Light', Helvetica, Arial, 'Lucida Grande',sans-serif; " +
                    "font-weight: 100;" +
                    "background-color: " + backgroundColor + "; color: " + foregroundColor + ";" +
                    "\">" +
                    "<center>" +
                    "<div style=\"font-size:13\"><p>Maximum likelihood and MAP estimation for BEAST models<br>" +
                    "Version " + versionString + ", " + VERSION.getDateString() + "</p>" +
                    "<hr><div style=\"font-size:11;\">Part of the BEAST package:" +
                    VERSION.getHTMLCredits() +
                    "</div></center></div></html>";

            String websiteURLString = "http://beast.community/BEAUti";
            String helpURLString = "http://beast.community/BEAUti";

            System.setProperty("BEAST & BEAUTi Version", VERSION.getVersion());

            TreeLimeApp app = new TreeLimeApp(nameString, aboutString, icon,
                    websiteURLString, helpURLString);
            app.setDocumentFrameFactory(new DocumentFrameFactory() {
                public DocumentFrame createDocumentFrame(Application app, MenuBarFactory menuBarFactory) {
                    TreeLimeFrame frame = new TreeLimeFrame(nameString);
                    frame.setIconImage(iconImage); // the window's icon on Windows and Linux
                    return frame;
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
