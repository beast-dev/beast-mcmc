/*
 * TracerApp.java
 *
 * Copyright (c) 2002-2013 Alexei Drummond, Andrew Rambaut and Marc Suchard
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

package dr.app.tracer.application;

import dr.app.util.OSType;
import dr.inference.trace.LogFileTraces;
import jam.framework.*;

import javax.swing.*;
import java.awt.*;
import java.io.File;
import java.util.Locale;

public class TracerApp extends MultiDocApplication {

    public TracerApp(String nameString, String aboutString, Icon icon,
                     String websiteURLString, String helpURLString) {
        super(new TracerMenuBarFactory(), nameString, aboutString, icon, websiteURLString, helpURLString);

        addPreferencesSection(new GeneralPreferencesSection());
    }

    private static boolean lafLoaded = false;

    // Main entry point
    static public void main(String[] args) {
        // There is a major issue with languages that use the comma as a decimal separator.
        // To ensure compatibility between programs in the package, enforce the US locale.
        Locale.setDefault(Locale.US);

        if (OSType.isMac()) {
            System.setProperty("apple.awt.graphics.UseQuartz", "true");
            System.setProperty("apple.awt.antialiasing", "true");
            System.setProperty("apple.awt.rendering", "VALUE_RENDER_QUALITY");

            System.setProperty("apple.laf.useScreenMenuBar", "true");
            System.setProperty("apple.awt.draggableWindowBackground", "true");
            System.setProperty("apple.awt.showGrowBox", "true");

            try {
                // set the Quaqua Look and Feel in the UIManager
                javax.swing.SwingUtilities.invokeAndWait(new Runnable() {
                    public void run() {
                        try {
                            UIManager.setLookAndFeel(
                                    "ch.randelshofer.quaqua.QuaquaLookAndFeel"
                            );

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

        if (!lafLoaded) {

            try {
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

        try {
            java.net.URL url = TracerApp.class.getResource("images/Tracer.png");
            Icon icon = null;

            if (url != null) {
                icon = new ImageIcon(url);
            }

            final String nameString = "Tracer";
            final String versionString = "v1.6.1pre";
            String aboutString = "<html><center><p>MCMC Trace Analysis Tool<br>" +
                    "Version " + versionString + ", 2003-2014</p>" +
                    "<p>by<br>" +

                    "Andrew Rambaut, Marc A. Suchard, Walter Xie and Alexei J. Drummond</p>" +

                    "<p>Institute of Evolutionary Biology, University of Edinburgh<br>" +
                    "<a href=\"mailto:a.rambaut@ed.ac.uk\">a.rambaut@ed.ac.uk</a></p>" +

                    "<p>Departments of Biomathematics, Biostatistics and Human Genetics, UCLA<br>" +
                    "<a href=\"mailto:msuchard@ucla.edu\">msuchard@ucla.edu</a></p>" +

                    "<p>Department of Computer Science, University of Auckland<br>" +
                    "<a href=\"mailto:alexei@cs.auckland.ac.nz\">alexei@cs.auckland.ac.nz</a></p>" +

                    "<p>Available from the BEAST site:<br>" +
                    "<a href=\"http://beast.bio.ed.ac.uk/\">http://beast.bio.ed.ac.uk/</a></p>" +
                    "<p>Source code distributed under the GNU LGPL:<br>" +
                    "<a href=\"http://github.com/beast-dev/beast-mcmc/\">http://github.com/beast-dev/beast-mcmc/</a></p>" +
                    "<p>Thanks for contributions to: Joseph Heled, Oliver Pybus & Benjamin Redelings</p>" +
                    "</center></html>";

            String websiteURLString = "http://beast.bio.ed.ac.uk/";
            String helpURLString = "http://beast.bio.ed.ac.uk/Tracer";

            TracerApp app = new TracerApp(nameString, aboutString, icon, websiteURLString, helpURLString);
            app.setDocumentFrameFactory(new DocumentFrameFactory() {
                public DocumentFrame createDocumentFrame(Application app, MenuBarFactory menuBarFactory) {
                    return new TracerFrame("Tracer");
                }
            });
            app.initialize();

            app.doNew();

            if (args.length > 0) {
                TracerFrame frame = (TracerFrame) app.getDefaultFrame();
                for (String fileName : args) {

                    File file = new File(fileName);
                    LogFileTraces[] traces = {new LogFileTraces(fileName, file)};

                    frame.processTraces(traces);
                }
            }
        } catch (Exception e) {
            JOptionPane.showMessageDialog(new JFrame(), "Fatal exception: " + e,
                    "Please report this to the authors",
                    JOptionPane.ERROR_MESSAGE);
            e.printStackTrace();
        }
    }
}