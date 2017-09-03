/*
 * TempestApp.java
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

package dr.app.tempest;

import dr.app.util.OSType;
import dr.util.Version;
import jam.framework.*;

import javax.swing.*;
import java.awt.*;

/**
 * @author Andrew Rambaut
 * @version $Id$
 */
public class TempEstApp extends MultiDocApplication {
    private final static Version version = new Version() {
        private static final String VERSION = "1.5.1";

        public String getVersion() {
            return VERSION;
        }

        public String getVersionString() {
            return "v" + VERSION;
        }

        public String getDateString() {
            return "2003-2017";
        }

        public String getBuildString() {
            return "GitHub 20160821";
        }

        public String[] getCredits() {
            return new String[0];
        }

        public String getHTMLCredits() {
            return "<p>by<br>" +
                    "Andrew Rambaut</p>" +
                    "<p>Institute of Evolutionary Biology, University of Edinburgh<br>" +
                    "<a href=\"mailto:a.rambaut@ed.ac.uk\">a.rambaut@ed.ac.uk</a></p>" +
                    "<p>Citation<br>" +
                    "<a href=\"http://dx.doi.org/10.1093/ve/vew007\">Rambaut, Lam, de Carvalho & Pybus (2016) Exploring the temporal structure of<br>" +
                    "heterochronous sequences using TempEst. <i>Virus Evolution</i> <b>2</b>: vew007</a></p>" +
                    "<p>Part of the BEAST package:<br>" +
                    "<a href=\"http://beast.community\">http://beast.community</a></p>";
        }

    };

    public TempEstApp(String nameString, String aboutString, Icon icon,
                      String websiteURLString, String helpURLString) {
        super(new TempestMenuBarFactory(), nameString, aboutString, icon, websiteURLString, helpURLString);
    }

    // Main entry point
    static public void main(String[] args) {


        if (OSType.isMac()) {
            System.setProperty("apple.laf.useScreenMenuBar","true");
            System.setProperty("apple.awt.showGrowBox","true");
            System.setProperty("apple.awt.graphics.UseQuartz","true");
            UIManager.put("SystemFont", new Font("Lucida Grande", Font.PLAIN, 13));
            UIManager.put("SmallSystemFont", new Font("Lucida Grande", Font.PLAIN, 11));
        }

        try {

            try {
                SwingUtilities.invokeAndWait(new Runnable() {
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

            java.net.URL url = TempEstApp.class.getResource("images/tempest.png");
            Icon icon = null;

            if (url != null) {
                icon = new ImageIcon(url);
            }

            final String nameString = "TempEst";
            final String versionString = version.getVersionString();
            String aboutString = "<html><center><p>Temporal Signal Estimator Tool<br>" +
                    "Version " + versionString + ", " + version.getDateString() + "</p>" +
                    version.getHTMLCredits() +
                    "</center></html>";

            String websiteURLString = "http://tree.bio.ed.ac.uk/";
            String helpURLString = "http://tree.bio.ed.ac.uk/software/tempest";

            TempEstApp app = new TempEstApp(nameString, aboutString, icon,
                    websiteURLString, helpURLString);
            app.setDocumentFrameFactory(new DocumentFrameFactory() {
                public DocumentFrame createDocumentFrame(Application app, MenuBarFactory menuBarFactory) {
                    return new TempestFrame(nameString);
                }
            });
            app.initialize();
            app.doOpen();
        } catch (Exception e) {
            JOptionPane.showMessageDialog(new JFrame(), "Fatal exception: " + e,
                    "Please report this to the authors",
                    JOptionPane.ERROR_MESSAGE);
            e.printStackTrace();
        }
    }

}