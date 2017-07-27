/*
 * TreeStatApp.java
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

package dr.app.treestat;

import dr.app.beast.BeastVersion;
import dr.app.util.OSType;
import dr.evolution.io.NexusImporter;
import dr.util.Version;
import jam.framework.SingleDocApplication;
import jam.mac.Utils;

import javax.swing.*;
import java.awt.*;
import java.util.Locale;

public class TreeStatApp extends SingleDocApplication {

    private final static Version version = new BeastVersion();

    public TreeStatApp(String nameString, String aboutString, Icon icon, String websiteURLString, String helpURLString) {
        super(nameString, aboutString, icon, websiteURLString, helpURLString);
    }

    // Main entry point
    static public void main(String[] args) {

        // There is a major issue with languages that use the comma as a decimal separator.
        // To ensure compatibility between programs in the package, enforce the US locale.
        Locale.setDefault(Locale.US);

        // don't display warnings from NexusImporter as they are only relevant to Beast-MCMC-1.x
        NexusImporter.setSuppressWarnings(true);


        if (OSType.isMac()) {
            if (Utils.getMacOSXVersion().startsWith("10.5")) {
                System.setProperty("apple.awt.brushMetalLook","true");
            }
            System.setProperty("apple.laf.useScreenMenuBar","true");
            System.setProperty("apple.awt.showGrowBox","true");
            UIManager.put("SystemFont", new Font("Lucida Grande", Font.PLAIN, 13));
            UIManager.put("SmallSystemFont", new Font("Lucida Grande", Font.PLAIN, 11));
        }

        try {

            UIManager.setLookAndFeel(UIManager.getSystemLookAndFeelClassName());

            java.net.URL url = TreeStatApp.class.getResource("images/TreeStat.png");
            Icon icon = null;

            if (url != null) {
                icon = new ImageIcon(url);
            }

            final String versionString = version.getVersionString();
            String nameString = "TreeStat " + versionString;
            String aboutString = "<html><center><p>Tree Statistic Calculation Tool<br>" +
                    versionString + ", " + version.getDateString() + "</p>" +
//                    "Version 1.2, 2005-2010</p>" +
                    "<p>by<br>" +
                    "Andrew Rambaut and Alexei J. Drummond</p>" +
                    "<p>Institute of Evolutionary Biology, University of Edinburgh<br>" +
                    "<a href=\"mailto:a.rambaut@ed.ac.uk\">a.rambaut@ed.ac.uk</a></p>" +
                    "<p>Department of Computer Science, University of Auckland<br>" +
                    "<a href=\"mailto:alexei@cs.auckland.ac.nz\">alexei@cs.auckland.ac.nz</a></p>" +
                    "<p>Visit the BEAST page:<br>" +
                    "<a href=\"http://beast.community\">http://beast.community</a></p>" +
                    "</center></html>";

            String websiteURLString = "http://beast.community";
            String helpURLString = "http://beast.community/TreeStat";

            TreeStatApp app = new TreeStatApp(nameString, aboutString, icon,
                    websiteURLString, helpURLString);
            app.setDocumentFrame(new TreeStatFrame(app, nameString));
            app.initialize();

        } catch (Exception e) {
            JOptionPane.showMessageDialog(new JFrame(), "Fatal exception: " + e,
                    "Please report this to the authors",
                    JOptionPane.ERROR_MESSAGE);
            e.printStackTrace();
        }
    }
}
