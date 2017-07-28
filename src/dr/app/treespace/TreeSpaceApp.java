/*
 * TreeSpaceApp.java
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

package dr.app.treespace;

import jam.framework.*;
import jam.mac.Utils;

import javax.swing.*;
import java.awt.*;
import java.io.File;

/**
 * @author Andrew Rambaut
 * @version $Id$
 */
public class TreeSpaceApp extends MultiDocApplication {

    public TreeSpaceApp(String nameString, String aboutString, Icon icon,
                        String websiteURLString, String helpURLString) {
        super(new TreeSpaceMenuBarFactory(), nameString, aboutString, icon, websiteURLString, helpURLString);
    }

    /**
     * In a departure from the standard UI, there is no "Open" command for this application
     * This overriden method should never be called and throw a RuntimeException if it is.
     *
     * @return the action
     */
    public Action getOpenAction() {
        throw new UnsupportedOperationException("getOpenAction is not supported");
    }

    // Main entry point
    static public void main(String[] args) {

        boolean lafLoaded = false;

        if (Utils.isMacOSX()) {
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

            }

            UIManager.put("SystemFont", new Font("Lucida Grande", Font.PLAIN, 13));
            UIManager.put("SmallSystemFont", new Font("Lucida Grande", Font.PLAIN, 11));
        }

        try {

            if (!lafLoaded) {
                UIManager.setLookAndFeel(UIManager.getSystemLookAndFeelClassName());
            }

            java.net.URL url = TreeSpaceApp.class.getResource("images/phylogeography.png");
            Icon icon = null;

            if (url != null) {
                icon = new ImageIcon(url);
            }

            final String nameString = "TreeSpace";
            final String versionString = "v1.0";
            String aboutString = "<html><div style=\"font-family:sans-serif;\"><center>" +
                    "<div style=\"font-size:12;\"><p>Phylogenetic Tree Space Exploration<br>" +
                    "Version " + versionString + ", 2009</p>" +
                    "<p>by Andrew Rambaut</p></div>" +
                    "</div></html>";

            String websiteURLString = "http://beast.community";
            String helpURLString = "http://beast.community/phylogeography";

            TreeSpaceApp app = new TreeSpaceApp(nameString, aboutString, icon,
                    websiteURLString, helpURLString);
            app.setDocumentFrameFactory(new DocumentFrameFactory() {
                public DocumentFrame createDocumentFrame(final Application application, final MenuBarFactory menuBarFactory) {
                    return new TreeSpaceFrame(nameString);
                }
            });
            app.initialize();
            app.doNew();

            if (args.length > 0) {
                for (String fileName : args) {
                    ((TreeSpaceFrame)app.getUpperDocumentFrame()).importDataFile(new File(fileName));
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
