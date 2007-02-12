/*
 * TreeStatApp.java
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

package dr.app.treestat;

import javax.swing.*;

import org.virion.jam.framework.SingleDocApplication;

public class TreeStatApp {
    public TreeStatApp() {
        try {

            UIManager.setLookAndFeel(UIManager.getSystemLookAndFeelClassName());

            java.net.URL url = TreeStatApp.class.getResource("images/TreeStat.png");
            Icon icon = null;

            if (url != null) {
                icon = new ImageIcon(url);
            }

            String nameString = "TreeStat";
            String aboutString = "Tree Statistic Calculation Tool\nVersion 1.1\n \nCopyright 2005-2006 Andrew Rambaut and Alexei Drummond\nUniversity of Oxford\nAll Rights Reserved.";


            SingleDocApplication app = new SingleDocApplication(nameString, aboutString, icon);


            TreeStatFrame frame = new TreeStatFrame(app, nameString);
            app.setDocumentFrame(frame);

        } catch (Exception e) {
            e.printStackTrace();
        }
    }

    // Main entry point
    static public void main(String[] args) {
        System.setProperty("com.apple.macos.useScreenMenuBar","true");
        System.setProperty("apple.laf.useScreenMenuBar","true");

        new TreeStatApp();
    }

}