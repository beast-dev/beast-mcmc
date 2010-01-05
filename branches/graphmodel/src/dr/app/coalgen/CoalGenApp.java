/*
 * CoalGenApp.java
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

package dr.app.coalgen;

import jam.framework.SingleDocApplication;

import javax.swing.*;

/**
 * @author Andrew Rambaut
 * @author Alexei Drummond
 * @version $Id$
 */
public class CoalGenApp {
    public CoalGenApp() {

        try {

            UIManager.setLookAndFeel(UIManager.getSystemLookAndFeelClassName());

            java.net.URL url = CoalGenApp.class.getResource("/images/CoalGen.png");
            Icon icon = null;

            if (url != null) {
                icon = new ImageIcon(url);
            }

            String nameString = "CoalGen";
            String aboutString = "Coalescent Tree Simulator\nVersion 2.0\n \nCopyright 2004-2009 Andrew Rambaut and Alexei Drummond\nAll Rights Reserved.";


            SingleDocApplication app = new SingleDocApplication(new CoalGenMenuFactory(), nameString, aboutString, icon);


            CoalGenFrame frame = new CoalGenFrame(nameString);
            app.setDocumentFrame(frame);

        } catch (Exception e) {
            e.printStackTrace();
        }
    }

    // Main entry point
    static public void main(String[] args) {
        System.setProperty("com.apple.macos.useScreenMenuBar", "true");
        System.setProperty("apple.laf.useScreenMenuBar", "true");

        new CoalGenApp();
    }
}