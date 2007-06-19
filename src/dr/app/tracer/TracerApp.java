/*
 * TracerApp.java
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

package dr.app.tracer;

import org.virion.jam.framework.*;
import org.virion.jam.util.IconUtils;
import javax.swing.*;

public class TracerApp {
	public TracerApp() {
		try {
			
			UIManager.setLookAndFeel(UIManager.getSystemLookAndFeelClassName());

            Icon icon = IconUtils.getIcon(TracerApp.class, "images/Tracer.png");
    		
			String nameString = "Tracer";
			String aboutString = "MCMC Trace Analysis Tool\nVersion 1.3\n \nCopyright 2003-2005 Andrew Rambaut and Alexei Drummond\nUniversity of Oxford\nAll Rights Reserved.";

			MultiDocApplication app = new MultiDocApplication(nameString, aboutString, icon);

            Application.getMenuBarFactory().registerMenuFactory(new AnalysisMenuFactory());

            app.setDocumentFrameFactory(new DocumentFrameFactory() {
                public DocumentFrame createDocumentFrame(Application app, MenuBarFactory menuBarFactory) {
                    return new TracerFrame("Tracer");
                }
            });

            app.initialize();
			app.doNew();
		} catch (Exception e) {
			e.printStackTrace();
		}
	}

	// Main entry point
	static public void main(String[] args) {
		System.setProperty("apple.laf.useScreenMenuBar","true");
		System.setProperty("apple.awt.showGrowBox","true");
		System.setProperty("apple.awt.antialiasing","on");
		System.setProperty("apple.awt.textantialiasing","on");
		
		new TracerApp();
	}
	
}