/*
 * BeautiApp.java
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

package dr.app.beauti;

import org.virion.jam.framework.SingleDocApplication;

import javax.swing.*;

/** 
 * @author			Andrew Rambaut	
 * @author			Alexei Drummond	
 * @version			$Id: BeautiApp.java,v 1.18 2006/09/09 16:07:05 rambaut Exp $
 */
public class BeautiApp extends SingleDocApplication {
	public BeautiApp(String nameString, String aboutString, Icon icon,
						String websiteURLString, String helpURLString) {
		super(nameString, aboutString, icon, websiteURLString, helpURLString);

		setDocumentFrame(new BeautiFrame(nameString));
	}

	// Main entry point
	static public void main(String[] args) {
		System.setProperty("com.apple.macos.useScreenMenuBar","true");
		System.setProperty("apple.laf.useScreenMenuBar","true");

		try {

			UIManager.setLookAndFeel(UIManager.getSystemLookAndFeelClassName());

			java.net.URL url = BeautiApp.class.getResource("images/beauti.png");
			Icon icon = null;

			if (url != null) {
            	icon = new ImageIcon(url);
    		}

			String nameString = "BEAUti";
			String aboutString = "Bayesian Evolutionary Analysis Utility\n" +
									"BEAST XML generation tool\n" +
									"Version 1.4\n \n" +
									"Copyright 2003-2006 Andrew Rambaut and Alexei Drummond\n" +
									"University of Oxford\n" +
									"All Rights Reserved.";
			String websiteURLString = "http://evolve.zoo.ox.ac.uk/beast/";
			String helpURLString = "http://evolve.zoo.ox.ac.uk/beast/help/";

			BeautiApp app = new BeautiApp(nameString, aboutString, icon,
											websiteURLString, helpURLString);
		} catch (Exception e) {
			JOptionPane.showMessageDialog(new JFrame(), "Fatal exception: " + e,
												"Please report this to the authors",
												JOptionPane.ERROR_MESSAGE);
			e.printStackTrace();
		}
	}

}