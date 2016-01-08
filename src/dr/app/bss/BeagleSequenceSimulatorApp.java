/*
 * BeagleSequenceSimulatorApp.java
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

package dr.app.bss;

import jam.framework.SingleDocApplication;

import java.awt.Font;
import java.awt.Image;
import java.awt.Toolkit;
import java.net.URL;
import java.util.Locale;

import javax.swing.ImageIcon;
import javax.swing.UIManager;
import javax.swing.UnsupportedLookAndFeelException;

/**
 * @author Filip Bielejec
 * @version $Id$
 */
public class BeagleSequenceSimulatorApp {

	public static final boolean DEBUG = false;
	public static final boolean VERBOSE = true;
	
	// Share those if neccessary
	public static final String SHORT_NAME = "\u03C0BUSS";
	public static final String LONG_NAME = "Parallel BEAST/BEAGLE Utility for Sequence Simulation";
	public static final String VERSION = "1.3.8rc";
	public static final String DATE = "2015";

	// Icons
	private Image beagleSequenceSimulatorImage;
	
	public BeagleSequenceSimulatorApp() throws ClassNotFoundException,
			InstantiationException, IllegalAccessException,
			UnsupportedLookAndFeelException {

		boolean lafLoaded = false;

		// Setup Look & Feel
		if (System.getProperty("os.name").toLowerCase().startsWith("mac os x")) {

			// Mac stuff
			System.setProperty("apple.awt.showGrowBox", "true");
			System.setProperty("apple.awt.brushMetalLook", "true");
			System.setProperty("apple.laf.useScreenMenuBar", "true");

			System.setProperty("apple.awt.graphics.UseQuartz", "true");
			System.setProperty("apple.awt.antialiasing", "true");
			System.setProperty("apple.awt.rendering", "VALUE_RENDER_QUALITY");

			System.setProperty("apple.laf.useScreenMenuBar", "true");
			System.setProperty("apple.awt.draggableWindowBackground", "true");
			System.setProperty("apple.awt.showGrowBox", "true");

			UIManager.put("SystemFont", new Font("Lucida Grande", Font.PLAIN,
					13));
			UIManager.put("SmallSystemFont", new Font("Lucida Grande",
					Font.PLAIN, 11));

			try {

				// UIManager.setLookAndFeel(UIManager
				// .getSystemLookAndFeelClassName());

				UIManager
						.setLookAndFeel("ch.randelshofer.quaqua.QuaquaLookAndFeel");
				lafLoaded = true;

			} catch (Exception e) {
				//
			}

		} else {

			try {

				// UIManager.setLookAndFeel(UIManager
				// .getSystemLookAndFeelClassName());

				UIManager
						.setLookAndFeel("com.sun.java.swing.plaf.nimbus.NimbusLookAndFeel");
				lafLoaded = true;

			} catch (Exception e) {
				//
			}

		}

		if (!lafLoaded) {

			try {

				UIManager.setLookAndFeel(UIManager
						.getSystemLookAndFeelClassName());
				System.out
						.println("Specified l&f not found. Loading system default l&f");

			} catch (Exception e) {

				e.printStackTrace();

			}
		}

		ImageIcon bssIcon = Utils.createImageIcon(Utils.BSS_ICON);
		SingleDocApplication app = new SingleDocApplication(new MenuBarFactory(), //
				LONG_NAME, //
				VERSION.concat(" ").concat(
				DATE), //
				bssIcon //
		);

		beagleSequenceSimulatorImage = CreateImage(Utils.BSS_ICON);
		MainFrame frame = new MainFrame(SHORT_NAME);
		frame.setIconImage(beagleSequenceSimulatorImage);
		app.setDocumentFrame(frame);

	}// END: Constructor

	public static void main(String args[]) {

		 Locale.setDefault(Locale.US);
		
		if (args.length > 0) {
			
			try {
			
				BeagleSequenceSimulatorConsoleApp app = new BeagleSequenceSimulatorConsoleApp();
				app.simulate(args);
				
			} catch (UnsupportedClassVersionError e) {

				Utils.handleException(e, "Your Java Runtime Environment is too old. Please update");

			}//END: try-catch block
				
		} else {

			try {

				Thread.setDefaultUncaughtExceptionHandler(new ExceptionHandler());
				new BeagleSequenceSimulatorApp();

				// sub-title for a software note
				System.out.println("Do the evolution baby!");

			} catch (UnsupportedClassVersionError e) {

				Utils.handleException(e, "Your Java Runtime Environment is too old. Please update");
				
			} catch (ClassNotFoundException e) {
				Utils.handleException(e, e.getMessage());
			} catch (InstantiationException e) {
				Utils.handleException(e, e.getMessage());
			} catch (IllegalAccessException e) {
				Utils.handleException(e, e.getMessage());
			} catch (UnsupportedLookAndFeelException e) {
				Utils.handleException(e, e.getMessage());
			}// END: try catch block

		}// END: command line check

	}// END: main

	private Image CreateImage(String path) {
		URL imgURL = this.getClass().getResource(path);
		Toolkit kit = Toolkit.getDefaultToolkit();
		Image img = kit.createImage(imgURL);

		if (img != null) {
			return img;
		} else {
			System.out.println("Couldn't find file: " + path + "\n");
			return null;
		}

	}// END: CreateImage

}// END: class

