/*
 * ConsoleApplication.java
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

package org.virion.jam.console;

import org.virion.jam.framework.Application;
import org.virion.jam.framework.MenuBarFactory;

import javax.swing.*;
import java.io.File;
import java.io.IOException;

public class ConsoleApplication extends Application {

	private ConsoleFrame consoleFrame = null;

    public ConsoleApplication(String nameString, String aboutString, Icon icon) throws IOException {
        this(new ConsoleMenuBarFactory(), nameString, aboutString, icon);
    }

    public ConsoleApplication(MenuBarFactory menuBarFactory, String nameString, String aboutString, Icon icon) throws IOException {

		super(menuBarFactory, nameString, aboutString, icon);

		consoleFrame = new ConsoleFrame();
		consoleFrame.initialize();
		consoleFrame.setVisible(true);

		// event handling
		consoleFrame.addWindowListener(new java.awt.event.WindowAdapter() {
			public void windowClosing(java.awt.event.WindowEvent e) {
				thisWindowClosing(e);
			}
		});
	}

    public void initialize() {
        // nothing to do...
    }

    protected JFrame getDefaultFrame() { return consoleFrame; }

	public void doNew() {
		throw new RuntimeException("A ConsoleApplication cannot do a New command");
	}

	public void doOpenFile(File file) {
		throw new RuntimeException("A ConsoleApplication cannot do an Open command");
	}

	public void doCloseWindow() {
		doQuit();
	}

	public void doQuit() {
		if (consoleFrame.requestClose()) {

			consoleFrame.setVisible(false);
			consoleFrame.dispose();
			System.exit(0);
		}
	}

    public void doPreferences() {
    }

    public void doStop() {
		doQuit();
    }

	// Close the window when the close box is clicked
	private void thisWindowClosing(java.awt.event.WindowEvent e) {
		doQuit();
	}

}