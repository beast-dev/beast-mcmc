/*
 * SingleDocApplication.java
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

package org.virion.jam.framework;

import javax.swing.*;
import java.io.File;

public class SingleDocApplication extends Application {

    private DocumentFrame documentFrame = null;

    public SingleDocApplication(String nameString, String aboutString, Icon icon) {

        super(new SingleDocMenuBarFactory(), nameString, aboutString, icon);
    }

    public SingleDocApplication(String nameString, String aboutString, Icon icon,
    							String websiteURLString, String helpURLString) {

        super(new SingleDocMenuBarFactory(), nameString, aboutString, icon, websiteURLString, helpURLString);
    }

    public SingleDocApplication(MenuBarFactory menuBarFactory, String nameString, String aboutString, Icon icon) {

        super(menuBarFactory, nameString, aboutString, icon);
    }

    public SingleDocApplication(MenuBarFactory menuBarFactory, String nameString, String aboutString, Icon icon,
    							String websiteURLString, String helpURLString) {

        super(menuBarFactory, nameString, aboutString, icon, websiteURLString, helpURLString);
    }

	public final void initialize() {
		// nothing to do...
	}

    public void setDocumentFrame(DocumentFrame documentFrame) {

        this.documentFrame = documentFrame;

        documentFrame.initialize();
        documentFrame.setVisible(true);

        // event handling
        documentFrame.addWindowListener(new java.awt.event.WindowAdapter() {
            public void windowClosing(java.awt.event.WindowEvent e) {
                thisWindowClosing(e);
            }
        });
    }

	protected JFrame getDefaultFrame() { return documentFrame; }

	protected String getDocumentExtension() { return ""; }

    public void doNew() {
        throw new RuntimeException("A SingleDocApplication cannot do a New command");
    }

    public void doOpenFile(File file) {
        documentFrame.openFile(file);
    }

    public void doCloseWindow() {
        doQuit();
    }

    public void doQuit() {
        if (documentFrame == null) {
            return;
        }
        if (documentFrame.requestClose()) {

            documentFrame.setVisible(false);
            documentFrame.dispose();
            System.exit(0);
        }
    }

    public void doPreferences() {
    }

    // Close the window when the close box is clicked
    private void thisWindowClosing(java.awt.event.WindowEvent e) {
        doQuit();
    }
}