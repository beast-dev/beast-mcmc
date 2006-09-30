/*
 * Utils.java
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

package org.virion.jam.app;

import java.awt.*;
import java.io.File;
import java.io.PrintWriter;
import java.io.StringWriter;

/**
 *
 * @author Alexei Drummond
 * @author Andrew Rambaut
 *
 * @version $Id: Utils.java,v 1.1.1.1 2006/07/16 13:17:56 rambaut Exp $
 */
public class Utils {

	
	private static javax.swing.JFileChooser SAVE_FILE_CHOOSER = null;
        
	public static String getLoadFileName(String message) {
		java.io.File file = getLoadFile(message);
		if (file == null) return null;
		return file.getAbsolutePath();
	}

	public static String getSaveFileName(String message) {
		java.io.File file = getSaveFile(message);
		if (file == null) return null;
		return file.getAbsolutePath();
	}

	public static File getLoadFile(String message) {
		// No file name in the arguments so throw up a dialog box...
		java.awt.Frame frame = new java.awt.Frame();
		java.awt.FileDialog chooser = new java.awt.FileDialog(frame, message, 
															java.awt.FileDialog.LOAD);
		chooser.show();
		if (chooser.getFile() == null) return null;
		java.io.File file = new java.io.File(chooser.getDirectory(), chooser.getFile());
		chooser.dispose();
		frame.dispose();
		
		return file;
	}

	public static File getSaveFile(String message) {
		// No file name in the arguments so throw up a dialog box...
		java.awt.Frame frame = new java.awt.Frame();
		java.awt.FileDialog chooser = new java.awt.FileDialog(frame, message,
															java.awt.FileDialog.SAVE);
        chooser.setVisible(true);
        java.io.File file = null;
        if(chooser.getDirectory() != null && chooser.getFile() != null)
            file = new java.io.File(chooser.getDirectory(), chooser.getFile());

        chooser.dispose();
		frame.dispose();

		return file;
	}

	/**
	 * This function takes a file name and an array of extensions (specified
	 * without the leading '.'). If the file name ends with one of the extensions
	 * then it is returned with this trimmed off. Otherwise the file name is
	 * return as it is.
	 * @return the trimmed filename
	 */
	public static String trimExtensions(String fileName, String[] extensions) {

		String newName = null;
		
		for (int i = 0; i < extensions.length; i++) {
			String ext = "." + extensions[i];
			if (fileName.endsWith(ext)) {
				newName = fileName.substring(0, fileName.length() - ext.length());
			}
		}
		
		if (newName == null) newName = fileName;
	
		return newName;
	}

	/**
	 * @return a named image from file or resource bundle.
	 */
	public static Image getImage(Object caller, String name) {
		
		java.net.URL url = caller.getClass().getResource(name);
		if (url != null) { 
			return Toolkit.getDefaultToolkit().createImage(url);
		} else {
			if (caller instanceof Component) {
				Component c = (Component)caller;
				Image i = c.createImage(100,20);
				Graphics g = c.getGraphics();
				g.drawString("Not found!", 1, 15);
				return i;
			} else return null;
		}
	}
       
	
	/*public static void showExceptionDialog(JFrame parent, Throwable e, String title) {
		
		JOptionPane pane = new JOptionPane(e.getMessage(), JOptionPane.ERROR_MESSAGE);
		JDialog dialog = pane.createDialog(parent, title);
		dialog.show();
	}

	public static void showExceptionDialog(JFrame parent, Throwable e) {
		
		JExceptionDialog d = new JExceptionDialog(parent, true, "Exception", getStackTrace(e));	
	}*/

	/**
	 * Return the stack trace of an exception as a string.
	 */
	private static String getStackTrace(Throwable e) {
		StringWriter sw = new StringWriter();
		PrintWriter pw = new PrintWriter(sw);
		e.printStackTrace(pw);
		return sw.toString();
	}

}
