/*
 * Application.java
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

import org.virion.jam.html.HTMLViewer;
import org.virion.jam.util.BrowserLauncher;

import javax.swing.*;
import java.awt.*;
import java.awt.event.ActionEvent;
import java.awt.print.PageFormat;
import java.awt.print.PrinterJob;
import java.io.*;

/*
 * @todo Implement a list of open windows
 * @todo Implement the recent files menu (persistance)
 * @todo Implement a preferences system
 */

public abstract class Application {

    private static MenuBarFactory menuBarFactory;
    private static Icon icon;
    private static String nameString;
    private static String aboutString;
    private static String websiteURLString;
    private static String helpURLString;

    private static Application application = null;

    private JMenu recentFileMenu = null;

    public static Application getApplication() {
        return application;
    }

    public static MenuBarFactory getMenuBarFactory() {
        return menuBarFactory;
    }

    public static Icon getIcon() {
        return icon;
    }

    public static String getNameString() {
        return nameString;
    }

    public static String getAboutString() {
        return aboutString;
    }

    public static String getWebsiteURLString() {
        return websiteURLString;
    }

    public static String getHelpURLString() {
        return helpURLString;
    }

    public Application(MenuBarFactory menuBarFactory, String nameString, String aboutString, Icon icon) {
    	this(menuBarFactory, nameString, aboutString, icon, null, null);
    }

    public Application(MenuBarFactory menuBarFactory, String nameString, String aboutString, Icon icon,
    					String websiteURLString, String helpURLString) {

        Application.menuBarFactory = menuBarFactory;
        Application.nameString = nameString;
        Application.aboutString = aboutString;
        Application.websiteURLString = websiteURLString;
        Application.helpURLString = helpURLString;
        Application.icon = icon;

        aboutAction = new AbstractAction("About " + nameString + "...") {
            /**
			 * 
			 */
			private static final long serialVersionUID = -5041909266347767945L;

			public void actionPerformed(ActionEvent ae) {
                doAbout();
            }
        };

        if (application != null) {
            throw new RuntimeException("Only on instance of Application is allowed");
        }
        application = this;

        if (org.virion.jam.mac.Utils.isMacOSX()) {
            org.virion.jam.mac.Utils.macOSXRegistration(application);
        }
    }

	public abstract void initialize();

	public void addMenuFactory(MenuFactory menuFactory) {
		getMenuBarFactory().registerMenuFactory(menuFactory);
	}

    public JMenu getRecentFileMenu() {
        if (recentFileMenu == null) {
            recentFileMenu = new JMenu("Recent Files");
        }
	    recentFileMenu.setEnabled(getOpenAction().isEnabled());
        return recentFileMenu;
    }

    public void addRecentFile(File file) {

        if (recentFileMenu != null) {
            if (recentFileMenu.getItemCount() == 20) {
                recentFileMenu.remove(19);
            }
            recentFileMenu.insert(file.toString(), 0);

            //menuBar.validate();
        }
    }

    protected abstract JFrame getDefaultFrame();

    public void doAbout() {
        AboutBox aboutBox = new AboutBox(getNameString(), getAboutString(), getIcon());
        //aboutBox.initialize();        //causes about frame to have the menu system from the main frame.
        aboutBox.setVisible(true);
    }

    public void doHelp() {
    	if (helpURLString != null) {
    		displayURL(helpURLString);
    	} else {
	        try {
	            InputStream in = getClass().getResourceAsStream("/help/application.help");
	            if (in == null) return;
	            Reader reader = new InputStreamReader(in);
	            StringWriter writer = new StringWriter();
	            int c;
	            while ((c = reader.read()) != -1) writer.write(c);
	            reader.close();
	            writer.close();
	            JFrame frame = new HTMLViewer(getNameString() + " Help", writer.toString());
	            frame.setVisible(true);
	        } catch (IOException ignore) {
	        }
		}
    }

    public void doWebsite() {
    	if (websiteURLString != null) {
    		displayURL(websiteURLString);
    	}
    }

    public void displayURL(String urlString) {
        try {
            BrowserLauncher.openURL(urlString);
        } catch (IOException ioe) {
            // do nothing
        }
    }

    public void doPageSetup() {
        PrinterJob printJob = PrinterJob.getPrinterJob();
        printJob.pageDialog(new PageFormat());
    }

    public void doOpen() {
        Frame frame = getDefaultFrame();
        if (frame == null) {
            frame = new JFrame();
        }

        FileDialog dialog = new FileDialog(frame,
                "Open Document",
                FileDialog.LOAD);
        dialog.setVisible(true);
        if (dialog.getFile() != null) {
            File file = new File(dialog.getDirectory(), dialog.getFile());
            doOpenFile(file);
            addRecentFile(file);
        }
    }

	public void doOpen(String fileName) {
        if (fileName != null && fileName.length() > 0) {
            File file = new File(fileName);
            doOpenFile(file);
            addRecentFile(file);
        }
	}

    public abstract void doNew();

    public abstract void doOpenFile(File file);

    public abstract void doQuit();

    public void doPreferences() {
    }

    public Action getNewAction() {
        return newAction;
    }

    public Action getOpenAction() {
        return openAction;
    }

    public Action getPageSetupAction() {
        return pageSetupAction;
    }

    public Action getExitAction() {
        return exitAction;
    }

    public Action getAboutAction() {
        return aboutAction;
    }

    public Action getPreferencesAction() {
        return preferencesAction;
    }

    public Action getHelpAction() {
        return helpAction;
    }

    public Action getWebsiteAction() {
        return websiteAction;
    }

    protected AbstractAction newAction = new AbstractAction("New") {
        /**
		 * 
		 */
		private static final long serialVersionUID = -1521224399481196083L;

		public void actionPerformed(ActionEvent ae) {
            doNew();
        }
    };

    protected AbstractAction openAction = new AbstractAction("Open...") {
        /**
		 * 
		 */
		private static final long serialVersionUID = -4447260642111669429L;

		public void actionPerformed(ActionEvent ae) {
            doOpen();
        }
    };

    protected AbstractAction pageSetupAction = new AbstractAction("Page Setup...") {
        /**
		 * 
		 */
		private static final long serialVersionUID = 2066709411506589740L;

		public void actionPerformed(ActionEvent ae) {
            doPageSetup();
        }
    };

    protected AbstractAction exitAction = new AbstractAction("Exit") {
        /**
		 * 
		 */
		private static final long serialVersionUID = -4248102898807999282L;

		public void actionPerformed(ActionEvent ae) {
            doQuit();
        }
    };

    protected AbstractAction aboutAction = null;

    protected AbstractAction preferencesAction = new AbstractAction("Preferences...") {
        /**
		 * 
		 */
		private static final long serialVersionUID = -6992111486368469010L;

		public void actionPerformed(ActionEvent ae) {
            doPreferences();
        }
    };

    protected AbstractAction helpAction = new AbstractAction("Help...") {
        /**
		 * 
		 */
		private static final long serialVersionUID = -3959121576921987360L;

		public void actionPerformed(ActionEvent ae) {
            doHelp();
        }
    };

    protected AbstractAction websiteAction = new AbstractAction("Website...") {
        /**
		 * 
		 */
		private static final long serialVersionUID = 837521123495914197L;

		public void actionPerformed(ActionEvent ae) {
            doWebsite();
        }
    };

}
