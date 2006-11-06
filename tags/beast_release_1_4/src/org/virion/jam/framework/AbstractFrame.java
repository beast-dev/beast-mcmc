/*
 * AbstractFrame.java
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

import org.virion.jam.util.PrintUtilities;

import javax.swing.*;
import java.awt.event.ActionEvent;
import java.awt.print.Printable;
import java.awt.print.PrinterException;
import java.awt.print.PrinterJob;
import java.awt.*;

public abstract class AbstractFrame extends JFrame implements Exportable {

    private JMenuBar menuBar = null;
    private boolean isDirty = false;

    public AbstractFrame() {
    }

    public final void initialize() {
        initializeComponents();
        if (menuBar == null) {
            menuBar = new JMenuBar();
            if(Application.getMenuBarFactory() == null)
                return;
            Application.getMenuBarFactory().populateMenuBar(menuBar, this);
        }
        setJMenuBar(menuBar);
    }

    protected abstract void initializeComponents();

    public final boolean isDirty() {
        return isDirty;
    }

    public final void setDirty() {
        getRootPane().putClientProperty("windowModified", Boolean.TRUE);
        this.isDirty = true;
    }

    public final void clearDirty() {
        getRootPane().putClientProperty("windowModified", Boolean.FALSE);
        this.isDirty = false;
    }

    public abstract boolean requestClose();

    public void doImport() {
        throw new RuntimeException("Not implemented in AbstractFrame - this must be overridden");
    }

    public void doExport() {
        throw new RuntimeException("Not implemented in AbstractFrame - this must be overridden");
    }

    public final void doPrint() {
        doPrint(false);
    }
    public final void doPrint(final boolean scaleIfDoesntImplementPrintable) {
        // create a separate thread to do this, since it locks up the CPU
        // for about five seconds, and it looks incredibly ugly when the main window
        // is only partially painted.
/*        Runnable runnable = new Runnable() {
            public void run() {*/
                final PrinterJob printJob = PrinterJob.getPrinterJob();

        JComponent component = getExportableComponent();
        if (component != null) {
            if (component instanceof Printable) {

                printJob.setPrintable((Printable) component);
                if (printJob.printDialog()) {
//                            RepaintManager.currentManager(component).paintDirtyRegions();
                    try {
                        printJob.print();
                    } catch (PrinterException pe) {
                                JOptionPane.showMessageDialog(AbstractFrame.this, "Printing error: " + pe,
                                "Error Printing",
                                JOptionPane.ERROR_MESSAGE);
                    }
                }
            } else {
                        if (scaleIfDoesntImplementPrintable)
                            PrintUtilities.printComponentScaled(component);
                        else
                PrintUtilities.printComponent(component);
            }
        } else {
                    JOptionPane.showMessageDialog(AbstractFrame.this, "Printing error: No panel provided to print",
                    "Error Printing",
                    JOptionPane.ERROR_MESSAGE);
        }
/*            }
        };
        Thread thread = new Thread(runnable);
        thread.setPriority(Thread.MIN_PRIORITY);
        thread.start();*/

    }

    public void doCloseWindow() {
        if (requestClose()) {
            setVisible(false);
            dispose();
        }
    }

    public void doZoomWindow() {
        // not supported by JDK1.3
        if (Toolkit.getDefaultToolkit().isFrameStateSupported(JFrame.MAXIMIZED_BOTH)) {
        	this.setExtendedState(JFrame.MAXIMIZED_BOTH);
        }
    }

    public void doMinimizeWindow() {
        // not supported by JDK1.3
        if (Toolkit.getDefaultToolkit().isFrameStateSupported(JFrame.ICONIFIED)) {
        	this.setExtendedState(JFrame.ICONIFIED);
        }
    }

    public void doCut() {
        throw new RuntimeException("Not implemented in AbstractFrame - this must be overridden");
    }

    public void doCopy() {
        throw new RuntimeException("Not implemented in AbstractFrame - this must be overridden");
    }

    public void doPaste() {
        throw new RuntimeException("Not implemented in AbstractFrame - this must be overridden");
    }

    public void doDelete() {
        throw new RuntimeException("Not implemented in AbstractFrame - this must be overridden");
    }

    public void doSelectAll() {
        throw new RuntimeException("Not implemented in AbstractFrame - this must be overridden");
    }

    public void doFind() {
        throw new RuntimeException("Not implemented in AbstractFrame - this must be overridden");
    }

    public Action getNewAction() {
        return Application.getApplication().getNewAction();
    }

    public Action getOpenAction() {
        return Application.getApplication().getOpenAction();
    }

    public Action getPageSetupAction() {
        return Application.getApplication().getPageSetupAction();
    }

    public Action getExitAction() {
        return Application.getApplication().getExitAction();
    }

    public Action getCloseWindowAction() {
        return closeWindowAction;
    }

    public Action getAboutAction() {
        return Application.getApplication().getAboutAction();
    }

    public Action getPreferencesAction() {
        return Application.getApplication().getPreferencesAction();
    }

    public Action getWebsiteAction() {
        return Application.getApplication().getWebsiteAction();
    }

    public Action getSaveAction() {
        return saveAction;
    }

    public Action getSaveAsAction() {
        return saveAsAction;
    }


    public Action getZoomWindowAction() {
        return zoomWindowAction;
    }

    public Action getMinimizeWindowAction() {
        return minimizeWindowAction;
    }

	public void setImportAction(Action importAction) {
		assert importAction != null: "Import Action already set for this frame";
	    this.importAction = importAction;
	}

    public void setExportAction(Action exportAction) {
	    assert exportAction != null: "Export Action already set for this frame";
        this.exportAction = exportAction;
    }

    public Action getImportAction() {
        return importAction;
    }

    public Action getExportAction() {
        return exportAction;
    }

    public Action getPrintAction() {
        return printAction;
    }

    public Action getCutAction() {
        return cutAction;
    }

    public Action getCopyAction() {
        return copyAction;
    }

    public Action getPasteAction() {
        return pasteAction;
    }

    public Action getDeleteAction() {
        return deleteAction;
    }

    public Action getSelectAllAction() {
        return selectAllAction;
    }

    public Action getFindAction() {
        return findAction;
    }

    /**
     * override this to provide a document specific help menu item
     */
    public Action getHelpAction() {
        return null;
    }

    private AbstractAction saveAction = new AbstractAction("Save") {
        /**
		 * 
		 */
		private static final long serialVersionUID = -4540478663138470260L;

		public void actionPerformed(ActionEvent ae) {
            // Do nothing.. This is just a dummy action - getSaveAction should be overriden to provide a proper action
        }
    };

    private AbstractAction saveAsAction = new AbstractAction("Save As...") {
        /**
		 * 
		 */
		private static final long serialVersionUID = -1584673587355509882L;

		public void actionPerformed(ActionEvent ae) {
            // Do nothing.. This is just a dummy action - getSaveAction should be overriden to provide a proper action
        }
    };

    private Action importAction = null;
    private Action exportAction = null;

    private AbstractAction printAction = new AbstractAction("Print...") {
        /**
		 * 
		 */
		private static final long serialVersionUID = 7608963909000084240L;

		public void actionPerformed(ActionEvent ae) {
            doPrint();
        }
    };

    protected AbstractAction closeWindowAction = new AbstractAction("Close") {
        /**
		 * 
		 */
		private static final long serialVersionUID = -8653835282951538865L;

		public void actionPerformed(ActionEvent ae) {
            doCloseWindow();
        }
    };

    private AbstractAction zoomWindowAction = new AbstractAction("Zoom Window") {
        /**
		 * 
		 */
		private static final long serialVersionUID = -5303318185278668551L;

		public void actionPerformed(ActionEvent ae) {
            doZoomWindow();
        }
    };

    private AbstractAction minimizeWindowAction = new AbstractAction("Minimize Window") {
        /**
		 * 
		 */
		private static final long serialVersionUID = 4122140938271674271L;

		public void actionPerformed(ActionEvent ae) {
            doMinimizeWindow();
        }
    };

    private AbstractAction cutAction = new AbstractAction("Cut") {
        /**
		 * 
		 */
		private static final long serialVersionUID = -389065906967196702L;

		public void actionPerformed(ActionEvent ae) {
            doCut();
        }
    };

    private AbstractAction copyAction = new AbstractAction("Copy") {
        /**
		 * 
		 */
		private static final long serialVersionUID = 1150131467885279315L;

		public void actionPerformed(ActionEvent ae) {
            doCopy();
        }
    };

    private AbstractAction pasteAction = new AbstractAction("Paste") {
        /**
		 * 
		 */
		private static final long serialVersionUID = 7399738765663916935L;

		public void actionPerformed(ActionEvent ae) {
            doPaste();
        }
    };

    private AbstractAction deleteAction = new AbstractAction("Delete") {
        /**
		 * 
		 */
		private static final long serialVersionUID = 59368389446966562L;

		public void actionPerformed(ActionEvent ae) {
            doDelete();
        }
    };

    private AbstractAction selectAllAction = new AbstractAction("Select All") {
        /**
		 * 
		 */
		private static final long serialVersionUID = 378627540775330105L;

		public void actionPerformed(ActionEvent ae) {
            doSelectAll();
        }
    };

    private AbstractAction findAction = new AbstractAction("Find...") {
        /**
		 * 
		 */
		private static final long serialVersionUID = 1693178283344457760L;

		public void actionPerformed(ActionEvent ae) {
            doFind();
        }
    };
}
