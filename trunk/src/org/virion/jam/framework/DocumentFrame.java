/*
 * DocumentFrame.java
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
import java.awt.*;
import java.awt.event.ActionEvent;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;

public abstract class DocumentFrame extends AbstractFrame {

    private File documentFile = null;

    public DocumentFrame() {
        super();
    }

    protected abstract void initializeComponents();

    protected abstract boolean readFromFile(File file) throws FileNotFoundException, IOException;

    protected abstract boolean writeToFile(File file) throws IOException;

    public final boolean hasFile() {
        return documentFile != null;
    }

    public final File getFile() {
        return documentFile;
    }

    public boolean requestClose() {
        if (isDirty()) {
            int option = JOptionPane.showConfirmDialog(this, "Do you wish to save?",
                    "Unsaved changes",
                    JOptionPane.YES_NO_CANCEL_OPTION,
                    JOptionPane.WARNING_MESSAGE);

            if (option == JOptionPane.YES_OPTION) {
                return doSave();
            } else if (option == JOptionPane.CANCEL_OPTION || option == -1) {
                return false;
            }
            return true;
        }
        return true;
    }

    public void openFile(File file) {

        try {
            if (readFromFile(file)) {

                clearDirty();
                documentFile = file;
            }
        } catch (FileNotFoundException fnfe) {
            JOptionPane.showMessageDialog(this, "Unable to open file: File not found",
                    "Unable to open file",
                    JOptionPane.ERROR_MESSAGE);
        } catch (IOException ioe) {
            JOptionPane.showMessageDialog(this, "Unable to read file: " + ioe,
                    "Unable to read file",
                    JOptionPane.ERROR_MESSAGE);
        }
    }

    public final boolean doSave() {
        if (!hasFile()) {
            return doSaveAs();
        } else {
            try {
                if (writeToFile(documentFile)) {

                    clearDirty();
                }
            } catch (IOException ioe) {
                JOptionPane.showMessageDialog(this, "Unable to save file: " + ioe,
                        "Unable to save file",
                        JOptionPane.ERROR_MESSAGE);
            }
        }
        return true;
    }

    public final boolean doSaveAs() {
        FileDialog dialog = new FileDialog(this,
                "Save Document As...",
                FileDialog.SAVE);

        dialog.setVisible(true);
        if (dialog.getFile() == null) {
            // the dialog was cancelled...
            return false;
        }

        File file = new File(dialog.getDirectory(), dialog.getFile());

        try {
            if (writeToFile(file)) {

                clearDirty();
                documentFile = file;
            }
        } catch (IOException ioe) {
            JOptionPane.showMessageDialog(this, "Unable to save file: " + ioe,
                    "Unable to save file",
                    JOptionPane.ERROR_MESSAGE);
        }

        return true;
    }

    public Action getSaveAction() {
        return saveAction;
    }

    public Action getSaveAsAction() {
        return saveAsAction;
    }

    private AbstractAction saveAction = new AbstractAction("Save") {
        /**
		 * 
		 */
		private static final long serialVersionUID = -7222084794535198182L;

		public void actionPerformed(ActionEvent ae) {
            doSave();
        }
    };

    private AbstractAction saveAsAction = new AbstractAction("Save As...") {
        /**
		 * 
		 */
		private static final long serialVersionUID = 201002809081635295L;

		public void actionPerformed(ActionEvent ae) {
            doSaveAs();
        }
    };

}
