/*
 * AuxilaryFrame.java
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

public class AuxilaryFrame extends AbstractFrame {

    private DocumentFrame documentFrame;
    private JPanel contentsPanel;

    public AuxilaryFrame(DocumentFrame documentFrame) {
        super();

        this.documentFrame = documentFrame;
        this.contentsPanel = null;
    }

    public AuxilaryFrame(DocumentFrame documentFrame,
                         JPanel contentsPanel) {
        super();

        this.documentFrame = documentFrame;
        setContentsPanel(contentsPanel);
    }

    public void setContentsPanel(JPanel contentsPanel) {
        this.contentsPanel = contentsPanel;
        getContentPane().add(contentsPanel);
        pack();
    }

    public DocumentFrame getDocumentFrame() {
        return documentFrame;
    }

    protected void initializeComponents() {
    }

    public boolean requestClose() {
        return true;
    }

    public JComponent getExportableComponent() {
        return contentsPanel;
    }

    public void doCloseWindow() {
        hide();
    }

    public Action getSaveAction() {
        return documentFrame.getSaveAction();
    }

    public Action getSaveAsAction() {
        return documentFrame.getSaveAsAction();
    }

}
