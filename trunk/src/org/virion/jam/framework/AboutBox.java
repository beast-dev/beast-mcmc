/*
 * AboutBox.java
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

import org.virion.jam.util.IconUtils;
import org.virion.jam.util.Utils;

import javax.swing.*;
import java.awt.*;
import java.awt.event.*;
import java.util.StringTokenizer;

public class AboutBox extends AbstractFrame {

    /**
     * Creates an AboutBox with a given title, message and icon
     * and centers it over the parent component.
     */
    public AboutBox(String title, String message, Icon icon) {
        super();

		if (icon != null) {
	        setIconImage(IconUtils.getBufferedImageFromIcon(icon));
		}

        addWindowListener(new WindowAdapter() {
            public void windowClosing(WindowEvent event) {
                close();
            }
        });
        addKeyListener(new KeyAdapter() {
            public void keyPressed(KeyEvent event) {
                if (event.getKeyCode() == KeyEvent.VK_ESCAPE
                        ||
                        event.getKeyCode() == KeyEvent.VK_ENTER) {
                    close();
                }
            }
        });

	    JPanel contentsPanel = new JPanel(new GridBagLayout());
        contentsPanel.setBorder(BorderFactory.createEmptyBorder(10, 10, 10, 10));
	    contentsPanel.setBackground(Color.white);

	    JLabel iconLabel = new JLabel(icon, JLabel.CENTER);
	    JLabel titleLabel = new JLabel(title, JLabel.CENTER);

        Font font = titleLabel.getFont();
        titleLabel.setFont(font.deriveFont(16.0f).deriveFont(Font.BOLD));

	    GridBagConstraints c = new GridBagConstraints();
	    c.gridwidth = GridBagConstraints.REMAINDER;
	    c.insets = new Insets(5,5,5,5);
	    contentsPanel.add(iconLabel, c);
	    contentsPanel.add(titleLabel, c);

        font = font.deriveFont(10.0f);
        StringTokenizer tokens = new StringTokenizer(message, "\n");
        while (tokens.hasMoreElements()) {
            String text = tokens.nextToken();
            JLabel messageLabel = new JLabel(text, JLabel.CENTER);
            messageLabel.setFont(font);
	        contentsPanel.add(messageLabel, c);
        }

        getSaveAction().setEnabled(false);
        getSaveAsAction().setEnabled(false);
	    getPrintAction().setEnabled(false);
	    getPageSetupAction().setEnabled(false);

        getCutAction().setEnabled(false);
        getCopyAction().setEnabled(false);
        getPasteAction().setEnabled(false);
        getDeleteAction().setEnabled(false);
        getSelectAllAction().setEnabled(false);
        getFindAction().setEnabled(false);

        getZoomWindowAction().setEnabled(false);
	    getMinimizeWindowAction().setEnabled(true);
	    getCloseWindowAction().setEnabled(true);

        getContentPane().add(contentsPanel);
        setDefaultCloseOperation(DISPOSE_ON_CLOSE);
        setResizable(false);
        pack();
        Utils.centerComponent(this, null);
	}

    /**
     * Sets the visibility to false and disposes the frame.
     */
    public void close() {
        setVisible(false);
        dispose();
    }

    protected void initializeComponents() {
    }

    public boolean requestClose() {
        return false;
    }

    public JComponent getExportableComponent() {
        return null;
    }
}


