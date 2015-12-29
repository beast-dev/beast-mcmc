/*
 * InputsPanel.java
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

package dr.app.coalgen;

import jam.framework.Exportable;
import jam.panels.OptionsPanel;

import javax.swing.*;
import java.awt.*;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.io.File;
import java.io.IOException;

/**
 * @author Andrew Rambaut
 * @author Alexei Drummond
 * @version $Id$
 */
public class InputsPanel extends JPanel implements Exportable {

    private final CoalGenFrame frame;
    private final CoalGenData data;

    private OptionsPanel optionPanel;

    private final JButton logFileButton = new JButton("Choose File...");
    private final JTextField logFileNameText = new JTextField("not selected", 16);

    private final JButton treesFileButton = new JButton("Choose File...");
    private final JTextField treesFileNameText = new JTextField("not selected", 16);


    public InputsPanel(final CoalGenFrame frame, final CoalGenData data) {

        super();

        this.frame = frame;
        this.data = data;

        setOpaque(false);
        setLayout(new BorderLayout());

        optionPanel = new OptionsPanel(12, 12, SwingConstants.CENTER);
        add(optionPanel, BorderLayout.NORTH);

        logFileButton.addActionListener(new ActionListener() {
            public void actionPerformed(ActionEvent ae) {
                FileDialog dialog = new FileDialog(frame,
                        "Select input log file...",
                        FileDialog.LOAD);

                dialog.setVisible(true);
                if (dialog.getFile() == null) {
                    // the dialog was cancelled...
                    return;
                }

                logFileButton.setEnabled(false);

                File file = new File(dialog.getDirectory(), dialog.getFile());
                try {
                    frame.readFromFile(file);
                } catch (IOException e) {
                    e.printStackTrace();
                }
            }
        });
        logFileNameText.setEditable(false);

        treesFileButton.addActionListener(new ActionListener() {
            public void actionPerformed(ActionEvent ae) {
                FileDialog dialog = new FileDialog(frame,
                        "Select input trees file...",
                        FileDialog.LOAD);

                dialog.setVisible(true);
                if (dialog.getFile() == null) {
                    // the dialog was cancelled...
                    return;
                }

                data.treesFile = new File(dialog.getDirectory(), dialog.getFile());
                treesFileNameText.setText(data.treesFile.getName());

                frame.fireTracesChanged();
            }
        });
        treesFileNameText.setEditable(false);

        JPanel panel1 = new JPanel(new BorderLayout(0, 0));
        panel1.setOpaque(false);
        panel1.add(logFileNameText, BorderLayout.CENTER);
        panel1.add(logFileButton, BorderLayout.EAST);
        optionPanel.addComponentWithLabel("Log File: ", panel1);

        JPanel panel2 = new JPanel(new BorderLayout(0, 0));
        panel2.setOpaque(false);
        panel2.add(treesFileNameText, BorderLayout.CENTER);
        panel2.add(treesFileButton, BorderLayout.EAST);
        optionPanel.addComponentWithLabel("Input Trees File: ", panel2);
    }

    public final void tracesChanged() {
        logFileButton.setEnabled(true);
        if (data.logFile != null) {
            logFileNameText.setText(data.logFile.getName());
        }
    }

    public void collectSettings() {

    }

    public JComponent getExportableComponent() {
        return this;
    }
}