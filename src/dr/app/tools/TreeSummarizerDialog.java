/*
 * TreeSummarizerDialog.java
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

package dr.app.tools;

import org.virion.jam.panels.OptionsPanel;
import org.virion.jam.components.WholeNumberField;

import javax.swing.*;
import javax.swing.border.EmptyBorder;
import java.io.File;
import java.awt.*;
import java.awt.event.ActionListener;
import java.awt.event.ActionEvent;

public class TreeSummarizerDialog {
    private JFrame frame;

    private OptionsPanel optionPanel;

    private WholeNumberField burninText = new WholeNumberField(0, Integer.MAX_VALUE);
    private JComboBox summaryTreeCombo = new JComboBox(new String[] { "Highest posterior frequency", "Maximum clade credibility" });

    private File inputFile = null;
    private File outputFile = null;

    public TreeSummarizerDialog(final JFrame frame) {
        this.frame = frame;

        optionPanel = new OptionsPanel(12, 12);

        this.frame = frame;

        JPanel panel = new JPanel(new BorderLayout());
        panel.setOpaque(false);

        burninText.setColumns(12);
        burninText.setValue(0);
        optionPanel.addComponentWithLabel("Burnin: ", burninText);
        optionPanel.addComponentWithLabel("Summary type: ", summaryTreeCombo);

        optionPanel.addSeparator();

        JButton inputFileButton = new JButton("Choose File...");
        final JTextField inputFileNameText = new JTextField("not selected", 16);

        inputFileButton.addActionListener( new ActionListener() {
            public void actionPerformed(ActionEvent ae) {
                FileDialog dialog = new FileDialog(frame,
                        "Select input tree file...",
                        FileDialog.LOAD);

                dialog.setVisible(true);
                if (dialog.getFile() == null) {
                    // the dialog was cancelled...
                    return;
                }

                inputFile = new File(dialog.getDirectory(), dialog.getFile());
                inputFileNameText.setText(inputFile.getName());

            }});
        inputFileNameText.setEditable(false);

        JPanel panel2 = new JPanel(new BorderLayout(0,0));
        panel2.add(inputFileNameText, BorderLayout.CENTER);
        panel2.add(inputFileButton, BorderLayout.EAST);
        optionPanel.addComponentWithLabel("Input Tree File: ", panel2);

        JButton outputFileButton = new JButton("Choose File...");
        final JTextField outputFileNameText = new JTextField("not selected", 16);

        outputFileButton.addActionListener( new ActionListener() {
            public void actionPerformed(ActionEvent ae) {
                FileDialog dialog = new FileDialog(frame,
                        "Select output file...",
                        FileDialog.SAVE);

                dialog.setVisible(true);
                if (dialog.getFile() == null) {
                    // the dialog was cancelled...
                    return;
                }

                outputFile = new File(dialog.getDirectory(), dialog.getFile());
                outputFileNameText.setText(outputFile.getName());

            }});
        outputFileNameText.setEditable(false);

        JPanel panel3 = new JPanel(new BorderLayout(0,0));
        panel3.add(outputFileNameText, BorderLayout.CENTER);
        panel3.add(outputFileButton, BorderLayout.EAST);
        optionPanel.addComponentWithLabel("Output File: ", panel3);
    }

    public boolean showDialog(String title) {

        JOptionPane optionPane = new JOptionPane(optionPanel,
                JOptionPane.PLAIN_MESSAGE,
                JOptionPane.OK_CANCEL_OPTION,
                null,
                new String[] { "Run", "Quit" },
                null);
        optionPane.setBorder(new EmptyBorder(12, 12, 12, 12));

        final JDialog dialog = optionPane.createDialog(frame, title);
        //dialog.setResizable(true);
        dialog.pack();

        dialog.setVisible(true);

        return optionPane.getValue().equals("Run");
    }

    public int getBurnin() {
        return burninText.getValue().intValue();
    }

    public String getInputFileName() {
        if (inputFile == null) return null;
        return inputFile.getPath();
    }

    public String getOutputFileName() {
        if (outputFile == null) return null;
        return outputFile.getPath();
    }

}
