/*
 * TreeAnnotatorDialog.java
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

package dr.app.beast;

import org.virion.jam.panels.OptionsPanel;
import org.virion.jam.components.WholeNumberField;
import org.virion.jam.components.RealNumberField;

import javax.swing.*;
import javax.swing.border.EmptyBorder;
import java.awt.*;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.event.ItemListener;
import java.awt.event.ItemEvent;
import java.io.File;

import dr.app.tools.TreeAnnotator;

public class BeastDialog {
	private JFrame frame;

	private OptionsPanel optionPanel;

    private WholeNumberField seedText = new WholeNumberField(1, Integer.MAX_VALUE);
    private JCheckBox beagleCheckBox = new JCheckBox("Use BEAGLE library");
    private JComboBox threadsCombo = new JComboBox(new Integer[] { 1, 2, 3, 4, 5, 6, 7, 8 });

	private File inputFile = null;

	public BeastDialog(final JFrame frame, final String titleString, final Icon icon) {
		this.frame = frame;

		optionPanel = new OptionsPanel(12, 12);

		this.frame = frame;

		JPanel panel = new JPanel(new BorderLayout());
		panel.setOpaque(false);

        final JLabel titleText = new JLabel(titleString);
        titleText.setIcon(icon);
        optionPanel.addSpanningComponent(titleText);
        titleText.setFont(new Font("sans-serif", 0, 12));
        
        final JButton inputFileButton = new JButton("Choose File...");
		final JTextField inputFileNameText = new JTextField("not selected", 16);

		inputFileButton.addActionListener( new ActionListener() {
			public void actionPerformed(ActionEvent ae) {
				FileDialog dialog = new FileDialog(frame,
						"Select target file...",
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

		JPanel panel1 = new JPanel(new BorderLayout(0,0));
		panel1.add(inputFileNameText, BorderLayout.CENTER);
		panel1.add(inputFileButton, BorderLayout.EAST);
		optionPanel.addComponentWithLabel("BEAST XML File: ", panel1);

        optionPanel.addSeparator();

        seedText.setColumns(12);
        optionPanel.addComponentWithLabel("Random number seed: ", seedText);

        optionPanel.addComponentWithLabel("Thread pool size: ", threadsCombo);

        optionPanel.addComponent(beagleCheckBox);
    

	}

	public boolean showDialog(String title, int seed) {

		JOptionPane optionPane = new JOptionPane(optionPanel,
				JOptionPane.PLAIN_MESSAGE,
				JOptionPane.OK_CANCEL_OPTION,
				null,
				new String[] { "Run", "Quit" },
				null);
		optionPane.setBorder(new EmptyBorder(12, 12, 12, 12));

        seedText.setValue(seed);

		final JDialog dialog = optionPane.createDialog(frame, title);
		//dialog.setResizable(true);
		dialog.pack();

		dialog.setVisible(true);

		return optionPane.getValue().equals("Run");
	}

    public int getSeed() {
        return seedText.getValue();
    }

    public boolean useBeagle() {
        return beagleCheckBox.isSelected();
    }

	public int getThreadPoolSize() {
		return (Integer)threadsCombo.getSelectedItem();
	}

    public File getInputFile() {
        return inputFile;
    }
}