/*
 * TreeAnnotatorDialog.java
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

package dr.app.tools;

import dr.app.gui.components.RealNumberField;
import dr.app.gui.FileDrop;
import jam.panels.OptionsPanel;
import dr.app.gui.components.WholeNumberField;

import javax.swing.*;
import javax.swing.border.Border;
import javax.swing.border.EmptyBorder;
import java.awt.*;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.event.ItemListener;
import java.awt.event.ItemEvent;
import java.io.File;

public class TreeAnnotatorDialog {
	private JFrame frame;

	private OptionsPanel optionPanel;

    private WholeNumberField burninStatesText = new WholeNumberField(0, Long.MAX_VALUE);
    private WholeNumberField burninTreesText = new WholeNumberField(0, Long.MAX_VALUE);
	private RealNumberField limitText = new RealNumberField(0.0, 1.0);

    private JComboBox summaryTreeCombo = new JComboBox(TreeAnnotator.Target.values());
    private JComboBox nodeHeightsCombo = new JComboBox(TreeAnnotator.HeightsSummary.values());

	private File targetFile = null;
	private File inputFile = null;
	private File outputFile = null;

	public TreeAnnotatorDialog(final JFrame frame) {
		this.frame = frame;

		optionPanel = new OptionsPanel(12, 12);

		this.frame = frame;

		JPanel panel = new JPanel(new BorderLayout());
		panel.setOpaque(false);

        final JRadioButton burninStatesRadio = new JRadioButton("Specify the burnin as the number of states");
        optionPanel.addSpanningComponent(burninStatesRadio);
        burninStatesText.setColumns(12);
        burninStatesText.setValue(0);
        burninStatesText.setToolTipText("<html>Specify the burnin as the number of states<br> in the MCMC chain</html>");
        final JLabel burninStatesLabel = optionPanel.addComponentWithLabel("Burnin (as states): ", burninStatesText);

        final JRadioButton burninTreesRadio = new JRadioButton("Specify the burnin as the number of trees");
        optionPanel.addSpanningComponent(burninTreesRadio);
        burninTreesText.setColumns(12);
        burninTreesText.setValue(0);
        burninTreesText.setToolTipText("<html>Specify the burnin as the number of trees<br> in the tree file</html>");
        final JLabel burninTreesLabel = optionPanel.addComponentWithLabel("Burnin (as trees): ", burninTreesText);

        ButtonGroup group = new ButtonGroup();
        group.add(burninStatesRadio);
        group.add(burninTreesRadio);

        //Register a listener for the radio buttons.
        burninStatesRadio.addItemListener(new ItemListener() {
            @Override
            public void itemStateChanged(ItemEvent e) {
                burninStatesLabel.setEnabled(burninStatesRadio.isSelected());
                burninStatesText.setEnabled(burninStatesRadio.isSelected());

                burninTreesLabel.setEnabled(!burninStatesRadio.isSelected());
                burninTreesText.setEnabled(!burninStatesRadio.isSelected());
            }
        });
        burninStatesRadio.setSelected(true);

        limitText.setColumns(12);
		limitText.setValue(0.0);
        limitText.setToolTipText("<html>Specify a lower limit on the posterior probability<br>" +
                "below which a clade will not be annotated</html>");
        optionPanel.addComponentWithLabel("Posterior probability limit: ", limitText);

        optionPanel.addComponentWithLabel("Target tree type: ", summaryTreeCombo);
        optionPanel.addComponentWithLabel("Node heights: ", nodeHeightsCombo);

        optionPanel.addSeparator();

        final JButton targetFileButton = new JButton("Choose File...");
		final JTextField targetFileNameText = new JTextField("not selected", 16);

		targetFileButton.addActionListener( new ActionListener() {
			public void actionPerformed(ActionEvent ae) {
				FileDialog dialog = new FileDialog(frame,
						"Select target file...",
						FileDialog.LOAD);

				dialog.setVisible(true);
				if (dialog.getFile() == null) {
					// the dialog was cancelled...
					return;
				}

				targetFile = new File(dialog.getDirectory(), dialog.getFile());
				targetFileNameText.setText(targetFile.getName());

			}});
		targetFileNameText.setEditable(false);

		JPanel panel1 = new JPanel(new BorderLayout(0,0));
		panel1.add(targetFileNameText, BorderLayout.CENTER);
		panel1.add(targetFileButton, BorderLayout.EAST);
		final JLabel label1 = optionPanel.addComponentWithLabel("Target Tree File: ", panel1);

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

        label1.setEnabled(false);
        targetFileNameText.setEnabled(false);
        targetFileButton.setEnabled(false);

        summaryTreeCombo.addItemListener(new ItemListener() {
            public void itemStateChanged(ItemEvent itemEvent) {
                boolean selected = summaryTreeCombo.getSelectedItem().equals(TreeAnnotator.Target.USER_TARGET_TREE);
                label1.setEnabled(selected);
                targetFileNameText.setEnabled(selected);
                targetFileButton.setEnabled(selected);
            }
        });

        JPanel panel2 = new JPanel(new BorderLayout(0,0));
		panel2.add(inputFileNameText, BorderLayout.CENTER);
		panel2.add(inputFileButton, BorderLayout.EAST);

        Color focusColor = UIManager.getColor("Focus.color");
        Border focusBorder = BorderFactory.createMatteBorder( 2, 2, 2, 2, focusColor );
        new FileDrop( null, inputFileNameText, focusBorder, new FileDrop.Listener()
        {   public void filesDropped( java.io.File[] files )
            {
                inputFile = files[0];
                inputFileNameText.setText(inputFile.getName());
            }   // end filesDropped
        }); // end FileDrop.Listener

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

    public int getBurninStates() {
        return burninStatesText.getValue();
    }

    public int getBurninTrees() {
        return burninTreesText.getValue();
    }

    public double getPosteriorLimit() {
		return limitText.getValue();
	}

    public TreeAnnotator.Target getTargetOption() {
        return (TreeAnnotator.Target)summaryTreeCombo.getSelectedItem();
    }

    public TreeAnnotator.HeightsSummary getHeightsOption() {
        return (TreeAnnotator.HeightsSummary)nodeHeightsCombo.getSelectedItem();
    }

    public String getTargetFileName() {
		if (targetFile == null) return null;
		return targetFile.getPath();
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
