/*
 * FrequencyModelEditor.java
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

package dr.app.bss;

import jam.panels.OptionsPanel;

import java.awt.BorderLayout;
import java.awt.Dimension;
import java.awt.Frame;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.event.ItemEvent;
import java.awt.event.ItemListener;

import javax.swing.JButton;
import javax.swing.JDialog;
import javax.swing.JFrame;
import javax.swing.JLabel;
import javax.swing.JPanel;
import javax.swing.JScrollPane;
import javax.swing.SwingConstants;
import javax.swing.SwingUtilities;

import dr.app.gui.components.RealNumberField;

/**
 * @author Filip Bielejec
 * @version $Id$
 */
public class FrequencyModelEditor {

	// Data
	private PartitionDataList dataList = null;
	private int row;
	
	// Settings
	private OptionsPanel optionPanel;
	private JScrollPane scrollPane;
	private DisabledItemsComboBox frequencyCombo;
	private RealNumberField[] frequencyParameterFields;
	
	//Buttons
	private JButton done;
	private JButton cancel;
	
	// Window
	private JDialog window;
	private Frame owner;
	
	public FrequencyModelEditor(PartitionDataList dataList, int row) {
		
		this.dataList = dataList;
		this.row = row;

		frequencyParameterFields = new RealNumberField[PartitionData.frequencyParameterNames.length];
		window = new JDialog(owner, "Setup base frequencies for partition " + (row + 1));
		optionPanel = new OptionsPanel(12, 12, SwingConstants.CENTER);
		
		scrollPane = new JScrollPane();
        optionPanel.setOpaque(false);
        scrollPane = new JScrollPane(optionPanel,
				JScrollPane.VERTICAL_SCROLLBAR_AS_NEEDED,
				JScrollPane.HORIZONTAL_SCROLLBAR_AS_NEEDED);
        scrollPane.getViewport().setOpaque(false);

		frequencyCombo = new DisabledItemsComboBox();
		
		int indexOf = 0;
		for (String frequencyModel : PartitionData.frequencyModels) {

			if (PartitionData.frequencyCompatibleDataTypes[indexOf] == dataList.get(row).dataTypeIndex) {

				frequencyCombo.addItem(frequencyModel, false);

			} else {
				
				frequencyCombo.addItem(frequencyModel, true);
				
			}// END: compatible check

			indexOf++;
		}// END: fill loop
		
		frequencyCombo.addItemListener(new ListenFrequencyCombo());

		for (int i = 0; i < PartitionData.frequencyParameterNames.length; i++) {
			frequencyParameterFields[i] = new RealNumberField(0.0, 1.0);
			frequencyParameterFields[i].setColumns(8);
			frequencyParameterFields[i].setValue(dataList.get(0).frequencyParameterValues[i]);
		}// END: fill loop

		setFrequencyArguments();
		
		// Buttons
		JPanel buttonsHolder = new JPanel();
		buttonsHolder.setOpaque(false);
		
		cancel = new JButton("Cancel", Utils.createImageIcon(Utils.CLOSE_ICON));
		cancel.addActionListener(new ListenCancel());
		buttonsHolder.add(cancel);
		
		done = new JButton("Done", Utils.createImageIcon(Utils.CHECK_ICON));
		done.addActionListener(new ListenOk());
		buttonsHolder.add(done);
		
		// Window
		owner = Utils.getActiveFrame();
		window.setLocationRelativeTo(owner);
		window.getContentPane().setLayout(new BorderLayout());
		window.getContentPane().add(scrollPane, BorderLayout.CENTER);
		window.getContentPane().add(buttonsHolder, BorderLayout.SOUTH);
		window.pack();
		
		//return to the previously chosen index on start
		frequencyCombo.setSelectedIndex(dataList.get(row).frequencyModelIndex);
		
	}//END: Constructor
	
	private void setFrequencyArguments() {

		optionPanel.removeAll();
		optionPanel.addComponents(new JLabel("Base frequencies:"), frequencyCombo);
		
		optionPanel.addSeparator();
		optionPanel.addComponentWithLabel("Set parameter values:", new JLabel());

		int index = frequencyCombo.getSelectedIndex();

		for (int i = 0; i < dataList.get(row).frequencyParameterIndices[index].length; i++) {

			int k = dataList.get(row).frequencyParameterIndices[index][i];

			JPanel panel = new JPanel(new BorderLayout(6, 6));
			panel.add(frequencyParameterFields[k], BorderLayout.WEST);
			panel.setOpaque(false);
			optionPanel.addComponentWithLabel(
					PartitionData.frequencyParameterNames[k] + ":",
					panel);

		}// END: indices loop
		
		window.validate();
		window.repaint();
	}// END: setFrequencyArguments

	private class ListenFrequencyCombo implements ItemListener {
		public void itemStateChanged(ItemEvent ie) {

			setFrequencyArguments();

		}// END: actionPerformed
	}// END: ListenFrequencyCombo

	public void collectSettings() {

		dataList.get(row).frequencyModelIndex = frequencyCombo.getSelectedIndex();
		for (int i = 0; i < PartitionData.frequencyParameterNames.length; i++) {

			dataList.get(row).frequencyParameterValues[i] = frequencyParameterFields[i].getValue();

		}// END: fill loop
	}// END: collectSettings
	
	//TODO: maybe check whether sum to 1?
	private class ListenOk implements ActionListener {
		public void actionPerformed(ActionEvent ev) {

			window.setVisible(false);
			collectSettings();
			
		}// END: actionPerformed
	}// END: ListenSaveLocationCoordinates
	
	private class ListenCancel implements ActionListener {
		public void actionPerformed(ActionEvent ev) {

			window.setVisible(false);
			
		}// END: actionPerformed
	}// END: ListenCancel
	
	public void showWindow() {
		window.setDefaultCloseOperation(JFrame.HIDE_ON_CLOSE);
		window.setSize(new Dimension(450, 400));
		window.setMinimumSize(new Dimension(100, 100));
		window.setResizable(true);
		window.setModal(true);
		window.setVisible(true);
	}// END: showWindow

	public void launch() {

		if (SwingUtilities.isEventDispatchThread()) {
			showWindow();
		} else {
			SwingUtilities.invokeLater(new Runnable() {
				public void run() {
					showWindow();
				}
			});
		}// END: edt check

	}// END: launch
	
}//END: class
