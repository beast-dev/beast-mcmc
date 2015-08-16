/*
 * ClockRateModelEditor.java
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
import javax.swing.JComboBox;
import javax.swing.JDialog;
import javax.swing.JFrame;
import javax.swing.JLabel;
import javax.swing.JPanel;
import javax.swing.SwingConstants;
import javax.swing.SwingUtilities;

import dr.app.gui.components.RealNumberField;

/**
 * @author Filip Bielejec
 * @version $Id$
 */
public class ClockRateModelEditor {

	// Data
	private PartitionDataList dataList = null;
	private int row;

	// Settings
	private OptionsPanel optionPanel;
	private JComboBox clockCombo;
	private RealNumberField[] clockParameterFields;
	
	// Buttons
	private JButton done;
	private JButton cancel;

	// Window
	private JDialog window;
	private Frame owner;

	public ClockRateModelEditor(PartitionDataList dataList, int row) {

		this.dataList = dataList;
		this.row = row;

		clockParameterFields = new RealNumberField[PartitionData.clockParameterNames.length];
		window = new JDialog(owner, "Setup clock rate model for partition "
				+ (row + 1));
		optionPanel = new OptionsPanel(12, 12, SwingConstants.CENTER);

		clockCombo = new JComboBox();
		clockCombo.setOpaque(false);

		for (String clockModel : PartitionData.clockModels) {
			clockCombo.addItem(clockModel);
		}// END: fill loop

		clockCombo.addItemListener(new ListenClockCombo());

		for (int i = 0; i < PartitionData.clockParameterNames.length; i++) {
			
			switch (i) {

			case 0:// clockrate
				clockParameterFields[i] = new RealNumberField(0.0, Double.MAX_VALUE);
				break;

			case 1: // ucld.mean
				clockParameterFields[i] = new RealNumberField(-Double.MAX_VALUE, Double.MAX_VALUE);
				break;

			case 2:// ucld.stdev
				clockParameterFields[i] = new RealNumberField(0.0, Double.MAX_VALUE);
				break;

			case 3: // ucld.offset
				clockParameterFields[i] = new RealNumberField(-Double.MAX_VALUE, Double.MAX_VALUE);
				break;

			case 4: // uced.mean
				clockParameterFields[i] = new RealNumberField(-Double.MAX_VALUE, Double.MAX_VALUE);
				break;

			case 5: // uced.offset
				clockParameterFields[i] = new RealNumberField(-Double.MAX_VALUE, Double.MAX_VALUE);
				break;

			case 6: // ig.mean
				clockParameterFields[i] = new RealNumberField(-Double.MAX_VALUE, Double.MAX_VALUE);
				break;

			case 7: // ig.stdev
				clockParameterFields[i] = new RealNumberField(0.0, Double.MAX_VALUE);
				break;

			case 8: // ig.offset
				clockParameterFields[i] = new RealNumberField(-Double.MAX_VALUE, Double.MAX_VALUE);
				break;

			default:
				clockParameterFields[i] = new RealNumberField();
			}// END: parameter switch
			
			clockParameterFields[i].setColumns(8);
			clockParameterFields[i]
					.setValue(dataList.get(row).clockParameterValues[i]);
		}// END: fill loop

		setClockArguments();

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
		window.getContentPane().add(optionPanel, BorderLayout.CENTER);
		window.getContentPane().add(buttonsHolder, BorderLayout.SOUTH);
		window.pack();
		
		// return to the previously chosen index on start
		clockCombo.setSelectedIndex(dataList.get(row).clockModelIndex);
		
	}// END: Constructor

	private void setClockArguments() {

		optionPanel.removeAll();
		optionPanel.addComponents(new JLabel("Clock rate model:"), clockCombo);
		optionPanel.addSeparator();
		optionPanel.addLabel("Set parameter values:");

		int index = clockCombo.getSelectedIndex();

		for (int i = 0; i < PartitionData.clockParameterIndices[index].length; i++) {

			int k = PartitionData.clockParameterIndices[index][i];

			if (k == 1 || k == 2) {
				clockParameterFields[k]
						.setToolTipText("Parameter in the real space");
			}//END: ucld.mean ucld.stdev check
			
			JPanel panel = new JPanel(new BorderLayout(6, 6));
			panel.add(clockParameterFields[k], BorderLayout.WEST);
			panel.setOpaque(false);
			optionPanel.addComponentWithLabel(
					PartitionData.clockParameterNames[k] + ":", panel);

		}// END: indices loop

		window.validate();
		window.repaint();
	}// END: setClockArguments

	private class ListenClockCombo implements ItemListener {
		public void itemStateChanged(ItemEvent ie) {

			setClockArguments();

		}// END: actionPerformed
	}// END: ListenClockCombo

	public void collectSettings() {

		dataList.get(row).clockModelIndex = clockCombo.getSelectedIndex();
		for (int i = 0; i < PartitionData.clockParameterNames.length; i++) {

			dataList.get(row).clockParameterValues[i] = clockParameterFields[i]
					.getValue();

		}// END: fill loop
	}// END: collectSettings

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
	
}// END: class
