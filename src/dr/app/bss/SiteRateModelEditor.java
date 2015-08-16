/*
 * SiteRateModelEditor.java
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
import javax.swing.JSpinner;
import javax.swing.SpinnerModel;
import javax.swing.SpinnerNumberModel;
import javax.swing.SwingConstants;
import javax.swing.SwingUtilities;
import javax.swing.text.BadLocationException;

import dr.app.gui.components.RealNumberField;

/**
 * @author Filip Bielejec
 * @version $Id$
 */
public class SiteRateModelEditor {

	// Data
	private PartitionDataList dataList = null;
	private int row;
	
	// Settings	
	private OptionsPanel optionPanel;
	private JComboBox siteCombo;
	private RealNumberField[] siteParameterFields;
    private JSpinner gammaCategoriesSpinner;
	
	//Buttons
	private JButton done;
	private JButton cancel;
	
	// Window
	private JDialog window;
	private Frame owner;
    
	public SiteRateModelEditor(PartitionDataList dataList, int row) throws NumberFormatException, BadLocationException {

		this.dataList = dataList;
		this.row = row;
		
		siteParameterFields = new RealNumberField[PartitionData.siteRateModelParameterNames.length];
		window = new JDialog(owner, "Setup site rate model for partition " + (row + 1));
		optionPanel = new OptionsPanel(12, 12, SwingConstants.CENTER);

		siteCombo = new JComboBox();
		siteCombo.setOpaque(false);

		for (String siteModel : PartitionData.siteRateModels) {
			siteCombo.addItem(siteModel);
		}// END: fill loop

		siteCombo.addItemListener(new ListenSiteCombo());

		for (int i = 0; i < PartitionData.siteRateModelParameterNames.length; i++) {
			
			switch (i) {

			case 0: // GammaCategories
				siteParameterFields[i] = new RealNumberField(1.0, Double.valueOf(Integer.MAX_VALUE));
				break;

			case 1: // Alpha
				siteParameterFields[i] = new RealNumberField(0.0, Double.MAX_VALUE);
				break;

			case 2: // Invariant sites proportion
				siteParameterFields[i] = new RealNumberField(0.0, 1.0);
				break;

			default:
				siteParameterFields[i] = new RealNumberField();

			}//END: parameter switch
			
			siteParameterFields[i].setColumns(8);
			siteParameterFields[i].setValue(dataList.get(0).siteRateModelParameterValues[i]);
		}// END: fill loop

		setSiteArguments();

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
		
		//return to the previously chosen index on start
		siteCombo.setSelectedIndex(dataList.get(row).siteRateModelIndex);
		
	}// END: Constructor

	private void setSiteArguments() throws NumberFormatException, BadLocationException {

		optionPanel.removeAll();
		optionPanel.addComponents(new JLabel("Site Rate model:"), siteCombo);
		optionPanel.addSeparator();
		optionPanel.addLabel("Set parameter values:");

		int index = siteCombo.getSelectedIndex();
		
		for (int i = 0; i < PartitionData.siteRateModelParameterIndices[index].length; i++) {

			if(index == 1 && i == 0) {
				
				int k = PartitionData.siteRateModelParameterIndices[index][i];
				
				Integer initValue = Integer.valueOf(siteParameterFields[k].getText(0, 1)); 
				Integer	min = 1;
				Integer max = 10;//Integer.MAX_VALUE;
				Integer step = 1;
				
				SpinnerModel model = new SpinnerNumberModel(initValue, min, max, step);
				gammaCategoriesSpinner = new JSpinner(model);
				
				JPanel panel = new JPanel(new BorderLayout(6, 6));
				panel.add(gammaCategoriesSpinner, BorderLayout.WEST);
				panel.setOpaque(false);
				optionPanel.addComponentWithLabel(
						PartitionData.siteRateModelParameterNames[k] + ":",
						panel);
				
			} else {
			
			int k = PartitionData.siteRateModelParameterIndices[index][i];

			JPanel panel = new JPanel(new BorderLayout(6, 6));
			panel.add(siteParameterFields[k], BorderLayout.WEST);
			panel.setOpaque(false);
			optionPanel.addComponentWithLabel(
					PartitionData.siteRateModelParameterNames[k] + ":",
					panel);

			}// END: gama categories field check
			
		}// END: indices loop
		
		window.validate();
		window.repaint();
	}// END: setSiteArguments

	private class ListenSiteCombo implements ItemListener {
		public void itemStateChanged(ItemEvent ie) {

			try {

				setSiteArguments();

			} catch (NumberFormatException e) {
				e.printStackTrace();
			} catch (BadLocationException e) {
				e.printStackTrace();
			}

		}// END: actionPerformed
	}// END: ListenSiteCombo

	public void collectSettings() {

		int index = siteCombo.getSelectedIndex();
		dataList.get(row).siteRateModelIndex = index;
		
		for (int i = 0; i < PartitionData.siteRateModelParameterNames.length; i++) {

			if(index == 1 && i == 0) { 
				
				dataList.get(row).siteRateModelParameterValues[i] = Double.valueOf(gammaCategoriesSpinner.getValue().toString()); 
						
			} else {
			
				dataList.get(0).siteRateModelParameterValues[i] = siteParameterFields[i].getValue();
			
			}// END: gama categories field check

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
