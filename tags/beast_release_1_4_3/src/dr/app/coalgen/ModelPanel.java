/*
 * ModelPanel.java
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

package dr.app.coalgen;

import org.virion.jam.panels.OptionsPanel;
import org.virion.jam.components.RealNumberField;
import org.virion.jam.components.WholeNumberField;
import org.virion.jam.framework.Exportable;

import javax.swing.*;
import java.awt.*;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.event.ItemEvent;
import java.awt.event.ItemListener;

/**
 * @author Andrew Rambaut
 * @author Alexei Drummond
 *
 * @version $Id$
 */
public class ModelPanel extends OptionsPanel implements Exportable {
	
	CoalGenFrame frame = null;
	CoalGenData data = null;
	
	private JComboBox demographicCombo; 

	private OptionsPanel optionPanel;

    private WholeNumberField replicatesField = new WholeNumberField(1, Integer.MAX_VALUE);

	private RealNumberField[] argumentFields = new RealNumberField[CoalGenData.argumentNames.length];
	private JCheckBox[] argumentCheckBoxes = new JCheckBox[CoalGenData.argumentNames.length]; 
	private JComboBox[] argumentCombos = new JComboBox[CoalGenData.argumentNames.length]; 
	

	public ModelPanel(CoalGenFrame frame, CoalGenData data) {
	
		super();
		
		this.frame = frame;
		this.data = data;
		
		setOpaque(false);
		setLayout(new BorderLayout());

        replicatesField.setColumns(8);

		demographicCombo = new JComboBox();
		demographicCombo.setOpaque(false);
		
		for (int i = 0; i < CoalGenData.demographicModels.length; i++) {
			demographicCombo.addItem(CoalGenData.demographicModels[i]);
		}
				
		demographicCombo.addItemListener(new ItemListener() {
			public void itemStateChanged(ItemEvent e) {
				setDemographicArguments();
			}});

		for (int i = 0; i < CoalGenData.argumentNames.length; i++) {
			argumentFields[i] = new RealNumberField();
			argumentFields[i].setColumns(8);
			argumentCheckBoxes[i] = new JCheckBox("From Trace:");
			argumentCheckBoxes[i].setOpaque(false);
			argumentCombos[i] = new JComboBox();
			argumentCombos[i].setEnabled(false);
			argumentCombos[i].setOpaque(false);
			argumentCheckBoxes[i].addActionListener(
				new ArgumentActionListener(argumentCheckBoxes[i], argumentFields[i], argumentCombos[i]));
		}
		
	  	optionPanel = new OptionsPanel();
		add(optionPanel, BorderLayout.CENTER);
		
		setDemographicArguments();
	}
	
	class ArgumentActionListener implements ActionListener {
		JCheckBox checkBox;
		RealNumberField field;
		JComboBox combo;
		
		ArgumentActionListener(JCheckBox checkBox, RealNumberField field, JComboBox combo) {
			this.checkBox = checkBox;
			this.field = field;
			this.combo = combo;
		}

		public void actionPerformed(ActionEvent ae) {
			field.setEnabled(!checkBox.isSelected());
			combo.setEnabled(checkBox.isSelected());
		}
	};
			
	private int findArgument(JComboBox comboBox, String argument) {
		for (int i = 0; i < comboBox.getItemCount(); i++) {
			String item = ((String)comboBox.getItemAt(i)).toLowerCase();
			if (item.indexOf(argument) != -1) return i;
		}
		return -1;
	}

	private void setDemographicArguments() {
 		optionPanel.removeAll();
 		
        optionPanel.addComponentWithLabel("Number of replicates (ignored if using a trace file):", replicatesField);

        optionPanel.addSeperator();

		optionPanel.addComponents(new JLabel("Demographic Model:"), demographicCombo);

        optionPanel.addSeperator();

 		optionPanel.addLabel("Select the parameter values (or obtain from a trace file):");

 		int demo = demographicCombo.getSelectedIndex();
 		
 		for (int i = 0; i < data.argumentIndices[demo].length; i++) {
 			int k = data.argumentIndices[demo][i];
 			
 			JPanel panel = new JPanel(new BorderLayout(6,6));
 			panel.add(argumentFields[k], BorderLayout.WEST);
 			panel.add(argumentCheckBoxes[k], BorderLayout.CENTER);
 			panel.add(argumentCombos[k], BorderLayout.EAST);
			panel.setOpaque(false);
			optionPanel.addComponentWithLabel(data.argumentNames[k] + ":", 
													panel);
		}
        validate();
        repaint();
	}
			
	public final void dataChanged() {
        replicatesField.setValue(data.replicateCount);
		demographicCombo.setSelectedIndex(data.demographicModel);

		if (data.traces != null) {
			for (int i = 0; i < CoalGenData.argumentNames.length; i++) {
				argumentCombos[i].removeAllItems();
				for (int j = 0; j < data.traces.getTraceCount(); j++) {
					String statistic = data.traces.getTraceName(j);
					argumentCombos[i].addItem(statistic);
				}
			
				int index = data.argumentTraces[i];
				
				for (int j = 0; j < data.argumentGuesses[i].length; j++) {
					if (index != -1) break;
					
					index = findArgument(argumentCombos[i], data.argumentGuesses[i][j]);
				}
				if (index == -1) {
					argumentCheckBoxes[i].setEnabled(false);
					index = 0;
				} else {
					argumentCheckBoxes[i].setEnabled(true);
				}
				
				argumentCombos[i].setSelectedIndex(index);
			}
		}
        for (int i = 0; i < CoalGenData.argumentNames.length; i++) {
            argumentFields[i].setValue(data.argumentValues[i]);
        }
        setDemographicArguments();
	}

	public final void updateData() {
        data.replicateCount = replicatesField.getValue().intValue();

		data.demographicModel = demographicCombo.getSelectedIndex();
		for (int i = 0; i < CoalGenData.argumentNames.length; i++) {
			data.argumentValues[i] = argumentFields[i].getValue().doubleValue();
			if (argumentCheckBoxes[i].isSelected()) {
				data.argumentTraces[i] = argumentCombos[i].getSelectedIndex();
			} else {
				data.argumentTraces[i] = -1;
			}
		}
	}

    public JComponent getExportableComponent() {
		return this;
	} 	
}
