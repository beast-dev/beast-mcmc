/*
 * SimulationPanel.java
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

import dr.app.gui.components.WholeNumberField;
import dr.evolution.alignment.SimpleAlignment;
import jam.framework.Exportable;
import jam.panels.OptionsPanel;

import javax.swing.*;

import java.awt.*;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.event.ItemEvent;
import java.awt.event.ItemListener;

/**
 * @author Filip Bielejec
 * @version $Id$
 */
@SuppressWarnings({ "rawtypes", "serial" })
public class SimulationPanel extends JPanel implements Exportable {

    private MainFrame frame;
    private PartitionDataList dataList;
    private OptionsPanel optionPanel;
    
    private WholeNumberField simulationsNumberField;
    private WholeNumberField startingSeedNumberField;

    // Buttons
    private JButton simulate;
    private JButton generateXML;

    // Check boxes
    private JCheckBox setSeed;
    private JCheckBox useParallel;
    private JCheckBox outputAncestralSequences;
    
    //Combo boxes
    private JComboBox outputFormat;
	private ComboBoxModel outputFormatModel;
    
    @SuppressWarnings({ "unchecked" })
	public SimulationPanel(final MainFrame frame,
                           final PartitionDataList dataList) {

        this.frame = frame;
        this.dataList = dataList;

        optionPanel = new OptionsPanel(12, 12, SwingConstants.CENTER);

        simulationsNumberField = new WholeNumberField(1, Integer.MAX_VALUE);
        simulationsNumberField.setColumns(10);
        simulationsNumberField.setValue(dataList.simulationsCount);
        optionPanel.addComponentWithLabel("Number of simulations:",
                simulationsNumberField);

        setSeed = new JCheckBox();
        setSeed.addItemListener(new SetSeedCheckBoxListener());
        setSeed.setSelected(dataList.setSeed);
        optionPanel.addComponentWithLabel("Set seed:", setSeed);

        startingSeedNumberField = new WholeNumberField(1, Long.MAX_VALUE);
        startingSeedNumberField.setColumns(10);
        startingSeedNumberField.setValue(dataList.startingSeed);
        startingSeedNumberField.setEnabled(dataList.setSeed);
        optionPanel.addComponentWithLabel("Starting seed:",
                startingSeedNumberField);

    	outputFormat = new JComboBox();
        optionPanel.addComponentWithLabel("Output format:", outputFormat);
        outputFormatModel = new DefaultComboBoxModel(SimpleAlignment.OutputType.values());
    	outputFormat.setModel(outputFormatModel);

        outputAncestralSequences = new JCheckBox();
        outputAncestralSequences.addItemListener(new outputAncestralSequencesCheckBoxListener());
        outputAncestralSequences.setSelected(dataList.useParallel);
        optionPanel.addComponentWithLabel("Output ancestral sequences:",
        		outputAncestralSequences);
    	
        useParallel = new JCheckBox();
        useParallel.addItemListener(new UseParallelCheckBoxListener());
        useParallel.setSelected(dataList.useParallel);
        optionPanel.addComponentWithLabel("Use parallel implementation:",
                useParallel);

        // Buttons holder
        JPanel buttonsHolder = new JPanel();
        buttonsHolder.setOpaque(false);

        // simulate button
        simulate = new JButton("Simulate",
                Utils.createImageIcon(Utils.BIOHAZARD_ICON));
        simulate.addActionListener(new ListenSimulate());
        buttonsHolder.add(simulate);

        generateXML = new JButton("Generate XML",
                Utils.createImageIcon(Utils.HAMMER_ICON));
        generateXML.addActionListener(new ListenGenerateXML());
        buttonsHolder.add(generateXML);

        setOpaque(false);
        setLayout(new BorderLayout());
        add(optionPanel, BorderLayout.NORTH);
        add(buttonsHolder, BorderLayout.SOUTH);

    }// END: SimulationPanel

    public final void collectSettings() {

        dataList.simulationsCount = simulationsNumberField.getValue();
        if (dataList.setSeed) {
            dataList.startingSeed = startingSeedNumberField.getValue();
        }

        dataList.outputFormat = SimpleAlignment.OutputType.parseFromString(
                outputFormat.getSelectedItem().toString());

    }// END: collectSettings

    private class SetSeedCheckBoxListener implements ItemListener {
        public void itemStateChanged(ItemEvent e) {

            if (setSeed.isSelected()) {
                startingSeedNumberField.setEnabled(true);
                dataList.setSeed = true;
            } else {
                startingSeedNumberField.setEnabled(false);
                dataList.setSeed = false;
            }

        }
    }// END: CheckBoxListener

    private class UseParallelCheckBoxListener implements ItemListener {
        public void itemStateChanged(ItemEvent e) {

            if (useParallel.isSelected()) {
                dataList.useParallel = true;
            } else {
                dataList.useParallel = false;
            }

        }
    }// END: CheckBoxListener

    private class outputAncestralSequencesCheckBoxListener implements ItemListener {
        public void itemStateChanged(ItemEvent e) {

            if (outputAncestralSequences.isSelected()) {
                dataList.outputAncestralSequences = true;
            } else {
                dataList.outputAncestralSequences = false;
            }

        }
    }// END: CheckBoxListener
    
    private class ListenSimulate implements ActionListener {
        public void actionPerformed(ActionEvent ev) {

            frame.doExport();

        }// END: actionPerformed
    }// END: ListenSimulate

    private class ListenGenerateXML implements ActionListener {
        public void actionPerformed(ActionEvent ev) {

            frame.doGenerateXML();

        }// END: actionPerformed
    }// END: ListenSaveLocationCoordinates

    public void setBusy() {
        simulate.setEnabled(false);
        generateXML.setEnabled(false);
    }// END: setBusy

    public void setIdle() {
        simulate.setEnabled(true);
        generateXML.setEnabled(true);
    }// END: setIdle

    public JComponent getExportableComponent() {
        return this;
    }// END: getExportableComponent

    public void updateSimulationPanel(PartitionDataList dataList) {
    	
    	setDataList(dataList);
    	
		// TODO: DOES NOT WORK
    	outputFormatModel.setSelectedItem(dataList.outputFormat);
    	outputFormat.setSelectedItem(dataList.outputFormat);
    	
    }//END: updateSimulationPanel
    
    public void setDataList(PartitionDataList dataList) {
        this.dataList = dataList;
    }// END: setDataList

}// END: class
