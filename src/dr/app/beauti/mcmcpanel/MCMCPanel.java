/*
 * MCMCPanel.java
 *
 * Copyright (C) 2002-2009 Alexei Drummond and Andrew Rambaut
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
 * BEAST is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with BEAST; if not, write to the
 * Free Software Foundation, Inc., 51 Franklin St, Fifth Floor,
 * Boston, MA  02110-1301  USA
 */

package dr.app.beauti.mcmcpanel;

import dr.app.beauti.BeautiFrame;
import dr.app.beauti.BeautiPanel;
import dr.app.beauti.options.BeautiOptions;
import dr.app.beauti.options.PartitionSubstitutionModel;
import dr.app.beauti.options.PartitionTreeModel;
import dr.evomodel.coalescent.GMRFFixedGridImportanceSampler;
import org.virion.jam.components.WholeNumberField;
import org.virion.jam.panels.OptionsPanel;

import javax.swing.*;
import javax.swing.event.ChangeEvent;
import javax.swing.event.ChangeListener;
import java.awt.*;
import java.awt.event.ActionEvent;
import java.awt.event.KeyEvent;
import java.util.List;

/**
 * @author Andrew Rambaut
 * @author Alexei Drummond
 * @version $Id: MCMCPanel.java,v 1.16 2006/09/05 13:29:34 rambaut Exp $
 */
public class MCMCPanel extends BeautiPanel {

    /**
     *
     */
    private static final long serialVersionUID = -3710586474593827540L;
    WholeNumberField chainLengthField = new WholeNumberField(1, Integer.MAX_VALUE);
    WholeNumberField echoEveryField = new WholeNumberField(1, Integer.MAX_VALUE);
    WholeNumberField logEveryField = new WholeNumberField(1, Integer.MAX_VALUE);

    JCheckBox samplePriorCheckBox = new JCheckBox("Sample from prior only - create empty alignment");

    public static final String fileNameStem = "untitled";
    JTextField fileNameStemField = new JTextField(fileNameStem);

    JTextField logFileNameField = new JTextField(fileNameStem + ".log");
    JTextField treeFileNameField = new JTextField(fileNameStem + "." + GMRFFixedGridImportanceSampler.TREE_FILE_NAME);

    JCheckBox mapTreeLogCheck = new JCheckBox("Create tree file containing the MAP tree:");
    JTextField mapTreeFileNameField = new JTextField("untitled.MAP.tree");

    JCheckBox substTreeLogCheck = new JCheckBox("Create tree log file with branch length in substitutions:");
    JTextField substTreeFileNameField = new JTextField("untitled(subst).trees");

    BeautiFrame frame = null;
    private final OptionsPanel optionsPanel;
    private BeautiOptions options;

    public MCMCPanel(BeautiFrame parent) {
        setLayout(new BorderLayout());

        optionsPanel = new OptionsPanel(12, 24);

        this.frame = parent;

        setOpaque(false);
        optionsPanel.setOpaque(false);

        chainLengthField.setValue(100000);
        chainLengthField.setColumns(10);
        optionsPanel.addComponentWithLabel("Length of chain:", chainLengthField);

        optionsPanel.addSeparator();

        echoEveryField.setValue(1000);
        echoEveryField.setColumns(10);
        optionsPanel.addComponentWithLabel("Echo state to screen every:", echoEveryField);

        logEveryField.setValue(100);
        logEveryField.setColumns(10);
        optionsPanel.addComponentWithLabel("Log parameters every:", logEveryField);

        optionsPanel.addSeparator();

        fileNameStemField.setColumns(32);
        optionsPanel.addComponentWithLabel("File name stem:", fileNameStemField);
        fileNameStemField.setEditable(true);
        fileNameStemField.addKeyListener(new java.awt.event.KeyListener() {
            public void keyTyped(KeyEvent e) {
//				options.fileNameStem = fileNameStemField.getText();
//            	setOptions(options);
//                frame.setDirty();				
            }

            public void keyPressed(KeyEvent e) {
            }

            public void keyReleased(KeyEvent e) {
                options.fileNameStem = fileNameStemField.getText();
                setOptions(options);
                frame.setDirty();
            }
        });

        optionsPanel.addSeparator();

        logFileNameField.setColumns(32);
        optionsPanel.addComponentWithLabel("Log file name:", logFileNameField);
        logFileNameField.setEditable(false);
        treeFileNameField.setColumns(32);
        optionsPanel.addComponentWithLabel("Trees file name:", treeFileNameField);
        treeFileNameField.setEditable(false);

//        addComponent(mapTreeLogCheck);
//        mapTreeLogCheck.setOpaque(false);
//        mapTreeLogCheck.addActionListener(new java.awt.event.ActionListener() {
//            public void actionPerformed(ActionEvent e) {
//                mapTreeFileNameField.setEnabled(mapTreeLogCheck.isSelected());
//            }
//        });
//
//        mapTreeFileNameField.setColumns(32);
//        addComponentWithLabel("MAP tree file name:", mapTreeFileNameField);

        optionsPanel.addComponent(substTreeLogCheck);
        substTreeLogCheck.setOpaque(false);
        substTreeLogCheck.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(ActionEvent e) {
            	options.substTreeLog = substTreeLogCheck.isSelected();
            	updateTreeFileNameList();
            	substTreeFileNameField.setEnabled(substTreeLogCheck.isSelected());
                if (substTreeLogCheck.isSelected()) {
                	substTreeFileNameField.setText(displayTreeList(options.substTreeFileName));
                } else {
                    substTreeFileNameField.setText("");
                }

                frame.setDirty();
            }
        });

        substTreeFileNameField.setColumns(32);
        substTreeFileNameField.setEditable(false);
        substTreeFileNameField.setEnabled(false);
        optionsPanel.addComponentWithLabel("Substitutions trees file name:", substTreeFileNameField);

        java.awt.event.KeyListener listener = new java.awt.event.KeyAdapter() {
            public void keyTyped(java.awt.event.KeyEvent ev) {
                frame.setDirty();
            }
        };

        optionsPanel.addSeparator();

        optionsPanel.addComponent(samplePriorCheckBox);
        samplePriorCheckBox.setOpaque(false);
        samplePriorCheckBox.addChangeListener(new ChangeListener() {
            public void stateChanged(ChangeEvent changeEvent) {
                frame.setDirty();
            }
        });

        chainLengthField.addKeyListener(listener);
        echoEveryField.addKeyListener(listener);
        logEveryField.addKeyListener(listener);
        fileNameStemField.addKeyListener(listener);
//        logFileNameField.addKeyListener(listener);
//        treeFileNameField.addKeyListener(listener);
        //mapTreeFileNameField.addKeyListener(listener);
        substTreeFileNameField.addKeyListener(listener);

        add(optionsPanel, BorderLayout.CENTER);

    }

    private void updateTreeFileNameList(){
    	options.treeFileName.clear();
    	options.substTreeFileName.clear();
    	
    	for (PartitionTreeModel tree : options.getPartitionTreeModels()) {
            String treeFN;
            if (options.substTreeLog) {
                treeFN = options.fileNameStem + "." + tree.getPrefix() + "(time)." + GMRFFixedGridImportanceSampler.TREE_FILE_NAME;
            } else {
                treeFN = options.fileNameStem + "." + tree.getPrefix() + GMRFFixedGridImportanceSampler.TREE_FILE_NAME; // stem.partitionName.tree
            }
            options.treeFileName.add(treeFN);
            
            if (options.substTreeLog) {
            	options.substTreeFileName.add(options.fileNameStem + "." + tree.getPrefix() 
            			+ "(subst)." + GMRFFixedGridImportanceSampler.TREE_FILE_NAME);
            }            
        }
    	
    	if (options.isSpeciesAnalysis()) {
    		options.treeFileName.add(options.fileNameStem + "." + options.SPECIES_TREE_FILE_NAME);
    		//TODO: species sub tree
    	}
    }
    
    private String displayTreeList(List<String> treeList) {
    	String text = "";
    	
    	for (String t : treeList) {
    		text = text + t;
    		if (treeList.indexOf(t) < treeList.size() - 1) {
    			text = text + "; ";
    		}
    	}
    	
    	return text;
    }
    
    public void setOptions(BeautiOptions options) {
        this.options = options;

        chainLengthField.setValue(options.chainLength);

        echoEveryField.setValue(options.echoEvery);
        logEveryField.setValue(options.logEvery);

        if (options.fileNameStem != null) {
            fileNameStemField.setText(options.fileNameStem);

            options.logFileName = options.fileNameStem + ".log";
            logFileNameField.setText(options.logFileName);

//            if (options.mapTreeFileName == null) {
//			    mapTreeFileNameField.setText(options.fileNameStem + ".MAP.tree");
//            } else {
//                mapTreeFileNameField.setText(options.mapTreeFileName);
//            }
            
            updateTreeFileNameList();
            treeFileNameField.setText(displayTreeList(options.treeFileName));
            
            if (options.substTreeLog) {
                substTreeFileNameField.setText(displayTreeList(options.substTreeFileName));
            } else {
                substTreeFileNameField.setText("");
            }

//            mapTreeLogCheck.setEnabled(true);
//            mapTreeLogCheck.setSelected(options.mapTreeLog);
//            mapTreeFileNameField.setEnabled(options.mapTreeLog);

            substTreeLogCheck.setEnabled(true);
            substTreeLogCheck.setSelected(options.substTreeLog);
            
        } else {
            fileNameStemField.setText(fileNameStem);
            fileNameStemField.setEnabled(false);
            logFileNameField.setText(fileNameStem + ".log");
            treeFileNameField.setText(fileNameStem + "." + GMRFFixedGridImportanceSampler.TREE_FILE_NAME);
//            mapTreeLogCheck.setEnabled(false);
//            mapTreeFileNameField.setEnabled(false);
//            mapTreeFileNameField.setText("untitled");
            substTreeLogCheck.setEnabled(false);
            substTreeFileNameField.setEnabled(false);
            substTreeFileNameField.setText("untitled");
        }

        samplePriorCheckBox.setSelected(options.samplePriorOnly);

        optionsPanel.validate();
        optionsPanel.repaint();
    }

    public void getOptions(BeautiOptions options) {
        options.chainLength = chainLengthField.getValue();

        options.echoEvery = echoEveryField.getValue();
        options.logEvery = logEveryField.getValue();

        options.fileNameStem = fileNameStemField.getText();
        options.logFileName = logFileNameField.getText();        

//        options.mapTreeLog = mapTreeLogCheck.isSelected();
//        options.mapTreeFileName = mapTreeFileNameField.getText();

        options.substTreeLog = substTreeLogCheck.isSelected();
        updateTreeFileNameList();

        options.samplePriorOnly = samplePriorCheckBox.isSelected();
    }

    public JComponent getExportableComponent() {
        return optionsPanel;
    }

}
