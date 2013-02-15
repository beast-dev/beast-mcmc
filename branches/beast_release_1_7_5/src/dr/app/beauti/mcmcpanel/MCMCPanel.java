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
import dr.app.beauti.options.PartitionTreeModel;
import dr.app.beauti.options.STARBEASTOptions;
import dr.app.gui.components.WholeNumberField;
import dr.app.util.OSType;
import dr.evolution.datatype.Microsatellite;
import jam.panels.OptionsPanel;

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

    private static final long serialVersionUID = -3710586474593827540L;

    WholeNumberField chainLengthField = new WholeNumberField(1, Integer.MAX_VALUE);
    WholeNumberField echoEveryField = new WholeNumberField(1, Integer.MAX_VALUE);
    WholeNumberField logEveryField = new WholeNumberField(1, Integer.MAX_VALUE);

    JCheckBox samplePriorCheckBox = new JCheckBox("Sample from prior only - create empty alignment");

    public static final String DEFAULT_FILE_NAME_STEM = "untitled";
    JTextField fileNameStemField = new JTextField(DEFAULT_FILE_NAME_STEM);

    private JCheckBox addTxt = new JCheckBox("Add .txt suffix");

    JTextField logFileNameField = new JTextField(DEFAULT_FILE_NAME_STEM + ".log");
    JTextField treeFileNameField = new JTextField(DEFAULT_FILE_NAME_STEM + "." + STARBEASTOptions.TREE_FILE_NAME);
//    JCheckBox allowOverwriteLogCheck = new JCheckBox("Allow to overwrite the existing log file");

//    JCheckBox mapTreeLogCheck = new JCheckBox("Create tree file containing the MAP tree:");
//    JTextField mapTreeFileNameField = new JTextField("untitled.MAP.tree");

    JCheckBox substTreeLogCheck = new JCheckBox("Create tree log file with branch length in substitutions:");
    JTextField substTreeFileNameField = new JTextField("untitled(subst).trees");

    JCheckBox operatorAnalaysisCheck = new JCheckBox("Create operator analysis file:");
    JTextField operatorAnalaysisFileNameField = new JTextField(DEFAULT_FILE_NAME_STEM + ".ops");

    BeautiFrame frame = null;
    private final OptionsPanel optionsPanel;
    private BeautiOptions options;

    public MCMCPanel(BeautiFrame parent) {
        setLayout(new BorderLayout());

        // Mac OS X components have more spacing round them already
        optionsPanel = new OptionsPanel(12, (OSType.isMac() ? 6 : 24));

        this.frame = parent;

        setOpaque(false);
        optionsPanel.setOpaque(false);

        chainLengthField.setValue(100000);
        chainLengthField.setColumns(10);
        optionsPanel.addComponentWithLabel("Length of chain:", chainLengthField);
        chainLengthField.addKeyListener(new java.awt.event.KeyListener() {
            public void keyTyped(KeyEvent e) {
            }

            public void keyPressed(KeyEvent e) {
            }

            public void keyReleased(KeyEvent e) {
                options.chainLength = chainLengthField.getValue();
                frame.setDirty();
            }
        });

        optionsPanel.addSeparator();

        echoEveryField.setValue(1000);
        echoEveryField.setColumns(10);
        optionsPanel.addComponentWithLabel("Echo state to screen every:", echoEveryField);
        echoEveryField.addKeyListener(new java.awt.event.KeyListener() {
            public void keyTyped(KeyEvent e) {
            }

            public void keyPressed(KeyEvent e) {
            }

            public void keyReleased(KeyEvent e) {
                options.echoEvery = echoEveryField.getValue();
                frame.setDirty();
            }
        });

        logEveryField.setValue(100);
        logEveryField.setColumns(10);
        optionsPanel.addComponentWithLabel("Log parameters every:", logEveryField);
        logEveryField.addKeyListener(new java.awt.event.KeyListener() {
            public void keyTyped(KeyEvent e) {
            }

            public void keyPressed(KeyEvent e) {
            }

            public void keyReleased(KeyEvent e) {
                options.logEvery = logEveryField.getValue();
                frame.setDirty();
            }
        });

        optionsPanel.addSeparator();

        fileNameStemField.setColumns(32);
        optionsPanel.addComponentWithLabel("File name stem:", fileNameStemField);
        fileNameStemField.setEditable(true);
        fileNameStemField.addKeyListener(new java.awt.event.KeyListener() {
            public void keyTyped(KeyEvent e) {
            }

            public void keyPressed(KeyEvent e) {
            }

            public void keyReleased(KeyEvent e) {
                options.fileNameStem = fileNameStemField.getText();
                updateOtherFileNames(options);
                frame.setDirty();
            }
        });

        optionsPanel.addComponent(addTxt);
        if (OSType.isWindows()) {
            addTxt.setSelected(true);
        } else {
            addTxt.setSelected(false);
        }
        addTxt.addChangeListener(new ChangeListener() {
            public void stateChanged(ChangeEvent changeEvent) {
                setOptions(options);
                frame.setDirty();
            }
        });

        optionsPanel.addSeparator();

        logFileNameField.setColumns(32);
        optionsPanel.addComponentWithLabel("Log file name:", logFileNameField);
        logFileNameField.setEditable(false);

//        optionsPanel.addComponent(allowOverwriteLogCheck);
//        allowOverwriteLogCheck.setSelected(false);
//        allowOverwriteLogCheck.addChangeListener(new ChangeListener() {
//            public void stateChanged(ChangeEvent changeEvent) {
//            	options.allowOverwriteLog = allowOverwriteLogCheck.isSelected();
//            }
//        });

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

        optionsPanel.addComponent(operatorAnalaysisCheck);
        operatorAnalaysisCheck.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(ActionEvent e) {
                options.operatorAnalysis = operatorAnalaysisCheck.isSelected();

                updateOtherFileNames(options);

                frame.setDirty();
            }
        });

        operatorAnalaysisFileNameField.setColumns(32);
        operatorAnalaysisFileNameField.setEditable(false);
        operatorAnalaysisFileNameField.setEnabled(false);
        optionsPanel.addComponentWithLabel("Operator analysis file name:", operatorAnalaysisFileNameField);

        optionsPanel.addSeparator();

        optionsPanel.addComponent(samplePriorCheckBox);
        samplePriorCheckBox.setOpaque(false);
        samplePriorCheckBox.addChangeListener(new ChangeListener() {
            public void stateChanged(ChangeEvent changeEvent) {
                frame.setDirty();
            }
        });

//        logFileNameField.addKeyListener(listener);
//        treeFileNameField.addKeyListener(listener);
        //mapTreeFileNameField.addKeyListener(listener);
//        substTreeFileNameField.addKeyListener(listener);

//        optionsPanel.setPreferredSize(new java.awt.Dimension(500, 600));

        JPanel panel = new JPanel(new FlowLayout(FlowLayout.CENTER));
        panel.add(optionsPanel, BorderLayout.CENTER);
        panel.setOpaque(false);

        JScrollPane scrollPane = new JScrollPane(panel, JScrollPane.VERTICAL_SCROLLBAR_AS_NEEDED, JScrollPane.HORIZONTAL_SCROLLBAR_NEVER);
        scrollPane.setOpaque(false);
        scrollPane.setBorder(null);
        scrollPane.getViewport().setOpaque(false);

        add(scrollPane, BorderLayout.CENTER);
    }

    private void updateTreeFileNameList() {
        options.treeFileName.clear();
        options.substTreeFileName.clear();
        String treeFN;

        for (PartitionTreeModel tree : options.getPartitionTreeModels()) {
            if (options.substTreeLog) {
                treeFN = options.fileNameStem + "." + tree.getPrefix() + "(time)." + STARBEASTOptions.TREE_FILE_NAME;
            } else {
                treeFN = options.fileNameStem + "." + tree.getPrefix() + STARBEASTOptions.TREE_FILE_NAME; // stem.partitionName.tree
            }
            if (addTxt.isSelected()) treeFN = treeFN + ".txt";
            options.treeFileName.add(treeFN);

            if (options.substTreeLog) {
                treeFN = options.fileNameStem + "." + tree.getPrefix() + "(subst)." + STARBEASTOptions.TREE_FILE_NAME;
                if (addTxt.isSelected()) treeFN = treeFN + ".txt";
                options.substTreeFileName.add(treeFN);
            }
        }

        if (options.useStarBEAST) {
            treeFN = options.fileNameStem + "." + options.starBEASTOptions.SPECIES_TREE_FILE_NAME;
            if (addTxt.isSelected()) treeFN = treeFN + ".txt";
            options.treeFileName.add(treeFN);
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
        } else {
            fileNameStemField.setText(DEFAULT_FILE_NAME_STEM);
            fileNameStemField.setEnabled(false);
        }

        operatorAnalaysisCheck.setSelected(options.operatorAnalysis);

        updateOtherFileNames(options);

        if (options.contains(Microsatellite.INSTANCE)) {
            samplePriorCheckBox.setSelected(false);
            samplePriorCheckBox.setVisible(false);
        } else {
            samplePriorCheckBox.setVisible(true);
            samplePriorCheckBox.setSelected(options.samplePriorOnly);
        }

        optionsPanel.validate();
        optionsPanel.repaint();
    }

    private void updateOtherFileNames(BeautiOptions options) {
        if (options.fileNameStem != null) {
//            fileNameStemField.setText(options.fileNameStem);

            options.logFileName = options.fileNameStem + ".log";
            if (addTxt.isSelected()) options.logFileName = options.logFileName + ".txt";
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

            options.operatorAnalysisFileName = options.fileNameStem + ".ops";
            if (addTxt.isSelected()) {
                options.operatorAnalysisFileName = options.operatorAnalysisFileName + ".txt";
            }
            operatorAnalaysisFileNameField.setEnabled(options.operatorAnalysis);
            if (options.operatorAnalysis) {
                operatorAnalaysisFileNameField.setText(options.operatorAnalysisFileName);
            } else {
                operatorAnalaysisFileNameField.setText("");
            }

//            mapTreeLogCheck.setEnabled(true);
//            mapTreeLogCheck.setSelected(options.mapTreeLog);
//            mapTreeFileNameField.setEnabled(options.mapTreeLog);

            substTreeLogCheck.setEnabled(true);
            substTreeLogCheck.setSelected(options.substTreeLog);

        } else {
//            fileNameStemField.setText(fileNameStem);
//            fileNameStemField.setEnabled(false);
            logFileNameField.setText(DEFAULT_FILE_NAME_STEM + ".log");
            treeFileNameField.setText(DEFAULT_FILE_NAME_STEM + "." + STARBEASTOptions.TREE_FILE_NAME);
//            mapTreeLogCheck.setEnabled(false);
//            mapTreeFileNameField.setEnabled(false);
//            mapTreeFileNameField.setText("untitled");
            substTreeLogCheck.setSelected(false);
            substTreeFileNameField.setEnabled(false);
            substTreeFileNameField.setText("");
            operatorAnalaysisCheck.setSelected(false);
            operatorAnalaysisFileNameField.setText("");
        }
    }

    public void getOptions(BeautiOptions options) {
        options.fileNameStem = fileNameStemField.getText();
        options.logFileName = logFileNameField.getText();

//        options.mapTreeLog = mapTreeLogCheck.isSelected();
//        options.mapTreeFileName = mapTreeFileNameField.getText();

        options.substTreeLog = substTreeLogCheck.isSelected();
        updateTreeFileNameList();

        options.operatorAnalysis = operatorAnalaysisCheck.isSelected();
        options.operatorAnalysisFileName = operatorAnalaysisFileNameField.getText();

        options.samplePriorOnly = samplePriorCheckBox.isSelected();
    }

    public JComponent getExportableComponent() {
        return optionsPanel;
    }

}
