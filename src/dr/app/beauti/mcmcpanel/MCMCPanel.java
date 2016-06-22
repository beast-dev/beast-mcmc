/*
 * MCMCPanel.java
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

package dr.app.beauti.mcmcpanel;


import dr.app.beagle.tools.Partition;
import dr.app.beauti.BeautiFrame;
import dr.app.beauti.BeautiPanel;
import dr.app.beauti.components.marginalLikelihoodEstimation.MLEDialog;
import dr.app.beauti.components.marginalLikelihoodEstimation.MLEGSSDialog;
import dr.app.beauti.components.marginalLikelihoodEstimation.MarginalLikelihoodEstimationOptions;
import dr.app.beauti.options.AbstractPartitionData;
import dr.app.beauti.options.BeautiOptions;
import dr.app.beauti.options.PartitionTreeModel;
import dr.app.beauti.options.STARBEASTOptions;
import dr.app.beauti.util.PanelUtils;
import dr.app.gui.components.WholeNumberField;
import dr.app.util.OSType;
import dr.evolution.datatype.DataType;
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
    JComboBox performMLECombo = new JComboBox(new String[] {"None", "path sampling/stepping-stone sampling", "generalized stepping-stone sampling"});
    //    JCheckBox performMLE = new JCheckBox("Perform marginal likelihood estimation (MLE) using path sampling/stepping-stone sampling");
    JButton buttonMLE = new JButton("Settings");
    //    JCheckBox performMLEGSS = new JCheckBox("Perform marginal likelihood estimation (MLE) using generalized stepping-stone sampling");
    //  JButton buttonMLEGSS = new JButton("Settings");

    public static final String DEFAULT_FILE_NAME_STEM = "untitled";
    JTextField fileNameStemField = new JTextField(DEFAULT_FILE_NAME_STEM);

    private JCheckBox addTxt = new JCheckBox("Add .txt suffix");

    JTextArea logFileNameField = new JTextArea(DEFAULT_FILE_NAME_STEM + ".log");
    JTextArea treeFileNameField = new JTextArea(DEFAULT_FILE_NAME_STEM + "." + STARBEASTOptions.TREE_FILE_NAME);
//    JCheckBox allowOverwriteLogCheck = new JCheckBox("Allow to overwrite the existing log file");

//    JCheckBox mapTreeLogCheck = new JCheckBox("Create tree file containing the MAP tree:");
//    JTextField mapTreeFileNameField = new JTextField("untitled.MAP.tree");

    JCheckBox substTreeLogCheck = new JCheckBox("Create tree log file with branch length in substitutions:");
    JTextArea substTreeFileNameField = new JTextArea("untitled(subst).trees");

    JCheckBox operatorAnalysisCheck = new JCheckBox("Create operator analysis file:");
    JTextArea operatorAnalysisFileNameField = new JTextArea(DEFAULT_FILE_NAME_STEM + ".ops");

    BeautiFrame frame = null;
    private final OptionsPanel optionsPanel;
    private BeautiOptions options;

    private MLEDialog mleDialog = null;
    private MLEGSSDialog mleGssDialog = null;
    private MarginalLikelihoodEstimationOptions mleOptions = new MarginalLikelihoodEstimationOptions();

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

        optionsPanel.addComponent(operatorAnalysisCheck);
        operatorAnalysisCheck.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(ActionEvent e) {
                options.operatorAnalysis = operatorAnalysisCheck.isSelected();

                updateOtherFileNames(options);

                frame.setDirty();
            }
        });

        operatorAnalysisFileNameField.setColumns(32);
        operatorAnalysisFileNameField.setEditable(false);
        operatorAnalysisFileNameField.setEnabled(false);
        optionsPanel.addComponentWithLabel("Operator analysis file name:", operatorAnalysisFileNameField);

        optionsPanel.addSeparator();

        optionsPanel.addComponent(samplePriorCheckBox);
        samplePriorCheckBox.setOpaque(false);
        samplePriorCheckBox.addChangeListener(new ChangeListener() {
            public void stateChanged(ChangeEvent changeEvent) {
                frame.setDirty();
            }
        });

        optionsPanel.addSeparator();

        JTextArea mleInfo = new JTextArea("Select the option below to perform marginal likelihood " +
                "estimation (MLE) using path sampling (PS) / stepping-stone sampling (SS) " +
                "or generalized stepping-stone sampling (GSS) which performs an additional" +
                "analysis after the standard MCMC chain has finished.");
        mleInfo.setColumns(50);
        PanelUtils.setupComponent(mleInfo);
        optionsPanel.addSpanningComponent(mleInfo);

        //add PS/SS button
        optionsPanel.addComponentWithLabel("Marginal likelihood estimation (MLE):", performMLECombo);

        optionsPanel.addComponent(buttonMLE);
        buttonMLE.setEnabled(false);
        performMLECombo.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(ActionEvent e) {
                if (performMLECombo.getSelectedIndex() == 1) {
                    mleOptions.performMLE = true;
                    mleOptions.performMLEGSS = false;
                    options.logCoalescentEventsStatistic = false;
                    buttonMLE.setEnabled(true);
                    updateMLEFileNameStem();
                } else if (performMLECombo.getSelectedIndex() == 2) {
                    // Generalized stepping-stone sampling
                    for (AbstractPartitionData partition : options.getDataPartitions()) {
                        if (partition.getDataType().getType() != DataType.NUCLEOTIDES) {
                            JOptionPane.showMessageDialog(frame,
                                    "Generalized stepping-stone sampling is not currently\n" +
                                            "compatible with substitution models other than those\n" +
                                            "for nucleotide data. \n\n" +
                                            "Use path sampling/stepping-stone sampling instead",
                                    "Warning",
                                    JOptionPane.WARNING_MESSAGE);

                            performMLECombo.setSelectedIndex(1);
                            return;

                        }
                    }
                    mleOptions.performMLE = false;
                    mleOptions.performMLEGSS = true;
                    //set to true because product of exponentials is the default option
                    options.logCoalescentEventsStatistic = true;
                    buttonMLE.setEnabled(true);
                    updateMLEFileNameStem();
                } else {
                    mleOptions.performMLE = false;
                    mleOptions.performMLEGSS = false;
                    mleOptions.printOperatorAnalysis = false;
                    options.logCoalescentEventsStatistic = false;
                    buttonMLE.setEnabled(false);
                }
            }
        });
        buttonMLE.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(ActionEvent e) {
                updateMLEFileNameStem();

                int result = -1;

                if (performMLECombo.getSelectedIndex() == 1) {
                    result  = mleDialog.showDialog();
                } else if (performMLECombo.getSelectedIndex() == 2) {
                    result = mleGssDialog.showDialog();
                }

                if (result == -1 || result == JOptionPane.CANCEL_OPTION) {
                    return;
                }

            }
        });

//        JTextArea mleGssInfo = new JTextArea("Select the option below to perform marginal likelihood " +
//                "estimation (MLE) using generalized stepping-stone sampling (GSS) which " +
//                "performs an additional analysis after the standard MCMC chain has finished.");
//        mleGssInfo.setColumns(50);
//        PanelUtils.setupComponent(mleGssInfo);
//        optionsPanel.addSpanningComponent(mleGssInfo);

//        //add GSS button
//        optionsPanel.addComponent(performMLEGSS);
//        //will be false by default
//        //options.performMLE = false; ??
//        optionsPanel.addComponent(buttonMLEGSS);
//        buttonMLEGSS.setEnabled(false);
//        performMLEGSS.addActionListener(new java.awt.event.ActionListener() {
//            public void actionPerformed(ActionEvent e) {
//                if (performMLEGSS.isSelected()) {
//                    mleOptions.performMLEGSS = true;
//                    //set to true because product of exponentials is the default option
//                    options.logCoalescentEventsStatistic = true;
//                    buttonMLEGSS.setEnabled(true);
//                    buttonMLE.setEnabled(false);
//                    performMLE.setEnabled(false);
//                    updateMLEFileNameStem();
//                } else {
//                    mleOptions.performMLEGSS = false;
//                    mleOptions.printOperatorAnalysis = false;
//                    options.logCoalescentEventsStatistic = false;
//                    buttonMLE.setEnabled(false);
//                    performMLE.setEnabled(true);
//                    buttonMLEGSS.setEnabled(false);
//                }
//            }
//        });
//        buttonMLEGSS.addActionListener(new java.awt.event.ActionListener() {
//            public void actionPerformed(ActionEvent e) {
//                updateMLEFileNameStem();
//
//                int result = mleGssDialog.showDialog();
//
//                if (result == -1 || result == JOptionPane.CANCEL_OPTION) {
//                    return;
//                }
//
//            }
//        });

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
        String treeFileName;

        for (PartitionTreeModel tree : options.getPartitionTreeModels()) {
            if (options.substTreeLog) {
                treeFileName = getTreeFileName(tree.getPrefix() + "(time).");
            } else {
                treeFileName = getTreeFileName(tree.getPrefix());
            }
            if (addTxt.isSelected()) treeFileName = treeFileName + ".txt";
            options.treeFileName.add(treeFileName);

            if (options.substTreeLog) {
                treeFileName = getTreeFileName(tree.getPrefix() + "(subst).");
                if (addTxt.isSelected()) treeFileName = treeFileName + ".txt";
                options.substTreeFileName.add(treeFileName);
            }
        }

        if (options.useStarBEAST) {
            treeFileName = options.fileNameStem + "." + options.starBEASTOptions.SPECIES_TREE_FILE_NAME;
            if (addTxt.isSelected()) treeFileName = treeFileName + ".txt";
            options.treeFileName.add(treeFileName);
            //TODO: species sub tree
        }
    }

    private String getTreeFileName(String treeName) {
        return options.fileNameStem + "." + treeName + STARBEASTOptions.TREE_FILE_NAME;
    }

    private String displayTreeList(List<String> treeList) {
        if (treeList.size() > 1) {
            return getTreeFileName("[tree name].");
        } else {
            return getTreeFileName("");
        }
    }

    public void setOptions(BeautiOptions options) {
        this.options = options;

        /*System.err.println("mleOptions: " + mleOptions);
        System.err.println("options.pathSteps: " + mleOptions.pathSteps);
        System.err.println("options.mleChainLength: " + mleOptions.mleChainLength);
        System.err.println("options.mleLogEvery: " + mleOptions.mleLogEvery);*/

        // get the MLE options
        mleOptions = (MarginalLikelihoodEstimationOptions)options.getComponentOptions(MarginalLikelihoodEstimationOptions.class);
        /*if (mleOptions == null) {
            mleOptions = new MarginalLikelihoodEstimationOptions();
        }*/

        /*System.err.println("mleOptions: " + mleOptions);
        System.err.println("options.pathSteps: " + mleOptions.pathSteps);
        System.err.println("options.mleChainLength: " + mleOptions.mleChainLength);
        System.err.println("options.mleLogEvery: " + mleOptions.mleLogEvery);*/

        if (mleDialog != null) {
            //mleDialog.setOptions(mleOptions);
        }
        if (mleGssDialog != null) {
            //mleGssDialog.setOptions(mleOptions);
        }

        chainLengthField.setValue(options.chainLength);
        echoEveryField.setValue(options.echoEvery);
        logEveryField.setValue(options.logEvery);

        if (options.fileNameStem != null) {
            fileNameStemField.setText(options.fileNameStem);
        } else {
            fileNameStemField.setText(DEFAULT_FILE_NAME_STEM);
            fileNameStemField.setEnabled(false);
        }

        operatorAnalysisCheck.setSelected(options.operatorAnalysis);

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
            operatorAnalysisFileNameField.setEnabled(options.operatorAnalysis);
            if (options.operatorAnalysis) {
                operatorAnalysisFileNameField.setText(options.operatorAnalysisFileName);
            } else {
                operatorAnalysisFileNameField.setText("");
            }

//            mapTreeLogCheck.setEnabled(true);
//            mapTreeLogCheck.setSelected(options.mapTreeLog);
//            mapTreeFileNameField.setEnabled(options.mapTreeLog);

            substTreeLogCheck.setEnabled(true);
            substTreeLogCheck.setSelected(options.substTreeLog);

            updateMLEFileNameStem();
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
            operatorAnalysisCheck.setSelected(false);
            operatorAnalysisFileNameField.setText("");
        }
    }

    private void updateMLEFileNameStem() {
        if (mleDialog == null) {
            mleDialog = new MLEDialog(frame, mleOptions);
        }
        if (mleGssDialog == null) {
            mleGssDialog = new MLEGSSDialog(frame, mleOptions, options);
        }
        mleDialog.setFilenameStem(options.fileNameStem, addTxt.isSelected());
        mleGssDialog.setFilenameStem(options.fileNameStem, addTxt.isSelected());
    }

    public void getOptions(BeautiOptions options) {
        options.fileNameStem = fileNameStemField.getText();
        options.logFileName = logFileNameField.getText();

//        options.mapTreeLog = mapTreeLogCheck.isSelected();
//        options.mapTreeFileName = mapTreeFileNameField.getText();

        options.substTreeLog = substTreeLogCheck.isSelected();
        updateTreeFileNameList();

        options.operatorAnalysis = operatorAnalysisCheck.isSelected();
        options.operatorAnalysisFileName = operatorAnalysisFileNameField.getText();

        options.samplePriorOnly = samplePriorCheckBox.isSelected();

        if (mleDialog != null) {
            //mleDialog.getOptions(mleOptions);
        }

        if (mleGssDialog != null) {
            //mleGssDialog.getOptions(mleOptions);
        }

    }

    public JComponent getExportableComponent() {
        return optionsPanel;
    }

}
