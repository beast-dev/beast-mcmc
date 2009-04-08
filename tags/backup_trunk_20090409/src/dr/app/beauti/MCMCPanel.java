/*
 * MCMCPanel.java
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

package dr.app.beauti;

import dr.app.beauti.options.BeautiOptions;
import org.virion.jam.components.WholeNumberField;
import org.virion.jam.framework.Exportable;
import org.virion.jam.panels.OptionsPanel;

import javax.swing.*;
import javax.swing.event.ChangeEvent;
import javax.swing.event.ChangeListener;
import java.awt.event.ActionEvent;
import java.awt.*;

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

    JTextField logFileNameField = new JTextField("untitled.log");
    JTextField treeFileNameField = new JTextField("untitled.trees");
    JCheckBox mapTreeLogCheck = new JCheckBox("Create tree file containing the MAP tree:");
    JTextField mapTreeFileNameField = new JTextField("untitled.MAP.tree");
    JCheckBox substTreeLogCheck = new JCheckBox("Create tree log file with branch length in substitutions:");
    JTextField substTreeFileNameField = new JTextField("untitled(subst).trees");

    BeautiFrame frame = null;
    private final OptionsPanel optionsPanel;

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

        logFileNameField.setColumns(32);
        optionsPanel.addComponentWithLabel("Log file name:", logFileNameField);
        treeFileNameField.setColumns(32);
        optionsPanel.addComponentWithLabel("Trees file name:", treeFileNameField);

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
                substTreeFileNameField.setEnabled(substTreeLogCheck.isSelected());
                frame.setDirty();
            }
        });

        substTreeFileNameField.setColumns(32);
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
        logFileNameField.addKeyListener(listener);
        treeFileNameField.addKeyListener(listener);
        //mapTreeFileNameField.addKeyListener(listener);
        substTreeFileNameField.addKeyListener(listener);

        add(optionsPanel, BorderLayout.CENTER);

    }

    public void setOptions(BeautiOptions options) {

        chainLengthField.setValue(options.chainLength);

        echoEveryField.setValue(options.echoEvery);
        logEveryField.setValue(options.logEvery);

        if (options.fileNameStem != null) {
            if (options.logFileName == null) {
                logFileNameField.setText(options.fileNameStem + ".log");
            } else {
                logFileNameField.setText(options.logFileName);
            }
            if (options.treeFileName == null) {
                treeFileNameField.setText(options.fileNameStem + ".trees");
            } else {
                treeFileNameField.setText(options.treeFileName);
            }
//            if (options.mapTreeFileName == null) {
//			    mapTreeFileNameField.setText(options.fileNameStem + ".MAP.tree");
//            } else {
//                mapTreeFileNameField.setText(options.mapTreeFileName);
//            }
            if (options.substTreeFileName == null) {
                substTreeFileNameField.setText(options.fileNameStem + "(subst).trees");
            } else {
                substTreeFileNameField.setText(options.substTreeFileName);
            }
            logFileNameField.setEnabled(true);
            treeFileNameField.setEnabled(true);

//            mapTreeLogCheck.setEnabled(true);
//            mapTreeLogCheck.setSelected(options.mapTreeLog);
//            mapTreeFileNameField.setEnabled(options.mapTreeLog);

            substTreeLogCheck.setEnabled(true);
            substTreeLogCheck.setSelected(options.substTreeLog);
            substTreeFileNameField.setEnabled(options.substTreeLog);
        } else {
            logFileNameField.setText("untitled");
            logFileNameField.setEnabled(false);
            treeFileNameField.setText("untitled");
            treeFileNameField.setEnabled(false);
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

        options.logFileName = logFileNameField.getText();
        options.treeFileName = treeFileNameField.getText();

//        options.mapTreeLog = mapTreeLogCheck.isSelected();
//        options.mapTreeFileName = mapTreeFileNameField.getText();

        options.substTreeLog = substTreeLogCheck.isSelected();
        options.substTreeFileName = substTreeFileNameField.getText();

        options.samplePriorOnly = samplePriorCheckBox.isSelected();
    }

    public JComponent getExportableComponent() {
        return optionsPanel;
    }

}
