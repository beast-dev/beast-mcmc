/*
 * SelectTraitDialog.java
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

package dr.app.beauti.datapanel;

import dr.app.beauti.options.TraitData;
import jam.panels.OptionsPanel;

import javax.swing.*;
import javax.swing.border.EmptyBorder;
import javax.swing.event.ListSelectionEvent;
import javax.swing.event.ListSelectionListener;
import java.awt.*;
import java.awt.event.ItemEvent;
import java.awt.event.ItemListener;
import java.util.Collection;
import java.util.List;

/**
 * @author Andrew Rambaut
 * @author Alexei Drummond
 * @version $Id: PriorDialog.java,v 1.4 2006/09/05 13:29:34 rambaut Exp $
 */
public class SelectTraitDialog {

    private JFrame frame;

    //    JComboBox traitCombo;
    JList traitList;
    DefaultListModel traitModel;
    private final JCheckBox copyCheck;
    JTextField nameField;
    JScrollPane scrollPane;
    private final JCheckBox independentBox;

    OptionsPanel optionPanel;

    public SelectTraitDialog(JFrame frame) {
        this.frame = frame;

        traitModel = new DefaultListModel();
        traitList = new JList(traitModel);
        scrollPane = new JScrollPane(traitList);

        copyCheck = new JCheckBox("Name trait partition:");
        nameField = new JTextField();
        nameField.setColumns(20);

        independentBox = new JCheckBox("New partition for each trait");

        optionPanel = new OptionsPanel(12, 12);

        copyCheck.addItemListener(
                new java.awt.event.ItemListener() {
                    public void itemStateChanged(java.awt.event.ItemEvent ev) {
                        nameField.setEnabled(copyCheck.isSelected());
                        if (getSelectedTraitCount() > 1) independentBox.setSelected(!copyCheck.isSelected());
                    }
                }
        );

        traitList.addListSelectionListener(new ListSelectionListener() {
            @Override
            public void valueChanged(ListSelectionEvent e) {
                if (getSelectedTraitCount() > 1) {
                    if (independentBox.isSelected() && copyCheck.isSelected()) independentBox.setSelected(false);
                    if (!(independentBox.isSelected() || copyCheck.isSelected())) copyCheck.setSelected(true);
                }
            }
        });

        independentBox.addItemListener(new ItemListener() {
            @Override
            public void itemStateChanged(ItemEvent e) {
                if (getSelectedTraitCount() > 1) copyCheck.setSelected(!(independentBox.isSelected()));
            }
        });

    }

    public int showDialog(Collection<TraitData> traits, String defaultName, Component parent, Boolean allowSelectTraits) {
        optionPanel.removeAll();
        nameField.setText(defaultName != null ? defaultName : "untitled_traits");

        if (traits == null || !allowSelectTraits) { // traits shouldn't be null
            optionPanel.addSpanningComponent(new JLabel("Create a new data partition using the selected trait(s)."));
            optionPanel.addComponentWithLabel("Name trait partition:", nameField);
            nameField.setEnabled(true);
            nameField.selectAll();
        } else {
            traitModel.removeAllElements();

            for (Object model : traits) {
                traitModel.addElement(model);
            }
            optionPanel.addSpanningComponent(new JLabel("Create a new data partition using the following trait(s)."));
            optionPanel.addComponentWithLabel("Trait(s):", scrollPane);
            optionPanel.addComponents(copyCheck, nameField);
            nameField.setEnabled(copyCheck.isSelected());
        }

        if (traits != null && traits.size() > 1) { // traits shouldn't be null
            independentBox.setSelected(false);
            optionPanel.addComponent(independentBox);
        }

        JOptionPane optionPane = new JOptionPane(optionPanel,
                JOptionPane.QUESTION_MESSAGE,
                JOptionPane.OK_CANCEL_OPTION,
                null,
                null,
                null);
        optionPane.setBorder(new EmptyBorder(12, 12, 12, 12));

        final JDialog dialog = optionPane.createDialog(frame, "Create New Partition");
        dialog.pack();

        int result;
        boolean isValid;

        do {
            dialog.setVisible(true);

            isValid = true;
            result = JOptionPane.CANCEL_OPTION;

            Integer value = (Integer) optionPane.getValue();
            if (value != null && value != -1) {
                result = value;
            }
            if (result != JOptionPane.CANCEL_OPTION) {
                String name = getName().trim();
                if (name.isEmpty() && (getMakeCopy() || traits == null)) {
                    isValid = false;
                    JOptionPane.showMessageDialog(parent, "Cannot have an empty partition name.", "No Partition Name", JOptionPane.ERROR_MESSAGE);
                }

                if (getTraits().size() == 0 && allowSelectTraits) {
                    isValid = false;
                    JOptionPane.showMessageDialog(parent, "Please select a trait(s).", "No Trait(s) Selected", JOptionPane.ERROR_MESSAGE);
                }

                if (getTraits().size() > 1 && !(getMakeCopy() || getForceIndependent())) { // pretty sure this can't happen, but warning just in case
                    isValid = false;
                    JOptionPane.showMessageDialog(parent,
                            "You have selected multiple traits.\nPlease either name the new partition or " +
                                    "check \"" + independentBox.getText() + "\"",
                            "Multiple Traits Selected Without Name", JOptionPane.ERROR_MESSAGE);
                }
            }
        } while (!isValid);

        return result;
    }

    public List<TraitData> getTraits() {
        return (List<TraitData>) traitList.getSelectedValuesList();
    }

    private int getSelectedTraitCount() {
        return traitList.getSelectedValuesList().size();
    }

    public boolean getMakeCopy() {
        return copyCheck.isSelected();
    }

    public void reset() {
        copyCheck.setSelected(false);
        independentBox.setSelected(false);
    }

    public boolean getForceIndependent() {
        return independentBox.isSelected();
    }

    public String getName() {
        return nameField.getText();
    }

}