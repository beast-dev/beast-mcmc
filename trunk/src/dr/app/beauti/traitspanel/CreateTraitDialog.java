/*
 * PriorDialog.java
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

package dr.app.beauti.traitspanel;

import dr.app.beauti.BeautiFrame;
import dr.app.beauti.options.STARBEASTOptions;
import dr.app.beauti.options.TraitData;
import dr.app.beauti.util.TextUtil;
import jam.panels.OptionsPanel;

import javax.swing.*;
import javax.swing.border.EmptyBorder;
import java.awt.*;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.event.ItemEvent;
import java.awt.event.ItemListener;

/**
 * @author Andrew Rambaut
 * @author Walter Xie
 * @version $Id: PriorDialog.java,v 1.4 2006/09/05 13:29:34 rambaut Exp $
 */
public class CreateTraitDialog {

    private final BeautiFrame frame;

    JTextField nameField;
//    JComboBox nameCombo;
    JComboBox typeCombo;
    private final JRadioButton createRadio = new JRadioButton("Create a new trait", true);
    private final JRadioButton importRadio = new JRadioButton("Import trait(s) from a mapping file", false);
    JButton exampleButton = new JButton("Show example of mapping file format");
    private final JCheckBox createTraitPartitionCheck = new JCheckBox("Create a corresponding data partition", true);

    public static final int OK_IMPORT = 10;

    OptionsPanel optionPanel;

    public CreateTraitDialog(final BeautiFrame frame, String traitName) {
        this.frame = frame;

        if (traitName == null) traitName = "Untitled";
        nameField = new JTextField(traitName);
        nameField.setColumns(20);

//        nameCombo = new JComboBox(TraitData.Traits.values());
        typeCombo = new JComboBox(TraitData.TraitType.values());

        optionPanel = new OptionsPanel(12, 12);

        ButtonGroup group = new ButtonGroup();
        group.add(createRadio);
        group.add(importRadio);
        ItemListener listener = new ItemListener() {
            public void itemStateChanged(ItemEvent e) {
                nameField.setEnabled(createRadio.isSelected());
                typeCombo.setEnabled(createRadio.isSelected());
                exampleButton.setEnabled(importRadio.isSelected());
            }
        };
        createRadio.addItemListener(listener);
        importRadio.addItemListener(listener);

        exampleButton.setEnabled(false);
        exampleButton.addActionListener(new ActionListener() {
            public void actionPerformed(ActionEvent ae) {
                JScrollPane scrollPane = TextUtil.createHTMLScrollPane(STARBEASTOptions.EXAMPLE_FORMAT, new Dimension(400,300));

                JOptionPane.showMessageDialog(frame, scrollPane,
                    "Example of mapping file format",
                    JOptionPane.PLAIN_MESSAGE);
            }
        });

        optionPanel.addComponent(importRadio);
        optionPanel.addComponent(exampleButton);
        optionPanel.addComponent(createRadio);
        optionPanel.addComponentWithLabel("Name:", nameField);
        optionPanel.addComponentWithLabel("Type:", typeCombo);
        optionPanel.addComponent(createTraitPartitionCheck);

    }

    public int showDialog() {

        JOptionPane optionPane = new JOptionPane(optionPanel,
                JOptionPane.QUESTION_MESSAGE,
                JOptionPane.OK_CANCEL_OPTION,
                null,
                null,
                null);
        optionPane.setBorder(new EmptyBorder(12, 12, 12, 12));

        boolean done;
        int result = JOptionPane.CANCEL_OPTION;

        do {
            final JDialog dialog = optionPane.createDialog(frame, "Create or Import Trait(s)");
            dialog.pack();

            dialog.setVisible(true);

            Integer value = (Integer) optionPane.getValue();
            if (value != null && value != -1) {
                result = value;
            }
            done = true;
            if (result == JOptionPane.OK_OPTION && createRadio.isSelected()) {
                done = frame.validateTraitName(getName());
            }
        } while (!done);

        if (importRadio.isSelected() && result == JOptionPane.OK_OPTION) result = OK_IMPORT;

        return result;
    }

    public String getName() {
        return nameField.getText();
//        return nameCombo.getSelectedItem().toString();
    }

    public TraitData.TraitType getType() {
        return (TraitData.TraitType) typeCombo.getSelectedItem();
    }

    public boolean createTraitPartition() {
        return createTraitPartitionCheck.isSelected();
    }

}