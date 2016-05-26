/*
 * SelectTreeDialog.java
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

import dr.app.beauti.options.PartitionTreeModel;
import jam.panels.OptionsPanel;

import javax.swing.*;
import javax.swing.border.EmptyBorder;
import java.awt.event.*;

/**
 * @author Andrew Rambaut
 * @version $Id:$
 */
public class SelectTreeDialog {

    private JFrame frame;

    JComboBox treeCombo;
    JCheckBox copyCheck;
    JTextField nameField;

    OptionsPanel optionPanel;

    public SelectTreeDialog(JFrame frame) {
        this.frame = frame;

        treeCombo = new JComboBox();

        copyCheck = new JCheckBox("Rename tree model partition to:");
        nameField = new JTextField();
        nameField.setColumns(20);
        nameField.setEnabled(false);

        optionPanel = new OptionsPanel(12, 12);
        optionPanel.addComponentWithLabel("Partition tree:", treeCombo);
        optionPanel.addComponents(copyCheck, nameField);

        copyCheck.addItemListener(
                new ItemListener() {
                    public void itemStateChanged(ItemEvent ev) {
                        nameField.setEnabled(copyCheck.isSelected());
                    }
                }
        );

    }

    public int showDialog(Object[] trees) {

        treeCombo.removeAllItems();
        for (Object tree : trees) {
            treeCombo.addItem(tree);
        }

        JOptionPane optionPane = new JOptionPane(optionPanel,
                JOptionPane.QUESTION_MESSAGE,
                JOptionPane.OK_CANCEL_OPTION,
                null,
                null,
                null);
        optionPane.setBorder(new EmptyBorder(12, 12, 12, 12));

        final JDialog dialog = optionPane.createDialog(frame, "Create New Partition Tree");
        dialog.pack();

        dialog.setVisible(true);

        int result = JOptionPane.CANCEL_OPTION;
        Integer value = (Integer) optionPane.getValue();
        if (value != null && value != -1) {
            result = value;
        }

        return result;
    }

    public PartitionTreeModel getTree() {
        return (PartitionTreeModel) treeCombo.getSelectedItem();
    }

    public boolean getMakeCopy() {
        return copyCheck.isSelected();
    }

    public String getName() {
        return nameField.getText();
    }

}