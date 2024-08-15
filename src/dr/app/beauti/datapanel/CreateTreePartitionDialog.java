/*
 * CreateTreePartitionDialog.java
 *
 * Copyright © 2002-2024 the BEAST Development Team
 * http://beast.community/about
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
 *
 */

package dr.app.beauti.datapanel;

import dr.app.beauti.options.TreeHolder;
import dr.evolution.tree.Tree;
import jam.panels.OptionsPanel;

import javax.swing.*;
import javax.swing.border.EmptyBorder;
import java.awt.*;
import java.util.Map;

/**
 * @author Andrew Rambaut
 * @author Alexei Drummond
 */
public class CreateTreePartitionDialog {

    private JFrame frame;

    //    JComboBox traitCombo;
    JComboBox treeCombo;
    JCheckBox renameCheck;
    JTextField nameField;

    OptionsPanel optionPanel;

    public CreateTreePartitionDialog(JFrame frame) {
        this.frame = frame;

        treeCombo = new JComboBox();

        optionPanel = new OptionsPanel(12, 12);
        optionPanel.addComponentWithLabel("Create partition from:", treeCombo);

        renameCheck = new JCheckBox("Name tree partition:");
        nameField = new JTextField();
        nameField.setColumns(20);
        nameField.setEnabled(false);

        optionPanel.addComponents(renameCheck, nameField);

        renameCheck.addItemListener(
                ev -> nameField.setEnabled(renameCheck.isSelected())
        );


    }

    public int showDialog(Map<String, TreeHolder> trees, String defaultName, Component parent) {
        treeCombo.removeAllItems();
        for (TreeHolder tree : trees.values()) {
            treeCombo.addItem(tree);
        }

        JOptionPane optionPane = new JOptionPane(optionPanel,
                JOptionPane.QUESTION_MESSAGE,
                JOptionPane.OK_CANCEL_OPTION,
                null,
                null,
                null);
        optionPane.setBorder(new EmptyBorder(12, 12, 12, 12));

        final JDialog dialog = optionPane.createDialog(frame, "Create New Partition from Tree");
        dialog.pack();

        dialog.setVisible(true);

        int result = JOptionPane.CANCEL_OPTION;
        Integer value = (Integer) optionPane.getValue();
        if (value != null && value != -1) {
            result = value;
        }

        return result;
    }

    public TreeHolder getTree() {
        return (TreeHolder) treeCombo.getSelectedItem();
    }

    public boolean getRenamePartition() {
        return renameCheck.isSelected();
    }

    public String getName() {
        return nameField.getText();
    }

}