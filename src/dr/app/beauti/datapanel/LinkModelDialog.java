/*
 * LinkModelDialog.java
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

import jam.panels.OptionsPanel;

import javax.swing.*;
import javax.swing.border.EmptyBorder;
import java.awt.event.ItemEvent;
import java.awt.event.ItemListener;

/**
 * Asks which model the selected partitions should share, and whether to rename it.
 *
 * @param <T> the kind of model being shared
 * @author Andrew Rambaut
 * @author Walter Xie
 */
public class LinkModelDialog<T> {

    private final JFrame frame;
    private final String title;

    private final JComboBox modelCombo;
    private final JCheckBox copyCheck;
    private final JTextField nameField;
    private final OptionsPanel optionPanel;

    /**
     * @param label    names the list of models (e.g. "Partition clock model:")
     * @param renaming names the rename option (e.g. "Rename clock model partition to:")
     * @param title    the dialog's title
     */
    public LinkModelDialog(JFrame frame, String label, String renaming, String title) {
        this.frame = frame;
        this.title = title;

        modelCombo = new JComboBox();

        copyCheck = new JCheckBox(renaming);
        nameField = new JTextField();
        nameField.setColumns(20);
        nameField.setEnabled(false);

        optionPanel = new OptionsPanel(12, 12);
        optionPanel.addComponentWithLabel(label, modelCombo);
        optionPanel.addComponents(copyCheck, nameField);

        copyCheck.addItemListener(new ItemListener() {
            public void itemStateChanged(ItemEvent ev) {
                nameField.setEnabled(copyCheck.isSelected());
            }
        });
    }

    @SuppressWarnings("unchecked")
    public int showDialog(Object[] models) {

        modelCombo.removeAllItems();
        for (Object model : models) {
            modelCombo.addItem(model);
        }

        JOptionPane optionPane = new JOptionPane(optionPanel,
                JOptionPane.QUESTION_MESSAGE,
                JOptionPane.OK_CANCEL_OPTION,
                null,
                null,
                null);
        optionPane.setBorder(new EmptyBorder(12, 12, 12, 12));

        final JDialog dialog = optionPane.createDialog(frame, title);
        dialog.pack();

        dialog.setVisible(true);

        // the dialog is modal, so it has been closed by the time this returns
        dialog.dispose();

        int result = JOptionPane.CANCEL_OPTION;
        Integer value = (Integer) optionPane.getValue();
        if (value != null && value != -1) {
            result = value;
        }

        return result;
    }

    @SuppressWarnings("unchecked")
    public T getModel() {
        return (T) modelCombo.getSelectedItem();
    }

    public boolean getMakeCopy() {
        return copyCheck.isSelected();
    }

    public String getName() {
        return nameField.getText();
    }
}
