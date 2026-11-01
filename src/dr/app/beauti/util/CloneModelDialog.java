/*
 * CloneModelDialog.java
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

package dr.app.beauti.util;

import jam.panels.OptionsPanel;

import javax.swing.*;
import javax.swing.border.EmptyBorder;
import java.util.List;

/**
 * Asks which of the selected models the others should take their settings from.
 *
 * @param <T> the kind of model being cloned
 * @author Andrew Rambaut
 */
public class CloneModelDialog<T> {

    private final JFrame frame;
    private final JComboBox sourceModelCombo;
    private final OptionsPanel optionPanel;

    /**
     * @param modelDescription what the models are, for the wording (e.g. "substitution model")
     */
    public CloneModelDialog(JFrame frame, String modelDescription) {
        this.frame = frame;

        String question = "<html>Select the " + modelDescription + " to act as a source<br>" +
                "to copy to the other selected models.</html>";

        sourceModelCombo = new JComboBox();
        PanelUtils.setupComponent(sourceModelCombo);
        sourceModelCombo.setToolTipText(question);

        optionPanel = new OptionsPanel(12, 12);
        optionPanel.addSpanningComponent(new JLabel(question));
        optionPanel.addComponentWithLabel("Source Model:", sourceModelCombo);
    }

    @SuppressWarnings("unchecked")
    public int showDialog(List<? extends T> sourceModels) {

        JOptionPane optionPane = new JOptionPane(optionPanel,
                JOptionPane.QUESTION_MESSAGE,
                JOptionPane.OK_CANCEL_OPTION,
                null,
                null,
                null);
        optionPane.setBorder(new EmptyBorder(12, 12, 12, 12));

        sourceModelCombo.removeAllItems();
        for (T model : sourceModels) {
            sourceModelCombo.addItem(model);
        }

        final JDialog dialog = optionPane.createDialog(frame, "Clone model settings");
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
    public T getSourceModel() {
        return (T) sourceModelCombo.getSelectedItem();
    }
}
