/*
 * CloneModelDialog.java
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

package dr.app.beauti.clockmodelspanel;

import dr.app.beauti.options.PartitionClockModel;
import dr.app.beauti.util.PanelUtils;
import jam.panels.OptionsPanel;

import javax.swing.*;
import javax.swing.border.EmptyBorder;
import java.util.List;

/**
 * @author Andrew Rambaut
 * @version $Id$
 */
public class CloneModelDialog {

    private JFrame frame;

    JComboBox sourceModelCombo;

    OptionsPanel optionPanel;

    public CloneModelDialog(JFrame frame) {
        this.frame = frame;

        sourceModelCombo = new JComboBox();
        PanelUtils.setupComponent(sourceModelCombo);
        sourceModelCombo
                .setToolTipText("<html>Select the substitution model to act as a source<br>to copy to the other selected models.</html>");

        optionPanel = new OptionsPanel(12, 12);
        optionPanel.addSpanningComponent(new JLabel("<html>Select the substitution model to act as a source<br>to copy to the other selected models.</html>"));
        optionPanel.addComponentWithLabel("Source Model:", sourceModelCombo);
    }

    public int showDialog(List<PartitionClockModel> sourceModels) {

        JOptionPane optionPane = new JOptionPane(optionPanel,
                JOptionPane.QUESTION_MESSAGE,
                JOptionPane.OK_CANCEL_OPTION,
                null,
                null,
                null);
        optionPane.setBorder(new EmptyBorder(12, 12, 12, 12));

        sourceModelCombo.removeAllItems();
        for (PartitionClockModel model : sourceModels) {
            sourceModelCombo.addItem(model);
        }

        final JDialog dialog = optionPane.createDialog(frame, "Clone model settings");
        dialog.pack();

        dialog.setVisible(true);

        int result = JOptionPane.CANCEL_OPTION;
        Integer value = (Integer) optionPane.getValue();
        if (value != null && value != -1) {
            result = value;
        }

        return result;
    }

    public PartitionClockModel getSourceModel() {
        return (PartitionClockModel)sourceModelCombo.getSelectedItem();
    }

}