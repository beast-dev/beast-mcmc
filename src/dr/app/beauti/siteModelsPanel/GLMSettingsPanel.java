/*
 * PriorSettingsPanel.java
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

package dr.app.beauti.siteModelsPanel;

import dr.app.util.OSType;
import jam.panels.OptionsPanel;

import javax.swing.*;
import java.awt.*;

/**
 * @author Andrew Rambaut
 */
public class GLMSettingsPanel extends JPanel {

    private final JFrame frame;
    private JDialog dialog;

    public GLMSettingsPanel(JFrame frame) {
        this.frame = frame;

    }

    /**
     * Set the trait to be controlled
     *                                                                                                       q
     * @param trait
     */
    public void setTrait(final String trait) {
//        this.parameter = parameter;
//
//        priorCombo = new JComboBox();
//        for (PriorType priorType : PriorType.getPriorTypes(parameter)) {
//            priorCombo.addItem(priorType);
//        }
//
//        if (parameter.priorType != null) {
//            priorCombo.setSelectedItem(parameter.priorType);
//        }
//
//        setupComponents(); // setArguments here
//
//
//        priorCombo.addItemListener(new ItemListener() {
//            public void itemStateChanged(ItemEvent e) {
//                setupComponents();
//                dialog.pack();
//                dialog.repaint();
//            }
//        });
//
//        for (PriorOptionsPanel optionsPanel : optionsPanels.values()) {
//            optionsPanel.removeAllListeners();
//            optionsPanel.addListener(new PriorOptionsPanel.Listener() {
//                public void optionsPanelChanged() {
//                    setupChart();
//                    dialog.pack();
//                    dialog.repaint();
//                }
//            });
//        }
    }

    private void setupComponents() {
        removeAll();

        OptionsPanel optionsPanel = new OptionsPanel(12, (OSType.isMac() ? 6 : 24));

        JPanel panel1 = new JPanel(new FlowLayout(FlowLayout.CENTER));
        panel1.add(optionsPanel);


        repaint();
    }

    public void setDialog(final JDialog dialog) {
        this.dialog = dialog;
    }
}