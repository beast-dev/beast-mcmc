/*
 * PriorDialog.java
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

package dr.app.beauti.sitemodelspanel;

import dr.app.beauti.BeautiFrame;
import dr.app.beauti.options.TraitData;
import dr.app.util.OSType;

import javax.swing.*;
import javax.swing.border.EmptyBorder;
import java.awt.*;

/**
 * @author Andrew Rambaut
 */
public class GLMSettingsDialog  {

    private final BeautiFrame frame;

    private final GLMSettingsPanel glmSettingsPanel;
    private TraitData trait;

    public GLMSettingsDialog(BeautiFrame frame) {
        this.frame = frame;

        glmSettingsPanel = new GLMSettingsPanel(frame);

    }

    /**
     * Set the trait to be controlled
     *
     * @param trait
     */
    public void setTrait(final TraitData trait) {
        this.trait = trait;
        glmSettingsPanel.setTrait(trait);
    }

    public int showDialog() {

        JPanel panel = new JPanel(new BorderLayout(0, 6));
        panel.add(new JLabel("Set GLM design for " + trait), BorderLayout.NORTH);
        panel.add(glmSettingsPanel, BorderLayout.CENTER);

        JScrollPane scrollPane = new JScrollPane(panel, JScrollPane.VERTICAL_SCROLLBAR_AS_NEEDED, JScrollPane.HORIZONTAL_SCROLLBAR_NEVER);
        scrollPane.setBorder(null);
        scrollPane.getViewport().setOpaque(false);

        Object importName = frame.getImportTraitsAction().getValue(Action.NAME);
        frame.getImportTraitsAction().putValue(Action.NAME, "Import Predictors...");

        JOptionPane optionPane = new JOptionPane(scrollPane,
                JOptionPane.PLAIN_MESSAGE,
                JOptionPane.OK_CANCEL_OPTION,
                null,
                null,
                null);
        optionPane.setBorder(new EmptyBorder(12, 12, 12, 12));

        final JDialog dialog = optionPane.createDialog(frame, "GLM settings for " + trait);

        glmSettingsPanel.setDialog(dialog);

        if (OSType.isMac()) {
            dialog.setMinimumSize(new Dimension(dialog.getBounds().width, 300));
        } else {
            Toolkit tk = Toolkit.getDefaultToolkit();
            Dimension d = tk.getScreenSize();
            if (d.height < 700 && glmSettingsPanel.getHeight() > 450) {
                dialog.setSize(new Dimension(glmSettingsPanel.getWidth() + 100, 550));
            } else {
                // setSize because optionsPanel is shrunk in dialog
                dialog.setSize(new Dimension(glmSettingsPanel.getWidth() + 100, glmSettingsPanel.getHeight() + 100));
            }

//            System.out.println("panel width = " + panel.getWidth());
//            System.out.println("panel height = " + panel.getHeight());
        }

        dialog.pack();
        dialog.setResizable(true);
        dialog.setVisible(true);

        int result = JOptionPane.CANCEL_OPTION;
        Integer value = (Integer) optionPane.getValue();
        if (value != null && value != -1) {
            result = value;
        }

        frame.getImportTraitsAction().putValue(Action.NAME, importName);

        return result;
    }
}