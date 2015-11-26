/*
 * SimulationsPanel.java
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

package dr.app.coalgen;

import jam.framework.Exportable;
import jam.panels.OptionsPanel;

import javax.swing.*;
import java.awt.*;
import dr.app.gui.components.WholeNumberField;

/**
 * @author Andrew Rambaut
 * @author Alexei Drummond
 * @version $Id$
 */
public class SimulationsPanel extends JPanel implements Exportable {

    private final CoalGenFrame frame;
    private final CoalGenData data;

    private OptionsPanel optionPanel;

    private final JLabel replicatesLabel = new JLabel("Number of replicates:");
    private final WholeNumberField replicatesField = new WholeNumberField(1, Integer.MAX_VALUE);


    public SimulationsPanel(final CoalGenFrame frame, final CoalGenData data) {

        super();

        this.frame = frame;
        this.data = data;

        setOpaque(false);
        setLayout(new BorderLayout());

        optionPanel = new OptionsPanel(12, 12, SwingConstants.CENTER);
        add(optionPanel, BorderLayout.NORTH);

        replicatesField.setColumns(8);
        replicatesField.setValue(data.replicateCount);

        optionPanel.addComponents(replicatesLabel, replicatesField);
    }

    public final void tracesChanged() {
        replicatesLabel.setEnabled(data.traces == null);
        replicatesField.setEnabled(data.traces == null);
        if (data.traces != null) {
            replicatesField.setValue(data.traces.getStateCount());
        }

        replicatesField.setValue(data.replicateCount);
    }

    public final void collectSettings() {
        data.replicateCount = replicatesField.getValue();
    }

    public JComponent getExportableComponent() {
        return this;
    }
}