/*
 * MCMCPanel.java
 *
 * Copyright (C) 2002-2009 Alexei Drummond and Andrew Rambaut
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
 * BEAST is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with BEAST; if not, write to the
 * Free Software Foundation, Inc., 51 Franklin St, Fifth Floor,
 * Boston, MA  02110-1301  USA
 */

package dr.app.phylogeography.spread;

import dr.app.beauti.BeautiFrame;
import dr.app.beauti.BeautiPanel;
import dr.app.beauti.options.BeautiOptions;
import dr.app.beauti.options.PartitionTreeModel;
import dr.app.util.OSType;
import dr.app.phylogeography.generator.Generator;
//import dr.evomodel.coalescent.GMRFFixedGridImportanceSampler;
import org.virion.jam.components.WholeNumberField;
import org.virion.jam.panels.OptionsPanel;

import javax.swing.*;
import javax.swing.event.ChangeEvent;
import javax.swing.event.ChangeListener;
import java.awt.*;
import java.awt.event.ActionEvent;
import java.awt.event.KeyEvent;
import java.util.List;

/**
 * @author Andrew Rambaut
 * @version $Id$
 */
public class TimelinePanel extends JPanel {

    private static final long serialVersionUID = -3710586474593827540L;

    private final SpreadFrame frame;

    private final SpreadDocument document;

    public TimelinePanel(final SpreadFrame parent, final SpreadDocument document) {

        this.frame = parent;
        this.document = document;

        setLayout(new BorderLayout());

        OptionsPanel optionsPanel = new OptionsPanel(12, 24);

        setOpaque(false);
        optionsPanel.setOpaque(false);

        optionsPanel.addSeparator();


        add(optionsPanel, BorderLayout.CENTER);

    }

}