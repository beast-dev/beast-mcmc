/*
 * OptionsPanel.java
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

package org.virion.jam.panels;

import javax.swing.*;
import java.awt.*;

/**
 * OptionsPanel.
 *
 * @author Andrew Rambaut
 * @version $Id: OptionsPanel.java,v 1.2 2006/09/09 18:16:16 rambaut Exp $
 */
public class OptionsPanel extends JPanel {

	/**
	 * 
	 */
	private static final long serialVersionUID = -7217075209325630197L;
	protected GridBagLayout gridbag = new GridBagLayout();
	private final int hGap;
	private final int vGap;

	public OptionsPanel() {
		this(4, 4);
	}

	public OptionsPanel(int hGap, int vGap) {
		this.hGap = hGap;
		this.vGap = vGap;

		setBorder(BorderFactory.createEmptyBorder(3,6,6,6));
		setLayout(gridbag);
		setOpaque(false);
	}

	public void addLabel(String text) {

		addSpanningComponent(new JLabel(text));
	}

	public void addSpanningComponent(JComponent comp) {

		adjustComponent(comp);

		GridBagConstraints c = new GridBagConstraints();
		c.insets = new Insets(vGap / 2, 0, vGap / 2, 0);
		c.weightx = 1.0;
		c.fill = GridBagConstraints.HORIZONTAL;

		c.gridwidth = GridBagConstraints.REMAINDER;

		// it is now up to the calling code to add a border
		//comp.setBorder(new EmptyBorder(6, 24, 6, 24));
		gridbag.setConstraints(comp, c);
		add(comp);
	}

    /**
     * This was a spelling mistake but remains here to avoid breaking
     * too many things. Marked as depriciated
     * @deprecated
     */
    public void addSeperator() {
        addSeparator();
    }

	public void addSeparator() {

		GridBagConstraints c = new GridBagConstraints();
		c.insets = new Insets(0,0,0,0);
		c.weightx = 1.0;
		c.fill = GridBagConstraints.HORIZONTAL;

		c.gridwidth = GridBagConstraints.REMAINDER;
		JSeparator separator = new JSeparator();
		adjustComponent(separator);
		separator.setOpaque(false);
		gridbag.setConstraints(separator, c);
		add(separator);
	}

	public void addComponent(JComponent comp) {
		addComponent(comp, false);
	}

	public void addComponent(JComponent comp, boolean fill) {
		JPanel panel = new JPanel();
		panel.setOpaque(false);
		addComponents(panel, false, comp, fill);
	}

	public JLabel addComponentWithLabel(String text, JComponent comp) {
		return addComponentWithLabel(text, comp, false);
	}

	public JLabel addComponentWithLabel(String text, JComponent comp, boolean fill) {

		JLabel label = new JLabel(text, javax.swing.SwingConstants.RIGHT);
		label.setLabelFor(comp);
		label.setOpaque(false);
		addComponents(label, false, comp, fill);

		return label;
	}

	public void addComponents(JComponent comp1, JComponent comp2) {
		addComponents(comp1, false, comp2, false);
	}

	public void addComponents(JComponent comp1, boolean fill1, JComponent comp2, boolean fill2) {

		adjustComponent(comp1);
		adjustComponent(comp2);

		GridBagConstraints c = new GridBagConstraints();
		c.weightx = 0.0;
		c.fill = fill1 ? GridBagConstraints.HORIZONTAL : GridBagConstraints.NONE;

		c.insets = new Insets(vGap / 2, 0, vGap / 2, hGap / 2);
		c.gridwidth = 1;
		c.anchor = GridBagConstraints.EAST;
		gridbag.setConstraints(comp1, c);
		add(comp1);

		c.weightx = 1.0;
		c.insets = new Insets(vGap / 2, hGap / 2, vGap / 2, 0);
		c.fill = fill2 ? GridBagConstraints.HORIZONTAL : GridBagConstraints.NONE;
		c.gridwidth = GridBagConstraints.REMAINDER;
		c.anchor = GridBagConstraints.WEST;
		gridbag.setConstraints(comp2, c);
		add(comp2);
	}

	private void adjustComponent(JComponent comp) {
        comp.putClientProperty("Quaqua.Component.visualMargin", new Insets(0,0,0,0));
//        if (comp.getFont() != null) {
//            comp.setFont(comp.getFont().deriveFont(11.0f));
//	}
}
}
