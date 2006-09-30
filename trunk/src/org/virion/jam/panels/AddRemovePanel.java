/*
 * AddRemovePanel.java
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

import org.virion.jam.util.IconUtils;

import javax.swing.*;
import java.awt.*;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.event.ItemEvent;
import java.awt.event.ItemListener;
import java.util.*;
import java.util.List;

/**
 * OptionsPanel.
 *
 * @author Andrew Rambaut
 * @version $Id: AddRemovePanel.java,v 1.1.1.1 2006/07/16 13:17:56 rambaut Exp $
 */


public abstract class AddRemovePanel extends JPanel {

    private Icon addIcon = null;
	private Icon addRolloverIcon = null;
	private Icon addPressedIcon = null;
    private Icon removeIcon = null;
	private Icon removeRolloverIcon = null;
	private Icon removePressedIcon = null;

    public AddRemovePanel() {
        this(null);
    }

    public AddRemovePanel(JPanel[] panels) {

        try {
            addIcon = IconUtils.getIcon(AddRemovePanel.class, "images/plusminus/plus.png");
	        addRolloverIcon = IconUtils.getIcon(AddRemovePanel.class, "images/plusminus/plusRollover.png");
	        addPressedIcon = IconUtils.getIcon(AddRemovePanel.class, "images/plusminus/plusPressed.png");
            removeIcon = IconUtils.getIcon(AddRemovePanel.class, "images/plusminus/minus.png");
	        removeRolloverIcon = IconUtils.getIcon(AddRemovePanel.class, "images/plusminus/minusRollover.png");
	        removePressedIcon = IconUtils.getIcon(AddRemovePanel.class, "images/plusminus/minusPressed.png");
        } catch (Exception e) {
            // do nothing
        }

        BoxLayout layout = new BoxLayout(this, BoxLayout.PAGE_AXIS);
        setLayout(layout);
        setBorder(BorderFactory.createEmptyBorder(4, 4, 4, 4));
        setOpaque(false);

        if (panels != null && panels.length > 0) {
            JPanel last = null;
            for (int i = 0; i < panels.length; i++) {
                JPanel contents = (JPanel)panels[i];
                addPanel(last, contents);
                last = contents;
            }
        } else {
            addPanel(null);
        }

        Dimension dim = super.getPreferredSize();
        dim.height += getComponent(0).getPreferredSize().getHeight();
        setMinimumSize(dim);
        setPreferredSize(dim);
    }

	public List getPanels() {
		return Collections.unmodifiableList(panels);
	}

    protected abstract JPanel createPanel();

    private void addPanel(JPanel previousPanel) {
        addPanel(previousPanel, createPanel());
    }

    private void addPanel(JPanel previousPanel, JPanel contents) {

	    RowPanel rowPanel = new RowPanel(contents);

        if (previousPanel != null) {
            int index = panels.indexOf(previousPanel);
            add(rowPanel, index + 1);
            panels.add(index + 1, contents);
        } else {
            add(rowPanel, 0);
            panels.add(0, contents);
        }
	    removeAction.setEnabled(panels.size() > 1);
        validate();
	    repaint();
    }

    private void removePanel(JPanel panel) {
        int index = panels.indexOf(panel);
        remove(index);
        panels.remove(index);
	    removeAction.setEnabled(panels.size() > 1);
        validate();
	    repaint();
    }

    class RowPanel extends JPanel {
		RowPanel(final JPanel contents) {

			setLayout(new GridBagLayout());

			setOpaque(true);
			setBorder(BorderFactory.createMatteBorder(1, 1, 1, 1, Color.lightGray));
			setBackground(new Color(0.0F, 0.0F, 0.0F, 0.05F));

			JButton addButton = new JButton("+");
			addButton.putClientProperty("JButton.buttonType", "toolbar");
//			addButton.setBorderPainted(false);
			addButton.setOpaque(false);
			if (addIcon != null) {
			    addButton.setIcon(addIcon);
				addButton.setPressedIcon(addPressedIcon);
				addButton.setRolloverIcon(addRolloverIcon);
			    addButton.setRolloverEnabled(true);
				addButton.setText(null);
			}
			addButton.addActionListener(new ActionListener() {
			    public void actionPerformed(ActionEvent ae) {
			        addPanel(contents);
			    }
			});
			addButton.setEnabled(true);

			JButton removeButton = new JButton(removeAction);
			removeButton.putClientProperty("JButton.buttonType", "toolbar");
//			removeButton.setBorderPainted(false);
			removeButton.setOpaque(false);
			if (removeIcon != null) {
			    removeButton.setIcon(removeIcon);
				removeButton.setPressedIcon(removePressedIcon);
				removeButton.setRolloverIcon(removeRolloverIcon);
			    removeButton.setRolloverEnabled(true);
			    removeButton.setText(null);
			}
			removeButton.addActionListener(new ActionListener() {
			    public void actionPerformed(ActionEvent ae) {
			        removePanel(contents);
			    }
			});

			GridBagConstraints c = new GridBagConstraints();
			c.fill = GridBagConstraints.NONE;
			c.weightx = 0.0;
			c.weighty = 0.0;
			c.anchor = GridBagConstraints.CENTER;
			c.insets = new Insets(1, 2, 0, 2);
			c.gridx = GridBagConstraints.RELATIVE;
			c.weightx = 1.0;
			c.fill = GridBagConstraints.HORIZONTAL;
            add(contents, c);
			c.fill = GridBagConstraints.NONE;
			c.weightx = 0.0;
			add(addButton, c);
			add(removeButton, c);
		}
    };

	AbstractAction removeAction = new AbstractAction("-") {
		  public void actionPerformed(ActionEvent ae) {
		  }
	  };


    private ArrayList panels = new ArrayList();


}
