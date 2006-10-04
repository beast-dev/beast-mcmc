/*
 * RulesPanel.java
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
 * @version $Id: RulesPanel.java,v 1.2 2006/09/09 18:16:16 rambaut Exp $
 */


public class RulesPanel extends JPanel {

    private Icon addIcon = null;
	private Icon addRolloverIcon = null;
	private Icon addPressedIcon = null;
    private Icon removeIcon = null;
	private Icon removeRolloverIcon = null;
	private Icon removePressedIcon = null;

    public RulesPanel(RuleModel ruleModel) {
        this.ruleModel = ruleModel;

        try {
            addIcon = IconUtils.getIcon(RulesPanel.class, "images/plusminus/plus.png");
	        addRolloverIcon = IconUtils.getIcon(RulesPanel.class, "images/plusminus/plusRollover.png");
	        addPressedIcon = IconUtils.getIcon(RulesPanel.class, "images/plusminus/plusPressed.png");
            removeIcon = IconUtils.getIcon(RulesPanel.class, "images/plusminus/minus.png");
	        removeRolloverIcon = IconUtils.getIcon(RulesPanel.class, "images/plusminus/minusRollover.png");
	        removePressedIcon = IconUtils.getIcon(RulesPanel.class, "images/plusminus/minusPressed.png");
        } catch (Exception e) {
            // do nothing
        }

        BoxLayout layout = new BoxLayout(this, BoxLayout.PAGE_AXIS);
        setLayout(layout);
        setBorder(BorderFactory.createEmptyBorder(4, 4, 4, 4));
        setOpaque(false);
        addRule(null);
        Dimension dim = super.getPreferredSize();
        dim.height += getComponent(0).getPreferredSize().getHeight();
        setMinimumSize(dim);
        setPreferredSize(dim);
    }

	public List getRules() {
		return Collections.unmodifiableList(rules);
	}

    private void addRule(Rule previousRule) {

        final DefaultRule rule = new DefaultRule();
	    RulePanel rulePanel = new RulePanel(rule);

        if (previousRule != null) {
            int index = rules.indexOf(previousRule);
            add(rulePanel, index + 1);
            rules.add(index + 1, rule);
        } else {
            add(rulePanel, 0);
            rules.add(0, rule);
        }
	    removeAction.setEnabled(rules.size() > 1);
        validate();
	    repaint();
    }

    private void removeRule(Rule rule) {
        int index = rules.indexOf(rule);
        remove(index);
        rules.remove(index);
	    removeAction.setEnabled(rules.size() > 1);
        validate();
	    repaint();
    }

    public interface Rule {
        Object getField();

        Object getCondition();

        Object getValue();
    };

    class DefaultRule implements Rule {
        JComboBox fieldCombo;
        JComboBox conditionCombo;
        JTextField valueText;

        public Object getField() {
            return fieldCombo.getSelectedItem();
        }

        public Object getCondition() {
            return conditionCombo.getSelectedItem();
        }

        public Object getValue() {
            return valueText.getText();
        }
    };

    class RulePanel extends JPanel {
		RulePanel(final DefaultRule rule) {

			setLayout(new GridBagLayout());

			setOpaque(true);
			setBorder(BorderFactory.createMatteBorder(1, 1, 1, 1, Color.lightGray));
			setBackground(new Color(0.0F, 0.0F, 0.0F, 0.05F));

			rule.fieldCombo = new JComboBox(ruleModel.getFields());
			rule.conditionCombo = new JComboBox(ruleModel.getConditions(rule.getField()));

			rule.fieldCombo.addItemListener(new ItemListener() {
			    public void itemStateChanged(ItemEvent ie) {
			        rule.conditionCombo.setModel(new DefaultComboBoxModel(ruleModel.getConditions(rule.getField())));
			    }
			});
			rule.valueText = new JTextField("");
			rule.valueText.setColumns(12);

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
			        addRule(rule);
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
			        removeRule(rule);
			    }
			});

			GridBagConstraints c = new GridBagConstraints();
			c.fill = GridBagConstraints.NONE;
			c.weightx = 0.0;
			c.weighty = 0.0;
			c.anchor = GridBagConstraints.CENTER;
			c.insets = new Insets(1, 2, 0, 2);
			c.gridx = GridBagConstraints.RELATIVE;
			add(rule.fieldCombo, c);
			add(rule.conditionCombo, c);
			c.weightx = 1.0;
			c.fill = GridBagConstraints.HORIZONTAL;
			add(rule.valueText, c);
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


    private RuleModel ruleModel;
    private ArrayList rules = new ArrayList();


}
