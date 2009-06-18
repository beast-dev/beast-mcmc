/*
 * PriorsPanel.java
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

package dr.app.beauti.treespanel;

import dr.app.beauti.BeautiFrame;
import dr.app.beauti.util.PanelUtils;
import dr.app.beauti.options.*;
import dr.app.tools.TemporalRooting;
import dr.evolution.alignment.Patterns;
import dr.evolution.distance.DistanceMatrix;
import dr.evolution.distance.F84DistanceMatrix;
import dr.evolution.tree.NeighborJoiningTree;
import dr.evolution.tree.Tree;
import dr.evolution.tree.UPGMATree;
import org.virion.jam.components.WholeNumberField;
import org.virion.jam.components.RealNumberField;
import org.virion.jam.panels.OptionsPanel;
import org.virion.jam.table.HeaderRenderer;
import org.virion.jam.table.TableEditorStopper;

import javax.swing.*;
import javax.swing.event.ListSelectionEvent;
import javax.swing.event.ListSelectionListener;
import javax.swing.plaf.BorderUIResource;
import javax.swing.table.AbstractTableModel;
import javax.swing.table.TableColumnModel;
import java.awt.*;
import java.awt.event.*;

/**
 * @author Andrew Rambaut
 * @author Alexei Drummond
 * @version $Id: PriorsPanel.java,v 1.9 2006/09/05 13:29:34 rambaut Exp $
 */
public class PartitionTreePanel extends OptionsPanel {

	private JComboBox startingTreeCombo = new JComboBox(StartingTreeType.values());
	private JComboBox userTreeCombo = new JComboBox();

	private BeautiFrame frame = null;
	private BeautiOptions options = null;

    private GenerateTreeAction generateTreeAction = new GenerateTreeAction();

	private GenerateTreeDialog generateTreeDialog = null;

    private final PartitionTree partitionTree;

    public PartitionTreePanel(PartitionTree partitionTree) {

		this.partitionTree = partitionTree;

		PanelUtils.setupComponent(startingTreeCombo);
		startingTreeCombo.addItemListener(
				new ItemListener() {
					public void itemStateChanged(ItemEvent ev) {
						fireTreePriorsChanged();
						setupPanel();
					}
				}
		);

		PanelUtils.setupComponent(userTreeCombo);
		userTreeCombo.addItemListener(
				new ItemListener() {
					public void itemStateChanged(ItemEvent ev) {
						fireTreePriorsChanged();
					}
				}
		);


		setOpaque(false);
		setLayout(new BorderLayout(0, 0));
		setBorder(new BorderUIResource.EmptyBorderUIResource(new Insets(12, 12, 12, 12)));

		setupPanel();
	}

	private void fireTreePriorsChanged() {
		if (!settingOptions) {
            frame.setDirty();
		}
	}

	private void generateTree() {
		if (generateTreeDialog == null) {
			generateTreeDialog = new GenerateTreeDialog(frame);
		}

		int result = generateTreeDialog.showDialog(options);
		if (result != JOptionPane.CANCEL_OPTION) {
			GenerateTreeDialog.MethodTypes methodType = generateTreeDialog.getMethodType();
			PartitionData partition = generateTreeDialog.getDataPartition();

			Patterns patterns = new Patterns(partition.getAlignment());
			DistanceMatrix distances = new F84DistanceMatrix(patterns);
			Tree tree;
			TemporalRooting temporalRooting;

			switch (methodType) {
				case NJ:
					tree = new NeighborJoiningTree(distances);
					temporalRooting = new TemporalRooting(tree);
					tree = temporalRooting.findRoot(tree);
					break;
				case UPGMA:
					tree = new UPGMATree(distances);
					temporalRooting = new TemporalRooting(tree);
					break;
				default:
					throw new IllegalArgumentException("unknown method type");
			}

			tree.setId(generateTreeDialog.getName());
			options.userTrees.add(tree);
		}

		fireTreePriorsChanged();

	}

	private void setupPanel() {

        removeAll();

		JPanel panel = new JPanel(new FlowLayout(FlowLayout.LEFT));
		panel.setOpaque(false);
		panel.add(startingTreeCombo);
		if (startingTreeCombo.getSelectedItem() == StartingTreeType.USER) {
			panel.add(new JLabel("  Select Tree:"));
			panel.add(userTreeCombo);
		}
		addComponentWithLabel("                          Starting Tree:", panel);

		addSeparator();

//		generateTreeAction.setEnabled(options != null && options.dataPartitions.size() > 0);

		validate();
		repaint();
	}

	private boolean settingOptions = false;

	public void setOptions(BeautiOptions options) {
		this.options = options;

		settingOptions = true;

		startingTreeCombo.setSelectedItem(options.startingTreeType);

		userTreeCombo.removeAllItems();
		if (options.userTrees.size() == 0) {
			userTreeCombo.addItem("no trees loaded");
			userTreeCombo.setEnabled(false);
		} else {
			for (Tree tree : options.userTrees) {
				userTreeCombo.addItem(tree.getId());
			}
			userTreeCombo.setEnabled(true);
		}

		setupPanel();

		settingOptions = false;

		validate();
		repaint();
	}

	public void getOptions(BeautiOptions options) {

		options.startingTreeType = (StartingTreeType) startingTreeCombo.getSelectedItem();
		options.userStartingTree = getSelectedUserTree();
	}

	private Tree getSelectedUserTree() {
		String treeId = (String) userTreeCombo.getSelectedItem();
		for (Tree tree : options.userTrees) {
			if (tree.getId().equals(treeId)) {
				return tree;
			}
		}
		return null;
	}

	public class GenerateTreeAction extends AbstractAction {
		public GenerateTreeAction() {
			super("Create Tree");
			setToolTipText("Create a NJ or UPGMA tree using a data partition");
		}

        public void actionPerformed(ActionEvent ae) {
			generateTree();
		}
	}


}