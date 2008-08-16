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
import dr.app.beauti.PanelUtils;
import dr.app.beauti.options.*;
import dr.evolution.tree.Tree;
import org.virion.jam.components.WholeNumberField;
import org.virion.jam.panels.OptionsPanel;
import org.virion.jam.panels.ActionPanel;
import org.virion.jam.table.HeaderRenderer;
import org.virion.jam.table.TableEditorStopper;

import javax.swing.*;
import javax.swing.event.ListSelectionListener;
import javax.swing.event.ListSelectionEvent;
import javax.swing.table.*;
import javax.swing.plaf.BorderUIResource;
import java.awt.*;
import java.awt.event.*;

/**
 * @author Andrew Rambaut
 * @author Alexei Drummond
 * @version $Id: PriorsPanel.java,v 1.9 2006/09/05 13:29:34 rambaut Exp $
 */
public class TreesPanel extends JPanel {

    private OptionsPanel treePriorPanel = new OptionsPanel();
    private JComboBox treePriorCombo;
    private JComboBox parameterizationCombo = new JComboBox(new String[]{
            "Growth Rate", "Doubling Time"});
    private JComboBox bayesianSkylineCombo = new JComboBox(new String[]{
            "Piecewise-constant", "Piecewise-linear"});
    private WholeNumberField groupCountField = new WholeNumberField(2, Integer.MAX_VALUE);

//    RealNumberField samplingProportionField = new RealNumberField(Double.MIN_VALUE, 1.0);

    private JComboBox startingTreeCombo = new JComboBox(StartingTreeType.values());
    private JComboBox userTreeCombo = new JComboBox();

    private TreeDisplayPanel treeDisplayPanel;

    private BeautiFrame frame = null;
    private BeautiOptions options = null;

    private JScrollPane scrollPane = new JScrollPane();
    private JTable treesTable = null;
    private TreesTableModel treesTableModel = null;

    public TreesPanel(BeautiFrame parent) {

        this.frame = parent;

        treesTableModel = new TreesTableModel();
        treesTable = new JTable(treesTableModel);

        treesTable.getTableHeader().setReorderingAllowed(false);
        treesTable.getTableHeader().setResizingAllowed(false);
        treesTable.getTableHeader().setDefaultRenderer(
                new HeaderRenderer(SwingConstants.LEFT, new Insets(0, 4, 0, 4)));

        final TableColumnModel model = treesTable.getColumnModel();
        final TableColumn tableColumn0 = model.getColumn(0);

        TableEditorStopper.ensureEditingStopWhenTableLosesFocus(treesTable);

        treesTable.getSelectionModel().setSelectionMode(ListSelectionModel.SINGLE_SELECTION);
        treesTable.getSelectionModel().addListSelectionListener(new ListSelectionListener() {
            public void valueChanged(ListSelectionEvent evt) {
                selectionChanged();
            }
        });

        scrollPane = new JScrollPane(treesTable,
                JScrollPane.VERTICAL_SCROLLBAR_ALWAYS,
                JScrollPane.HORIZONTAL_SCROLLBAR_NEVER);
        scrollPane.setOpaque(false);

        ActionPanel actionPanel1 = new ActionPanel(false);
        //actionPanel1.setAddAction(addTreeAction);
        //actionPanel1.setRemoveAction(removeTreeAction);

        JPanel controlPanel1 = new JPanel(new FlowLayout(FlowLayout.LEFT));
        controlPanel1.setOpaque(false);
        controlPanel1.add(actionPanel1);
        java.awt.event.ItemListener listener = new java.awt.event.ItemListener() {
            public void itemStateChanged(java.awt.event.ItemEvent ev) {
                fireTreePriorsChanged();
            }
        };

        treePriorCombo = new JComboBox(TreePrior.values());

        PanelUtils.setupComponent(treePriorCombo);
        treePriorCombo.addItemListener(
                new java.awt.event.ItemListener() {
                    public void itemStateChanged(java.awt.event.ItemEvent ev) {
                        fireTreePriorsChanged();
                        setupPanel();
                    }
                }
        );

        KeyListener keyListener = new KeyAdapter() {
            public void keyTyped(java.awt.event.KeyEvent ev) {
                if (ev.getKeyCode() == KeyEvent.VK_ENTER) {
                    fireTreePriorsChanged();
                }
            }
        };

        groupCountField.addKeyListener(keyListener);
//        samplingProportionField.addKeyListener(keyListener);

        FocusListener focusListener = new FocusAdapter() {
            public void focusLost(FocusEvent focusEvent) {
                fireTreePriorsChanged();
            }
        };
        groupCountField.addFocusListener(focusListener);
//        samplingProportionField.addFocusListener(focusListener);

        PanelUtils.setupComponent(parameterizationCombo);
        parameterizationCombo.addItemListener(listener);

        PanelUtils.setupComponent(bayesianSkylineCombo);
        bayesianSkylineCombo.addItemListener(listener);

        PanelUtils.setupComponent(startingTreeCombo);
        startingTreeCombo.addItemListener(
                new java.awt.event.ItemListener() {
                    public void itemStateChanged(java.awt.event.ItemEvent ev) {
                        fireTreePriorsChanged();
                        setupPanel();
                    }
                }
        );

        PanelUtils.setupComponent(userTreeCombo);
        userTreeCombo.addItemListener(
                new java.awt.event.ItemListener() {
                    public void itemStateChanged(java.awt.event.ItemEvent ev) {
                        fireTreePriorsChanged();
                    }
                }
        );

        JPanel panel1 = new JPanel(new BorderLayout(0, 0));
        panel1.setOpaque(false);
        panel1.add(scrollPane, BorderLayout.CENTER);
        panel1.add(controlPanel1, BorderLayout.SOUTH);

        treeDisplayPanel = new TreeDisplayPanel(parent);

        JSplitPane splitPane = new JSplitPane(JSplitPane.HORIZONTAL_SPLIT, panel1, treeDisplayPanel);
        splitPane.setDividerLocation(180);
        splitPane.setContinuousLayout(true);
        splitPane.setBorder(BorderFactory.createEmptyBorder());
        splitPane.setOpaque(false);

        setOpaque(false);
        setLayout(new BorderLayout(0, 0));
        setBorder(new BorderUIResource.EmptyBorderUIResource(new java.awt.Insets(12, 12, 12, 12)));

        treePriorPanel.setBorder(null);
        add(treePriorPanel, BorderLayout.NORTH);

        add(splitPane, BorderLayout.CENTER);

        setupPanel();
    }

    private void fireTreePriorsChanged() {
        if (!settingOptions) {
            frame.treePriorsChanged();
        }
    }

    private void selectionChanged() {
        int selRow = treesTable.getSelectedRow();
        treeDisplayPanel.setTree(options.trees.get(selRow));
    }

    private void setupPanel() {

        treePriorPanel.removeAll();

        treePriorPanel.addComponentWithLabel("Tree Prior:", treePriorCombo);
        if (treePriorCombo.getSelectedItem() == TreePrior.EXPONENTIAL ||
                treePriorCombo.getSelectedItem() == TreePrior.LOGISTIC ||
                treePriorCombo.getSelectedItem() == TreePrior.EXPANSION) {
            treePriorPanel.addComponentWithLabel("Parameterization for growth:", parameterizationCombo);
        } else if (treePriorCombo.getSelectedItem() == TreePrior.SKYLINE) {
            groupCountField.setColumns(6);
            treePriorPanel.addComponentWithLabel("Number of groups:", groupCountField);
            treePriorPanel.addComponentWithLabel("Skyline Model:", bayesianSkylineCombo);
        } else if (treePriorCombo.getSelectedItem() == TreePrior.BIRTH_DEATH) {
//            samplingProportionField.setColumns(8);
//            treePriorPanel.addComponentWithLabel("Proportion of taxa sampled:", samplingProportionField);
        }

        treePriorPanel.addSeparator();

        JPanel panel = new JPanel(new FlowLayout(FlowLayout.LEFT));
        panel.setOpaque(false);
        panel.add(startingTreeCombo);
        if (startingTreeCombo.getSelectedItem() == StartingTreeType.USER) {
            panel.add(new JLabel("  Select Tree:"));
            panel.add(userTreeCombo);
        }
        treePriorPanel.addComponentWithLabel("                          Starting Tree:", panel);

        treePriorPanel.addSeparator();

        treesTableModel.fireTableDataChanged();

        validate();
        repaint();
    }

    private boolean settingOptions = false;

    public void setOptions(BeautiOptions options) {
        this.options = options;

        settingOptions = true;

        treePriorCombo.setSelectedItem(options.nodeHeightPrior);

        groupCountField.setValue(options.skylineGroupCount);
        //samplingProportionField.setValue(options.birthDeathSamplingProportion);

        parameterizationCombo.setSelectedIndex(options.parameterization);
        bayesianSkylineCombo.setSelectedIndex(options.skylineModel);

        startingTreeCombo.setSelectedItem(options.startingTreeType);

        userTreeCombo.removeAllItems();
        if (options.trees.size() == 0) {
            userTreeCombo.addItem("no trees loaded");
            userTreeCombo.setEnabled(false);
        } else {
            for (Tree tree : options.trees) {
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
        options.nodeHeightPrior = (TreePrior) treePriorCombo.getSelectedItem();

        if (options.nodeHeightPrior == TreePrior.SKYLINE) {
            Integer groupCount = groupCountField.getValue();
            if (groupCount != null) {
                options.skylineGroupCount = groupCount;
            } else {
                options.skylineGroupCount = 5;
            }
        } else if (options.nodeHeightPrior == TreePrior.BIRTH_DEATH) {
//            Double samplingProportion = samplingProportionField.getValue();
//            if (samplingProportion != null) {
//                options.birthDeathSamplingProportion = samplingProportion;
//            } else {
//                options.birthDeathSamplingProportion = 1.0;
//            }
        }

        options.parameterization = parameterizationCombo.getSelectedIndex();
        options.skylineModel = bayesianSkylineCombo.getSelectedIndex();

        options.startingTreeType = (StartingTreeType)startingTreeCombo.getSelectedItem();
        options.userStartingTree = getSelectedUserTree();
    }

    private Tree getSelectedUserTree() {
        String treeId = (String)userTreeCombo.getSelectedItem();
        for (Tree tree : options.trees) {
            if (tree.getId().equals(treeId)) {
                return tree;
            }
        }
        return null;
    }

    class TreesTableModel extends AbstractTableModel {

         private static final long serialVersionUID = -6707994233020715574L;
         String[] columnNames = {"Trees"};

         public TreesTableModel() {
         }

         public int getColumnCount() {
             return columnNames.length;
         }

         public int getRowCount() {
             if (options == null) return 0;
             return options.trees.size();
         }

         public Object getValueAt(int row, int col) {
             Tree tree = options.trees.get(row);
             switch (col) {
                 case 0:
                     return tree.getId();
                 default:
                     throw new IllegalArgumentException("unknown column, " + col);
             }
         }

         public void setValueAt(Object aValue, int row, int col) {
             Tree tree = options.trees.get(row);
             switch (col) {
                 case 0:
                     String name = ((String) aValue).trim();
                     if (name.length() > 0) {
                         tree.setId(name);
                     }
                     break;
             }
             fireTreePriorsChanged();
         }

         public boolean isCellEditable(int row, int col) {
             boolean editable;

             switch (col) {
                 case 0:// name
                     editable = true;
                     break;
                 default:
                     editable = false;
             }

             return editable;
         }


         public String getColumnName(int column) {
             return columnNames[column];
         }

         public Class getColumnClass(int c) {
             if (getRowCount() == 0) {
                 return Object.class;
             }
             return getValueAt(0, c).getClass();
         }

         public String toString() {
             StringBuffer buffer = new StringBuffer();

             buffer.append(getColumnName(0));
             for (int j = 1; j < getColumnCount(); j++) {
                 buffer.append("\t");
                 buffer.append(getColumnName(j));
             }
             buffer.append("\n");

             for (int i = 0; i < getRowCount(); i++) {
                 buffer.append(getValueAt(i, 0));
                 for (int j = 1; j < getColumnCount(); j++) {
                     buffer.append("\t");
                     buffer.append(getValueAt(i, j));
                 }
                 buffer.append("\n");
             }

             return buffer.toString();
         }
     }

     public class UnlinkModelsAction extends AbstractAction {
         public UnlinkModelsAction() {
             super("Unlink Models");
             setToolTipText("Use this tool to use a different model for each selected data partition");
         }

         public void actionPerformed(ActionEvent ae) {
             //unlinkModels();
         }
     }


     public class LinkModelsAction extends AbstractAction {
         public LinkModelsAction() {
             super("Link Models");
             setToolTipText("Use this tool to set all the selected partitions to the same model");
         }

         public void actionPerformed(ActionEvent ae) {
             //linkModels();
         }
     }
}
