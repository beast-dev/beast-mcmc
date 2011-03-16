/*
 * ModelPanel.java
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
import dr.app.beauti.BeautiPanel;
import dr.app.beauti.options.AbstractPartitionData;
import dr.app.beauti.options.BeautiOptions;
import dr.app.beauti.options.PartitionTreeModel;
import dr.app.beauti.options.PartitionTreePrior;
import dr.app.beauti.types.TreePriorType;
import dr.app.gui.table.TableEditorStopper;
import dr.evolution.datatype.Microsatellite;
import jam.framework.Exportable;
import jam.panels.OptionsPanel;
import jam.table.TableRenderer;

import javax.swing.*;
import javax.swing.border.TitledBorder;
import javax.swing.event.ListSelectionEvent;
import javax.swing.event.ListSelectionListener;
import javax.swing.plaf.BorderUIResource;
import javax.swing.table.AbstractTableModel;
import javax.swing.table.TableColumn;
import javax.swing.table.TableColumnModel;
import java.awt.*;
import java.awt.event.ItemEvent;
import java.awt.event.ItemListener;
import java.util.HashMap;
import java.util.Map;


/**
 * @author Andrew Rambaut
 * @author Walter Xie
 * @version $Id:$
 */
public class TreesPanel extends BeautiPanel implements Exportable {

    public final static boolean DEBUG = false;

    private static final long serialVersionUID = 2778103564318492601L;

//    private JComboBox userTreeCombo = new JComboBox();
//    private JButton button;

//    private CreateTreeAction createTreeAction = new CreateTreeAction();
//    private TreeDisplayPanel treeDisplayPanel;

    private BeautiFrame frame = null;
    private BeautiOptions options = null;

    private JTable treesTable = null;
    private TreesTableModel treesTableModel = null;

    //    private GenerateTreeDialog generateTreeDialog = null;
    private boolean settingOptions = false;
//    boolean hasAlignment = false;

    public JCheckBox linkTreePriorCheck = new JCheckBox("Link tree prior for all trees");

    JPanel treeModelPanelParent;
    //    private OptionsPanel currentTreeModel = new OptionsPanel();
    public PartitionTreeModel currentTreeModel = null;
    TitledBorder treeModelBorder;
    Map<PartitionTreeModel, PartitionTreeModelPanel> treeModelPanels = new HashMap<PartitionTreeModel, PartitionTreeModelPanel>();
    //  private OptionsPanel treePriorPanel = new OptionsPanel();
    JPanel treePriorPanelParent;
    //    public PartitionTreePrior currentTreePrior = null;
    TitledBorder treePriorBorder;
    Map<PartitionTreePrior, OptionsPanel> treePriorPanels = new HashMap<PartitionTreePrior, OptionsPanel>();

    private boolean isCheckedTipDate = false;

    public TreesPanel(BeautiFrame parent, Action removeTreeAction) {
        super();
        this.frame = parent;

        treesTableModel = new TreesTableModel();
        treesTable = new JTable(treesTableModel);

        treesTable.getTableHeader().setReorderingAllowed(false);
        treesTable.getTableHeader().setResizingAllowed(false);
//        treesTable.getTableHeader().setDefaultRenderer(
//                new HeaderRenderer(SwingConstants.LEFT, new Insets(0, 4, 0, 4)));

        final TableColumnModel model = treesTable.getColumnModel();
        final TableColumn tableColumn0 = model.getColumn(0);
        tableColumn0.setCellRenderer(new ModelsTableCellRenderer(SwingConstants.LEFT, new Insets(0, 4, 0, 4)));

        TableEditorStopper.ensureEditingStopWhenTableLosesFocus(treesTable);

        treesTable.getSelectionModel().setSelectionMode(ListSelectionModel.SINGLE_SELECTION);
        treesTable.getSelectionModel().addListSelectionListener(new ListSelectionListener() {
            public void valueChanged(ListSelectionEvent evt) {
                selectionChanged();
            }
        });

        JScrollPane scrollPane = new JScrollPane(treesTable,
                JScrollPane.VERTICAL_SCROLLBAR_ALWAYS,
                JScrollPane.HORIZONTAL_SCROLLBAR_NEVER);
        scrollPane.setOpaque(false);

        JPanel controlPanel1 = new JPanel(new FlowLayout(FlowLayout.LEFT));
        controlPanel1.setOpaque(false);

        setCurrentModelAndPrior(null);

        JPanel panel1 = new JPanel(new BorderLayout(0, 0));
        panel1.setOpaque(false);
        panel1.add(scrollPane, BorderLayout.CENTER);

        treeModelPanelParent = new JPanel(new FlowLayout(FlowLayout.LEFT));
        treeModelPanelParent.setOpaque(false);
        treeModelBorder = new javax.swing.border.TitledBorder("Tree Model");
        treeModelPanelParent.setBorder(treeModelBorder);

        JScrollPane scrollPane2 = new JScrollPane(treeModelPanelParent, JScrollPane.VERTICAL_SCROLLBAR_AS_NEEDED, JScrollPane.HORIZONTAL_SCROLLBAR_NEVER);
        scrollPane2.setOpaque(false);
        scrollPane2.setBorder(null);
        scrollPane2.getViewport().setOpaque(false);

        JSplitPane splitPane = new JSplitPane(JSplitPane.HORIZONTAL_SPLIT, panel1, scrollPane2);
        splitPane.setDividerLocation(180);
        splitPane.setContinuousLayout(true);
        splitPane.setBorder(BorderFactory.createEmptyBorder());
        splitPane.setOpaque(false);

        treePriorPanelParent = new JPanel(new FlowLayout(FlowLayout.LEFT));
        treePriorPanelParent.setOpaque(false);
        treePriorBorder = new TitledBorder("Tree Prior");
        treePriorPanelParent.setBorder(treePriorBorder);

        JPanel panel3 = new JPanel(new FlowLayout(FlowLayout.LEFT));
        panel3.setOpaque(false);
        linkTreePriorCheck.setEnabled(true);
        linkTreePriorCheck.setSelected(true);
        linkTreePriorCheck.addItemListener(new ItemListener() {
            public void itemStateChanged(ItemEvent ev) {
                updateShareSameTreePriorChanged();
            }
        });
        linkTreePriorCheck.setToolTipText("Decide whether to use one tree prior for all trees");
        panel3.add(linkTreePriorCheck);

        JPanel panel4 = new JPanel(new BorderLayout());
        panel4.setOpaque(false);
        panel4.add(panel3, BorderLayout.NORTH);
        panel4.add(treePriorPanelParent, BorderLayout.CENTER);

        setOpaque(false);
        setLayout(new BorderLayout(0, 0));
        setBorder(new BorderUIResource.EmptyBorderUIResource(new Insets(12, 12, 12, 12)));

        add(panel4, BorderLayout.NORTH);
        add(splitPane, BorderLayout.CENTER);
    }

    public void updateShareSameTreePriorChanged() {
//    	options.shareSameTreePrior = linkTreePriorCheck.isSelected();
//    	fireShareSameTreePriorChanged ();

        if (linkTreePriorCheck.isSelected()) {
            options.linkTreePriors(currentTreeModel.getPartitionTreePrior());
        } else {
            options.unLinkTreePriors();
//    		currentTreePrior = currentTreeModel.getPartitionTreePrior();
        }
//    	setCurrentModelAndPrior(currentTreeModel);
        updateTreePriorBorder();
    }

//    private void fireShareSameTreePriorChanged() {
//    	linkTreePriorCheck.setSelected(options.shareSameTreePrior);
//    	if (options.shareSameTreePrior) {
//    		options.activedSameTreePrior = currentTreePrior;
//    		// keep previous prior for reuse
//    	} else {
//    		// reuse previous prior
//    		setCurrentModelAndPrior(currentTreeModel);
//    	}
//    	updateTreePriorBorder();
//    }

    private void updateTreePriorBorder() {
        if (options.useStarBEAST) {
            treePriorBorder.setTitle("Species tree prior used to start all gene tree models");
        } else if (options.isShareSameTreePrior()) {
            treePriorBorder.setTitle("Tree prior shared by all tree models");
        } else {
            treePriorBorder.setTitle("Tree Prior - " + currentTreeModel.getPartitionTreePrior().getName());
        }
        repaint();
    }

    private void fireTreePriorsChanged() {
        options.updatePartitionAllLinks();
        frame.setDirty();
    }

    public void updatePriorPanelForSpeciesAnalysis() {
        linkTreePriorCheck.setSelected(true);
        updateShareSameTreePriorChanged();

        if (currentTreeModel.getPartitionTreePrior() != null) treePriorPanelParent.removeAll();

        if (options.getPartitionTreePriors().size() > 0) {
            OptionsPanel p;
            if (options.useStarBEAST) {
                linkTreePriorCheck.setEnabled(false);

                options.getPartitionTreePriors().get(0).setNodeHeightPrior(TreePriorType.SPECIES_YULE);

//                options.clockModelOptions.fixRateOfFirstClockPartition(); // fix 1st partition

                p = new SpeciesTreesPanel(options.getPartitionTreePriors().get(0));

            } else {
                if (options.hasData() && options.contains(Microsatellite.INSTANCE)) {
                    linkTreePriorCheck.setEnabled(false);
                } else {
                    linkTreePriorCheck.setEnabled(true);
                }

                options.getPartitionTreePriors().get(0).setNodeHeightPrior(TreePriorType.CONSTANT);

                p = new PartitionTreePriorPanel(options.getPartitionTreePriors().get(0), this);

            }

            treePriorPanels.put(options.getPartitionTreePriors().get(0), p);
            updateTreePriorBorder();
            treePriorPanelParent.add(p);
        }
        repaint();
    }

    private void selectionChanged() {
        int selRow = treesTable.getSelectedRow();
        if (selRow >= 0) {
            PartitionTreeModel ptm = options.getPartitionTreeModels().get(selRow);
            setCurrentModelAndPrior(ptm);
            //TODO            treeDisplayPanel.setTree(options.userTrees.get(selRow));
        }
    }

    /**
     * Sets the current model that this model panel is displaying
     *
     * @param model the new model to display
     */
    private void setCurrentModelAndPrior(PartitionTreeModel model) {

        if (model != null) {
            if (currentTreeModel != null) treeModelPanelParent.removeAll();

            PartitionTreeModelPanel panel = treeModelPanels.get(model);
            if (panel == null) {
                panel = new PartitionTreeModelPanel(frame, model, options);
                treeModelPanels.put(model, panel);
            }

            currentTreeModel = model;
            treeModelBorder.setTitle("Tree Model - " + model.getName());
            treeModelPanelParent.add(panel);

            // ++++++++++++++ PartitionTreePrior ++++++++++++++++
            if (model.getPartitionTreePrior() != null) treePriorPanelParent.removeAll();

            PartitionTreePrior prior = model.getPartitionTreePrior();

            OptionsPanel panel1 = treePriorPanels.get(prior);
            if (panel1 == null) {
                panel1 = new PartitionTreePriorPanel(prior, this);
                treePriorPanels.put(prior, panel1);
            }

//            currentTreePrior = prior;
            updateTreePriorBorder();
            treePriorPanelParent.add(panel1);

            repaint();
        } else {
            //TODO
        }
    }

    private void resetPanel() {
        if (!options.hasData()) {
            currentTreeModel = null;
            treeModelPanels.clear();
//            currentTreePrior = null;
            treePriorPanels.clear();

            treeModelPanelParent.removeAll();
            treeModelBorder.setTitle("Tree Model");
            treePriorPanelParent.removeAll();
            treePriorBorder.setTitle("Tree Prior");

            return;
        }
    }

    public void setOptions(BeautiOptions options) {
        this.options = options;

        resetPanel();

        settingOptions = true;

        linkTreePriorCheck.setEnabled(options.getPartitionTreeModels().size() > 1 && (!options.contains(Microsatellite.INSTANCE)) );
        linkTreePriorCheck.setSelected(options.isShareSameTreePrior()); // important

        for (PartitionTreeModel model : treeModelPanels.keySet()) {
            if (model != null) {
                treeModelPanels.get(model).setOptions();
            }
        }

        for (PartitionTreePrior prior : options.getPartitionTreePriors()) {
            if (options.useStarBEAST) {
                SpeciesTreesPanel ptpp = (SpeciesTreesPanel) treePriorPanels.get(prior);
                if (ptpp != null) {
                    ptpp.setOptions();
                }
            } else {
                PartitionTreePriorPanel ptpp = (PartitionTreePriorPanel) treePriorPanels.get(prior);
                if (ptpp != null) {
                    ptpp.setOptions();

                    if (options.contains(Microsatellite.INSTANCE)) {
                        ptpp.setMicrosatelliteTreePrior();
                    } else if (options.clockModelOptions.isTipCalibrated()) {
                        ptpp.removeCertainPriorFromTreePriorCombo();
                    } else {
                        ptpp.recoveryTreePriorCombo();
                    }
                    ptpp.repaint();
                }
            }
        }

        settingOptions = false;

        int selRow = treesTable.getSelectedRow();
        treesTableModel.fireTableDataChanged();
        if (options.getPartitionTreeModels().size() > 0) {
            if (selRow < 0) {
                selRow = 0;
            }
            treesTable.getSelectionModel().setSelectionInterval(selRow, selRow);
        }

        if (currentTreeModel == null && options.getPartitionTreeModels().size() > 0) {
            treesTable.getSelectionModel().setSelectionInterval(0, 0);
        }

//        fireShareSameTreePriorChanged();

        validate();
        repaint();
    }

    public void getOptions(BeautiOptions options) {
        if (settingOptions) return;

//    	Set<PartitionTreeModel> models = treeModelPanels.keySet();
//
//        for (PartitionTreeModel model : models) {
//        	treeModelPanels.get(model).getOptions(options);
//     	}

        for (PartitionTreePrior prior : options.getPartitionTreePriors()) {
            if (options.useStarBEAST) {
                SpeciesTreesPanel ptpp = (SpeciesTreesPanel) treePriorPanels.get(prior);
                if (ptpp != null) {
                    ptpp.getOptions();
                }
            } else {
                PartitionTreePriorPanel ptpp = (PartitionTreePriorPanel) treePriorPanels.get(prior);
                if (ptpp != null) {
                    ptpp.getOptions();
                }
            }
        }
    }

    public BeautiFrame getFrame() {
        return frame;
    }

    public JComponent getExportableComponent() {
//        return treeDisplayPanel;
        return this;
    }

    private boolean isUsed(int row) {
        PartitionTreeModel model = options.getPartitionTreeModels().get(row);
        for (AbstractPartitionData partition : options.dataPartitions) {
            if (partition.getPartitionTreeModel() == model) {
                return true;
            }
        }
        return false;
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
            return options.getPartitionTreeModels().size();
        }

        public Object getValueAt(int row, int col) {
            PartitionTreeModel model = options.getPartitionTreeModels().get(row);
            switch (col) {
                case 0:
                    return model.getName();
                default:
                    throw new IllegalArgumentException("unknown column, " + col);
            }
        }

        public void setValueAt(Object aValue, int row, int col) {
//            Tree tree = options.userTrees.get(row);
            switch (col) {
                case 0:
                    String name = ((String) aValue).trim();
                    if (name.length() > 0) {
                        PartitionTreeModel model = options.getPartitionTreeModels().get(row);
                        model.setName(name);
                        // keep tree prior name same as tree model name
                        PartitionTreePrior prior = model.getPartitionTreePrior();
                        prior.setName(name);
                        fireTreePriorsChanged();
                    }
                    break;
                default:
                    throw new IllegalArgumentException("unknown column, " + col);
            }
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

    class ModelsTableCellRenderer extends TableRenderer {

        public ModelsTableCellRenderer(int alignment, Insets insets) {
            super(alignment, insets);
        }

        public Component getTableCellRendererComponent(JTable aTable,
                                                       Object value,
                                                       boolean aIsSelected,
                                                       boolean aHasFocus,
                                                       int aRow, int aColumn) {

            if (value == null) return this;

            Component renderer = super.getTableCellRendererComponent(aTable,
                    value,
                    aIsSelected,
                    aHasFocus,
                    aRow, aColumn);

            if (!isUsed(aRow))
                renderer.setForeground(Color.gray);
            else
                renderer.setForeground(Color.black);
            return this;
        }

    }

//    Action addTreeAction = new AbstractAction("+") {
//        public void actionPerformed(ActionEvent ae) {
//        	createTree();
//        }
//    };

//    public class CreateTreeAction extends AbstractAction {
//        public CreateTreeAction() {
//            super("Create Tree");
//            setToolTipText("Create a NJ or UPGMA tree using a data partition");
//        }
//
//        public void actionPerformed(ActionEvent ae) {
//            createTree();
//        }
//    }
}