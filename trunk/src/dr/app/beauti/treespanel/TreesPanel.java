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

import dr.app.beauti.util.PanelUtils;
import dr.app.beauti.components.SequenceErrorModelComponentOptions;
import dr.app.beauti.*;
import dr.app.beauti.options.*;
import dr.app.tools.TemporalRooting;
import dr.evolution.alignment.Patterns;
import dr.evolution.distance.DistanceMatrix;
import dr.evolution.distance.F84DistanceMatrix;
import dr.evolution.tree.NeighborJoiningTree;
import dr.evolution.tree.Tree;
import dr.evolution.tree.UPGMATree;

import org.virion.jam.framework.Exportable;
import org.virion.jam.panels.ActionPanel;
import org.virion.jam.table.HeaderRenderer;
import org.virion.jam.table.TableEditorStopper;
import org.virion.jam.table.TableRenderer;

import javax.swing.*;
import javax.swing.border.TitledBorder;
import javax.swing.event.ListSelectionEvent;
import javax.swing.event.ListSelectionListener;
import javax.swing.plaf.BorderUIResource;
import javax.swing.table.AbstractTableModel;
import javax.swing.table.TableColumn;
import javax.swing.table.TableColumnModel;
import java.awt.*;
import java.awt.event.*;
import java.util.EnumSet;
import java.util.HashMap;
import java.util.Map;
import java.util.Set;


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
    private TreeDisplayPanel treeDisplayPanel;

    private BeautiFrame frame = null;
    private BeautiOptions options = null;

    private JScrollPane scrollPane = new JScrollPane();
    private JTable treesTable = null;
    private TreesTableModel treesTableModel = null;

    private GenerateTreeDialog generateTreeDialog = null;    
    private boolean settingOptions = false;
//    boolean hasAlignment = false;        
    
    JPanel treeModelPanelParent;
//    private OptionsPanel currentTreeModel = new OptionsPanel();
    PartitionTreeModel currentTreeModel = null;
    TitledBorder treeModelBorder;
    Map<PartitionTreeModel, PartitionTreeModelPanel> treeModelPanels = new HashMap<PartitionTreeModel, PartitionTreeModelPanel>();
    //  private OptionsPanel treePriorPanel = new OptionsPanel();
    JPanel treePriorPanelParent;
    PartitionTreePrior currentTreePrior = null;
    TitledBorder treePriorBorder;
    Map<PartitionTreePrior, PartitionTreePriorPanel> treePriorPanels = new HashMap<PartitionTreePrior, PartitionTreePriorPanel>();

    
    TitledBorder treeBorder;

    // Overall model parameters ////////////////////////////////////////////////////////////////////////
    private boolean isCheckedTipDate = false;
    

    ////////////////////////////////////////////////////////////////////////////////////////////////////

	SequenceErrorModelComponentOptions comp;

    public TreesPanel(BeautiFrame parent, Action removeTreeAction) {
    	super();
        this.frame = parent;

        treesTableModel = new TreesTableModel();
        treesTable = new JTable(treesTableModel);

        treesTable.getTableHeader().setReorderingAllowed(false);
        treesTable.getTableHeader().setResizingAllowed(false);
        treesTable.getTableHeader().setDefaultRenderer(
                new HeaderRenderer(SwingConstants.LEFT, new Insets(0, 4, 0, 4)));

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

        scrollPane = new JScrollPane(treesTable,
                JScrollPane.VERTICAL_SCROLLBAR_ALWAYS,
                JScrollPane.HORIZONTAL_SCROLLBAR_NEVER);
        scrollPane.setOpaque(false);

//        JToolBar toolBar1 = new JToolBar();
//        toolBar1.setFloatable(false);
//        toolBar1.setOpaque(false);
//        toolBar1.setLayout(new FlowLayout(FlowLayout.LEFT, 0, 0));
//        button = new JButton(createTreeAction);
//        createTreeAction.setEnabled(true);
//        PanelUtils.setupComponent(button);
//        toolBar1.add(button);

//        button = new JButton(linkModelAction);
//        linkModelAction.setEnabled(false);
//        PanelUtils.setupComponent(button);
//        toolBar1.add(button);

        ActionPanel actionPanel1 = new ActionPanel(false);
        actionPanel1.setAddAction(addTreeAction);
        actionPanel1.setRemoveAction(removeTreeAction);

        JPanel controlPanel1 = new JPanel(new FlowLayout(FlowLayout.LEFT));
        controlPanel1.setOpaque(false);
        controlPanel1.add(actionPanel1);
    
        setCurrentModel(null);
        setCurrentPrior(null);

        JPanel panel1 = new JPanel(new BorderLayout(0, 0));
        panel1.setOpaque(false);
        panel1.add(scrollPane, BorderLayout.CENTER);
        panel1.add(controlPanel1, BorderLayout.SOUTH);
        
        JPanel panel2 = new JPanel(new BorderLayout(0, 0));        
        panel2.setOpaque(false);
        
        treeModelPanelParent = new JPanel(new FlowLayout(FlowLayout.LEFT));
        treeModelPanelParent.setOpaque(false);
        treeModelBorder = new TitledBorder("Tree Model");
        treeModelPanelParent.setBorder(treeModelBorder);
//        currentTreeModel.setBorder(null);
               
        panel2.add(treeModelPanelParent, BorderLayout.NORTH);
        
        treeDisplayPanel = new TreeDisplayPanel(parent);
        panel2.add(treeDisplayPanel, BorderLayout.CENTER);        
        
        JSplitPane splitPane = new JSplitPane(JSplitPane.HORIZONTAL_SPLIT, panel1, panel2);
        splitPane.setDividerLocation(180);
        splitPane.setContinuousLayout(true);
        splitPane.setBorder(BorderFactory.createEmptyBorder());
        splitPane.setOpaque(false);
                       
        setOpaque(false);
        setLayout(new BorderLayout(0, 0));
        setBorder(new BorderUIResource.EmptyBorderUIResource(new Insets(12, 12, 12, 12)));
        
        treePriorPanelParent = new JPanel(new FlowLayout(FlowLayout.LEFT));
        treePriorPanelParent.setOpaque(false);
        treePriorBorder = new TitledBorder("Tree Prior");
        treePriorPanelParent.setBorder(treePriorBorder);
        
        add(treePriorPanelParent, BorderLayout.SOUTH);
        add(splitPane, BorderLayout.CENTER);

        comp = new SequenceErrorModelComponentOptions ();
    }

    private void fireTreePriorsChanged() {
        if (!settingOptions) {
            frame.setDirty();
        }
    }
    public void removeSelection() {
        int selRow = treesTable.getSelectedRow();
        if (!isUsed(selRow)) {
            PartitionTreeModel model = options.getPartitionTreeModels().get(selRow);
            options.getPartitionTreeModels().remove(model);
        }

        treesTableModel.fireTableDataChanged();
        int n = options.getPartitionTreeModels().size();
        if (selRow >= n) {
            selRow--;
        }
        treesTable.getSelectionModel().setSelectionInterval(selRow, selRow);
        if (n == 0) {
            setCurrentModel(null);
        }

        fireTreePriorsChanged();
    }
        
    private void selectionChanged() {
        int selRow = treesTable.getSelectedRow();
        if (selRow >= 0) {
        	setCurrentModel(options.getPartitionTreeModels().get(selRow));
        	setCurrentPrior(options.getPartitionTreePriors().get(selRow));
//TODO            treeDisplayPanel.setTree(options.userTrees.get(selRow));            
            frame.modelSelectionChanged(!isUsed(selRow));
        } else {
        	setCurrentModel(null);
        	setCurrentPrior(null);
            treeDisplayPanel.setTree(null);
        }
    }

    private void createTree() {
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
            treesTableModel.fireTableDataChanged();
            int row = options.userTrees.size() - 1;
            treesTable.getSelectionModel().setSelectionInterval(row, row);
        }

        fireTreePriorsChanged();

    }

    
    /**
     * Sets the current model that this model panel is displaying
     *
     * @param model the new model to display
     */
    private void setCurrentModel(PartitionTreeModel model) {

        if (model != null) {
            if (currentTreeModel != null) treeModelPanelParent.removeAll();

            PartitionTreeModelPanel panel = treeModelPanels.get(model);
            if (panel == null) {
                panel = new PartitionTreeModelPanel(model);
                treeModelPanels.put(model, panel);
            }

            currentTreeModel = model;
            treeModelPanelParent.add(panel);
            
            repaint();
        }
    }
    
    private void setCurrentPrior(PartitionTreePrior prior) {

        if (prior != null) {
            if (currentTreePrior != null) treePriorPanelParent.removeAll();

            PartitionTreePriorPanel panel = treePriorPanels.get(prior);
            if (panel == null) {
                panel = new PartitionTreePriorPanel(prior);
                treePriorPanels.put(prior, panel);
            }

            currentTreePrior = prior;
            treePriorPanelParent.add(panel);

            repaint();
        }
    }
    
    public void setCheckedTipDate(boolean isCheckedTipDate) {
		this.isCheckedTipDate = isCheckedTipDate;
	}
    
    public void setOptions(BeautiOptions options) {
        this.options = options;

        settingOptions = true;

        Set<PartitionTreeModel> models = treeModelPanels.keySet();        
        
        for (PartitionTreeModel model : models) {
        	treeModelPanels.get(model).setOptions(options);
     	}
        
        Set<PartitionTreePrior> priors = treePriorPanels.keySet();
        
        for (PartitionTreePrior prior : priors) {
        	PartitionTreePriorPanel ptpp = treePriorPanels.get(prior);
        	
        	ptpp.setOptions(options);        	
        	
			if (isCheckedTipDate) { 
        		ptpp.treePriorCombo.removeItem(TreePrior.YULE);
        		ptpp.treePriorCombo.removeItem(TreePrior.BIRTH_DEATH);
        	} else {
        		ptpp.treePriorCombo = new JComboBox(EnumSet.range(TreePrior.CONSTANT, TreePrior.BIRTH_DEATH).toArray());    		
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

        validate();
        repaint();
    }

    public void getOptions(BeautiOptions options) {
    	Set<PartitionTreeModel> models = treeModelPanels.keySet();        
        
        for (PartitionTreeModel model : models) {
        	treeModelPanels.get(model).setOptions(options);
     	}
        
        Set<PartitionTreePrior> priors = treePriorPanels.keySet();
        
        for (PartitionTreePrior prior : priors) {
        	treePriorPanels.get(prior).getOptions(options);
     	}
    }
    
 
    
    public JComponent getExportableComponent() {
        return treeDisplayPanel;
    }
    
    private boolean isUsed(int row) {
        PartitionTreeModel model = options.getPartitionTreeModels().get(row);
        for (PartitionData partition : options.dataPartitions) {
            if (partition.getPartitionTreeModel() == model) {
                return true;
            }
        }
        return false;
    }

    class TreesTableModel extends AbstractTableModel {

        private static final long serialVersionUID = -6707994233020715574L;
        String[] columnNames = {"Tree(s)"};

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
    
    Action addTreeAction = new AbstractAction("+") {
        public void actionPerformed(ActionEvent ae) {
        	createTree();
        }
    };
    
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