/*
 * ModelsPanel.java
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

package dr.app.beauti.siteModelsPanel;

import dr.app.beauti.BeautiFrame;
import dr.app.beauti.BeautiPanel;
import dr.app.beauti.components.sequenceerror.SequenceErrorModelComponentOptions;
import dr.app.beauti.options.AbstractPartitionData;
import dr.app.beauti.options.BeautiOptions;
import dr.app.beauti.options.PartitionSubstitutionModel;
import dr.app.beauti.types.SequenceErrorType;
import dr.app.beauti.util.PanelUtils;
import dr.app.gui.table.TableEditorStopper;
import dr.evolution.datatype.DataType;
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
import java.util.EnumSet;
import java.util.HashMap;
import java.util.Map;
import java.util.logging.Logger;

/**
 * @author Andrew Rambaut
 * @author Alexei Drummond
 * @version $Id: ModelPanel.java,v 1.17 2006/09/05 13:29:34 rambaut Exp $
 */
public class SiteModelsPanel extends BeautiPanel implements Exportable {

    public final static boolean DEBUG = false;

    private static final long serialVersionUID = 2778103564318492601L;

    private static final int MINIMUM_TABLE_WIDTH = 140;

    JTable modelTable = null;
    ModelTableModel modelTableModel = null;
    BeautiOptions options = null;

    JPanel modelPanelParent;
    PartitionSubstitutionModel currentModel = null;
    Map<PartitionSubstitutionModel, PartitionModelPanel> modelPanels = new HashMap<PartitionSubstitutionModel, PartitionModelPanel>();
    TitledBorder modelBorder;

    BeautiFrame frame = null;
//    CreateModelDialog createModelDialog = null;
    boolean settingOptions = false;

    public SiteModelsPanel(BeautiFrame parent, Action removeModelAction) {

        super();

        this.frame = parent;

        modelTableModel = new ModelTableModel();
        modelTable = new JTable(modelTableModel);

        modelTable.getTableHeader().setReorderingAllowed(false);
        modelTable.getTableHeader().setResizingAllowed(false);
//        modelTable.getTableHeader().setDefaultRenderer(
//                new HeaderRenderer(SwingConstants.LEFT, new Insets(0, 4, 0, 4)));

        final TableColumnModel model = modelTable.getColumnModel();
        final TableColumn tableColumn0 = model.getColumn(0);
        tableColumn0.setCellRenderer(new ModelsTableCellRenderer(SwingConstants.LEFT, new Insets(0, 4, 0, 4)));

        TableEditorStopper.ensureEditingStopWhenTableLosesFocus(modelTable);

        modelTable.getSelectionModel().setSelectionMode(ListSelectionModel.SINGLE_SELECTION);
        modelTable.getSelectionModel().addListSelectionListener(new ListSelectionListener() {
            public void valueChanged(ListSelectionEvent evt) {
                selectionChanged();
            }
        });

        JScrollPane scrollPane = new JScrollPane(modelTable,
                JScrollPane.VERTICAL_SCROLLBAR_ALWAYS, JScrollPane.HORIZONTAL_SCROLLBAR_AS_NEEDED);
        scrollPane.setOpaque(false);

//        ActionPanel actionPanel1 = new ActionPanel(false);
//        actionPanel1.setAddAction(addModelAction);
//        actionPanel1.setRemoveAction(removeModelAction);

        JPanel controlPanel1 = new JPanel(new FlowLayout(FlowLayout.LEFT));
        controlPanel1.setOpaque(false);
//        controlPanel1.add(actionPanel1);

        JPanel panel = new JPanel(new BorderLayout(0, 0));
        panel.setOpaque(false);
        panel.add(scrollPane, BorderLayout.CENTER);
//        panel.add(controlPanel1, BorderLayout.SOUTH);
        panel.setMinimumSize(new Dimension(MINIMUM_TABLE_WIDTH, 0));

        modelPanelParent = new JPanel(new FlowLayout(FlowLayout.CENTER));
        modelPanelParent.setOpaque(false);
        modelBorder = new TitledBorder("Substitution Model");
        modelPanelParent.setBorder(modelBorder);

        setCurrentModel(null);

        JScrollPane scrollPane2 = new JScrollPane(modelPanelParent, JScrollPane.VERTICAL_SCROLLBAR_AS_NEEDED, JScrollPane.HORIZONTAL_SCROLLBAR_NEVER);
        scrollPane2.setOpaque(false);
        scrollPane2.setBorder(null);
        scrollPane2.getViewport().setOpaque(false);

        JSplitPane splitPane = new JSplitPane(JSplitPane.HORIZONTAL_SPLIT, panel, scrollPane2);
        splitPane.setDividerLocation(MINIMUM_TABLE_WIDTH);
        splitPane.setContinuousLayout(true);
        splitPane.setBorder(BorderFactory.createEmptyBorder());
        splitPane.setOpaque(false);

        setOpaque(false);
        setBorder(new BorderUIResource.EmptyBorderUIResource(new Insets(12, 12, 12, 12)));
        setLayout(new BorderLayout(0, 0));
        add(splitPane, BorderLayout.CENTER);
    }

    private void resetPanel() {
    	if (!options.hasData()) {
    		currentModel = null;
    		modelPanels.clear();
    		modelPanelParent.removeAll();
    		modelBorder.setTitle("Substitution Model");

//            if (currentDiscreteTraitOption != null) {
//                this.remove(d_splitPane);
//                currentDiscreteTraitOption = null;
//            }
        	return;
        }
    }

    public void setOptions(BeautiOptions options) {

        if (DEBUG) {
            Logger.getLogger("dr.app.beauti").info("ModelsPanel.setOptions");
        }

        this.options = options;

        resetPanel();

        settingOptions = true;

        int selRow = modelTable.getSelectedRow();
        modelTableModel.fireTableDataChanged();
        if (options.getPartitionSubstitutionModels().size() > 0) {
            if (selRow < 0) {
                selRow = 0;
            }
            modelTable.getSelectionModel().setSelectionInterval(selRow, selRow);
        }

        if (currentModel == null && options.getPartitionSubstitutionModels().size() > 0) {
            modelTable.getSelectionModel().setSelectionInterval(0, 0);
        }

        settingOptions = false;

        validate();
        repaint();
    }

    public void getOptions(BeautiOptions options) {
    }

    private void fireModelsChanged() {
        options.updatePartitionAllLinks();
        frame.setDirty();
    }

    private void selectionChanged() {
        int selRow = modelTable.getSelectedRow();

        if (selRow >= options.getPartitionSubstitutionModels().size()) {
            selRow = 0;
            modelTable.getSelectionModel().setSelectionInterval(selRow, selRow);
        }

        if (selRow >= 0) {
            setCurrentModel(options.getPartitionSubstitutionModels().get(selRow));
//            frame.modelSelectionChanged(!isUsed(selRow));
        }
    }

    /**
     * Sets the current model that this model panel is displaying
     *
     * @param model the new model to display
     */
    private void setCurrentModel(PartitionSubstitutionModel model) {
        if (model != null) {
            if (currentModel != null) modelPanelParent.removeAll();

            PartitionModelPanel panel = modelPanels.get(model);
            if (panel == null) {
                panel = new PartitionModelPanel(model);
                modelPanels.put(model, panel);
            }

            currentModel = model;
            panel.setOptions();
            modelPanelParent.add(panel);

            updateBorder();
        }
    }

    private void updateBorder() {

        String title;

        switch (currentModel.getDataType().getType()) {
            case DataType.NUCLEOTIDES:
                title = "Nucleotide";
                break;
            case DataType.AMINO_ACIDS:
                title = "Amino Acid";
                break;
            case DataType.TWO_STATES:
                title = "Binary";
                break;
            case DataType.GENERAL:
                title = "Discrete Traits";
                break;
            case DataType.CONTINUOUS:
                title = "Continuous Traits";
                break;
            case DataType.MICRO_SAT:
                title = "Microsatellite";
                break;
            default:
                throw new IllegalArgumentException("Unsupported data type");

        }
        modelBorder.setTitle(title + " Substitution Model - " + currentModel.getName());
        repaint();
    }

    private boolean isUsed(int row) {
        PartitionSubstitutionModel model = options.getPartitionSubstitutionModels().get(row);
        for (AbstractPartitionData partition : options.dataPartitions) {
            if (partition.getPartitionSubstitutionModel() == model) {
                return true;
            }
        }
        return false;
    }

    public JComponent getExportableComponent() {
        return this;
    }

    class ModelTableModel extends AbstractTableModel {

        /**
         *
         */
        private static final long serialVersionUID = -6707994233020715574L;
        String[] columnNames = {"Substitution Model"};

        public ModelTableModel() {
        }

        public int getColumnCount() {
            return columnNames.length;
        }

        public int getRowCount() {
            if (options == null) return 0;
            return options.getPartitionSubstitutionModels().size();
        }

        public Object getValueAt(int row, int col) {
            PartitionSubstitutionModel model = options.getPartitionSubstitutionModels().get(row);
            switch (col) {
                case 0:
                    return model.getName();
                default:
                    throw new IllegalArgumentException("unknown column, " + col);
            }
        }

        public boolean isCellEditable(int row, int col) {
            return true;
        }

        public void setValueAt(Object value, int row, int col) {
            String name = ((String) value).trim();
            if (name.length() > 0) {
                PartitionSubstitutionModel model = options.getPartitionSubstitutionModels().get(row);
                model.setName(name); //TODO: update every same model in diff PD?
                updateBorder();
                fireModelsChanged();
            }
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

//    Action addModelAction = new AbstractAction("+") {
//        public void actionPerformed(ActionEvent ae) {
//            createModel();
//        }
//    };
}