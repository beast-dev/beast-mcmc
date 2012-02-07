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

package dr.app.beauti.clockModelsPanel;

import dr.app.beauti.BeautiFrame;
import dr.app.beauti.BeautiPanel;
import dr.app.beauti.ComboBoxRenderer;
import dr.app.beauti.options.*;
import dr.app.beauti.types.ClockType;
import dr.app.beauti.util.PanelUtils;
import dr.app.gui.components.RealNumberField;
import dr.app.gui.table.RealNumberCellEditor;
import dr.app.gui.table.TableEditorStopper;
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
import java.awt.*;
import java.awt.event.ItemEvent;
import java.awt.event.ItemListener;
import java.util.EnumSet;
import java.util.HashMap;
import java.util.Map;

/**
 * @author Andrew Rambaut
 * @author Alexei Drummond
 * @author Walter Xie
 * @version $Id: ModelPanel.java,v 1.17 2006/09/05 13:29:34 rambaut Exp $
 * @deprecated
 */
public class ClockModelsPanel extends BeautiPanel implements Exportable {

    public final static boolean DEBUG = false;

    private static final long serialVersionUID = 2778103564318492601L;

    private static final int MINIMUM_TABLE_WIDTH = 140;

    private final String[] columnToolTips = {null, "Molecular clock model",
            "Decide whether to estimate molecular clock model",
            "Provide the rate if it is fixed"};

    JTable clockTable = null;
    ClockTableModel clockTableModel = null;
    BeautiOptions options = null;

    JPanel modelPanelParent;
    PartitionClockModel currentModel = null;
    Map<PartitionClockModel, PartitionClockModelPanel> modelPanels = new HashMap<PartitionClockModel, PartitionClockModelPanel>();
    TitledBorder modelBorder;

    JCheckBox fixedMeanRateCheck = new JCheckBox("Fix mean rate of molecular clock model to: ");
    RealNumberField meanRateField = new RealNumberField(Double.MIN_VALUE, Double.MAX_VALUE);

    BeautiFrame frame = null;
//    CreateModelDialog createModelDialog = null;
    boolean settingOptions = false;


    public ClockModelsPanel(BeautiFrame parent) {

        super();

        this.frame = parent;

        clockTableModel = new ClockTableModel();
        clockTable = new JTable(clockTableModel); // {
            //Implement table header tool tips.
//            protected JTableHeader createDefaultTableHeader() {
//                return new JTableHeader(columnModel) {
//                    public String getToolTipText(MouseEvent e) {
//                        Point p = e.getPoint();
//                        int index = columnModel.getColumnIndexAtX(p.x);
//                        int realIndex = columnModel.getColumn(index).getModelIndex();
//                        return columnToolTips[realIndex];
//                    }
//                };
//            }
//        };

        initTable(clockTable);

        clockTable.getSelectionModel().addListSelectionListener(new ListSelectionListener() {
            public void valueChanged(ListSelectionEvent evt) {
                selectionChanged();
            }
        });

        JScrollPane scrollPane = new JScrollPane(clockTable,
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
        splitPane.setDividerLocation(0.5);
        splitPane.setContinuousLayout(true);
        splitPane.setBorder(BorderFactory.createEmptyBorder());
        splitPane.setOpaque(false);

        PanelUtils.setupComponent(fixedMeanRateCheck);
        fixedMeanRateCheck.setSelected(false); // default to FixRateType.ESTIMATE
        fixedMeanRateCheck.addItemListener(new ItemListener() {
            public void itemStateChanged(ItemEvent ev) {
                meanRateField.setEnabled(fixedMeanRateCheck.isSelected());
//                if (fixedMeanRateCheck.isSelected()) {
//                    options.clockModelOptions.fixMeanRate();
//                } else {
//                    options.clockModelOptions.fixRateOfFirstClockPartition();
//                }

                clockTableModel.fireTableDataChanged();
                fireModelsChanged();
            }
        });
        fixedMeanRateCheck.setToolTipText("<html>Select this option to fix the mean substitution rate,<br>"
                + "rather than try to infer it. If this option is turned off, then<br>"
                + "either the sequences should have dates or the tree should have<br>"
                + "sufficient calibration informations specified as priors.<br>"
                + "In addition, it is only available for multi-clock partitions." + "</html>");// TODO Alexei

        PanelUtils.setupComponent(meanRateField);
        meanRateField.setEnabled(fixedMeanRateCheck.isSelected());
        meanRateField.setValue(1.0);
        meanRateField.addKeyListener(new java.awt.event.KeyAdapter() {
            public void keyTyped(java.awt.event.KeyEvent ev) {
                frame.setDirty();
            }
        });
        meanRateField.setToolTipText("<html>Enter the fixed mean rate here.</html>");
        meanRateField.setColumns(10);

        OptionsPanel panel2 = new OptionsPanel(12, 12);
        panel2.addComponents(fixedMeanRateCheck, meanRateField);

        setOpaque(false);
        setBorder(new BorderUIResource.EmptyBorderUIResource(new Insets(12, 12, 12, 12)));
        setLayout(new BorderLayout(0, 0));
        add(splitPane, BorderLayout.CENTER);
        add(panel2, BorderLayout.SOUTH);
    }

    private void initTable(JTable clockTable){
//        clockTable.setAutoResizeMode(JTable.AUTO_RESIZE_OFF);
        clockTable.getTableHeader().setReorderingAllowed(false);
//        clockTable.getTableHeader().setDefaultRenderer(new HeaderRenderer(SwingConstants.LEFT, new Insets(0, 4, 0, 4)));

        TableColumn col = clockTable.getColumnModel().getColumn(0);
        col.setCellRenderer(new ClockTableCellRenderer(SwingConstants.LEFT, new Insets(0, 4, 0, 4)));
//        col.setMinWidth(80);

        col = clockTable.getColumnModel().getColumn(1);
        ComboBoxRenderer comboBoxRenderer = new ComboBoxRenderer();
        comboBoxRenderer.putClientProperty("JComboBox.isTableCellEditor", Boolean.TRUE);
//        col.setMinWidth(260);

        col = clockTable.getColumnModel().getColumn(2);
        col.setMinWidth(60);
        col.setMaxWidth(60);

        col = clockTable.getColumnModel().getColumn(3);
        col.setCellRenderer(new ClockTableCellRenderer(SwingConstants.LEFT, new Insets(0, 4, 0, 4)));
        col.setCellEditor(new RealNumberCellEditor(0, Double.POSITIVE_INFINITY));
//        col.setMinWidth(80);

        TableEditorStopper.ensureEditingStopWhenTableLosesFocus(clockTable);
    }

    private void modelsChanged() {
        TableColumn col = clockTable.getColumnModel().getColumn(1);
        col.setCellEditor(new DefaultCellEditor(new JComboBox(EnumSet.range(ClockType.STRICT_CLOCK, ClockType.RANDOM_LOCAL_CLOCK).toArray())));
    }

    private void resetPanel() {
        if (!options.hasData()) {
            currentModel = null;
            modelPanels.clear();
            modelPanelParent.removeAll();
            modelBorder.setTitle("Clock Model");

            return;
        }
    }

    public void setOptions(BeautiOptions options) {

        this.options = options;

        resetPanel();

        settingOptions = true;

//        fixedMeanRateCheck.setSelected(options.clockModelOptions.getRateOptionClockModel() == FixRateType.FIX_MEAN);
//        fixedMeanRateCheck.setEnabled(!(options.clockModelOptions.getRateOptionClockModel() == FixRateType.TIP_CALIBRATED
//                || options.clockModelOptions.getRateOptionClockModel() == FixRateType.NODE_CALIBRATED
//                || options.clockModelOptions.getRateOptionClockModel() == FixRateType.RATE_CALIBRATED));
//        meanRateField.setValue(options.clockModelOptions.getMeanRelativeRate());

        int selRow = clockTable.getSelectedRow();
        clockTableModel.fireTableDataChanged();
        if (options.getPartitionSubstitutionModels().size() > 0) {
            if (selRow < 0) {
                selRow = 0;
            }
            clockTable.getSelectionModel().setSelectionInterval(selRow, selRow);
        }

        if (currentModel == null && options.getPartitionClockModels().size() > 0) {
            clockTable.getSelectionModel().setSelectionInterval(0, 0);
        }

        modelsChanged();

        settingOptions = false;
    }

    public void getOptions(BeautiOptions options) {
        if (settingOptions) return;

//        options.clockModelOptions.setMeanRelativeRate(meanRateField.getValue());
    }


    private void fireModelsChanged() {
        options.updatePartitionAllLinks();
        frame.setDirty();
    }

    private void selectionChanged() {
        int selRow = clockTable.getSelectedRow();

        if (selRow >= options.getPartitionClockModels().size()) {
            selRow = 0;
            clockTable.getSelectionModel().setSelectionInterval(selRow, selRow);
        }

        if (selRow >= 0) {
            setCurrentModel(options.getPartitionClockModels().get(selRow));
//            frame.modelSelectionChanged(!isUsed(selRow));
        }
    }

    /**
     * Sets the current model that this model panel is displaying
     *
     * @param model the new model to display
     */
    private void setCurrentModel(PartitionClockModel model) {

        if (model != null) {
            if (currentModel != null) modelPanelParent.removeAll();

            PartitionClockModelPanel panel = modelPanels.get(model);
            if (panel == null) {
                panel = new PartitionClockModelPanel(model);
                modelPanels.put(model, panel);
            }

            currentModel = model;
            modelPanelParent.add(panel);

            updateBorder();
        }
    }

    private void updateBorder() {
        modelBorder.setTitle("Clock Model - " + currentModel.getName());
        repaint();
    }

    private boolean isUsed(int row) {
        PartitionClockModel model = options.getPartitionClockModels().get(row);
        for (AbstractPartitionData partition : options.dataPartitions) {
            if (partition.getPartitionClockModel() == model) {
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
        String[] columnNames = {"Clock Model"};

        public ModelTableModel() {
        }

        public int getColumnCount() {
            return columnNames.length;
        }

        public int getRowCount() {
            if (options == null) return 0;
            return options.getPartitionClockModels().size();
        }

        public Object getValueAt(int row, int col) {
            PartitionClockModel model = options.getPartitionClockModels().get(row);
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
                PartitionClockModel model = options.getPartitionClockModels().get(row);
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


    class ClockTableModel extends AbstractTableModel {
        private static final long serialVersionUID = -2852144669936634910L;

//        String[] columnNames = {"Clock Model Name", "Molecular Clock Model"};
        String[] columnNames = {"Name", "Model", "Estimate", "Rate"};

        public ClockTableModel() {
        }

        public int getColumnCount() {
//        	if (estimateRelatieRateCheck.isSelected()) {
//        		return columnNames2.length;
//        	} else {
            return columnNames.length;
//        	}
        }

        public int getRowCount() {
            if (options == null) return 0;
            if (options.getPartitionClockModels().size() < 2) {
                fixedMeanRateCheck.setEnabled(false);
            } else {
                fixedMeanRateCheck.setEnabled(true);
            }
            return options.getPartitionClockModels().size();
        }

        public Object getValueAt(int row, int col) {
            PartitionClockModel model = options.getPartitionClockModels().get(row);
            switch (col) {
                case 0:
                    return model.getName();
                case 1:
                    return model.getClockType();
                case 2:
                    return model.isEstimatedRate();
                case 3:
                    return model.getRate();
            }
            return null;
        }

        public void setValueAt(Object aValue, int row, int col) {
            PartitionClockModel model = options.getPartitionClockModels().get(row);
            switch (col) {
                case 0:
                    String name = ((String) aValue).trim();
                    if (name.length() > 0) {
                        model.setName(name);
                    }
                    break;
                case 1:
                    model.setClockType((ClockType) aValue);
                    break;
                case 2:
                    model.setEstimatedRate((Boolean) aValue);
//                    if (options.clockModelOptions.getRateOptionClockModel() == FixRateType.RElATIVE_TO) {
//                        if (!options.clockModelOptions.validateRelativeTo()) {
//                            JOptionPane.showMessageDialog(frame, "It must have at least one clock rate to be fixed !",
//                                    "Validation Of Relative To ?th Rate", JOptionPane.WARNING_MESSAGE);
//                            model.setEstimatedRate(false);
//                        }
//                    }
                    break;
                case 3:
                    model.setRate((Double) aValue, true);
                    options.selectParameters();
                    break;
                default:
                    throw new IllegalArgumentException("unknown column, " + col);
            }
            fireModelsChanged();
        }

        public boolean isCellEditable(int row, int col) {
            switch (col) {
                case 2:// Check box
                    return !fixedMeanRateCheck.isSelected();
                case 3:
                    return !fixedMeanRateCheck.isSelected() && !((Boolean) getValueAt(row, 2));
                default:
                    return true;
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


    class ClockTableCellRenderer extends TableRenderer {

        public ClockTableCellRenderer(int alignment, Insets insets) {
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

            if (fixedMeanRateCheck.isSelected() && aColumn > 1) {
                renderer.setForeground(Color.gray);
            } else if (!fixedMeanRateCheck.isSelected() && aColumn == 3 && (Boolean) aTable.getValueAt(aRow, 2)) {
                renderer.setForeground(Color.gray);
            } else {
                renderer.setForeground(Color.black);
            }

            return this;
        }

    }



//    Action addModelAction = new AbstractAction("+") {
//        public void actionPerformed(ActionEvent ae) {
//            createModel();
//        }
//    };
}