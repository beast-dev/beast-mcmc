/*
 * ClockModelPanel.java
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
import dr.app.beauti.options.BeautiOptions;
import dr.app.beauti.options.PartitionClockModel;
import dr.app.beauti.types.ClockType;
import dr.app.beauti.types.FixRateType;
import dr.app.beauti.util.PanelUtils;
import dr.app.gui.components.RealNumberField;
import dr.app.gui.table.RealNumberCellEditor;
import dr.app.gui.table.TableEditorStopper;
import jam.framework.Exportable;
import jam.panels.OptionsPanel;
import jam.table.TableRenderer;

import javax.swing.*;
import javax.swing.border.TitledBorder;
import javax.swing.plaf.BorderUIResource;
import javax.swing.table.AbstractTableModel;
import javax.swing.table.JTableHeader;
import javax.swing.table.TableColumn;
import java.awt.*;
import java.awt.event.ItemEvent;
import java.awt.event.ItemListener;
import java.awt.event.MouseEvent;
import java.util.EnumSet;

/**
 * @author Andrew Rambaut
 * @author Alexei Drummond
 * @author Walter Xie
 * @version $Id: ClockModelPanel.java,v 1.17 2006/09/05 13:29:34 rambaut Exp $
 */
public class OldClockModelsPanel extends BeautiPanel implements Exportable {

    private static final long serialVersionUID = 2945922234432540027L;
    private final String[] columnToolTips = {"Name", "Clock model",
            "Decide whether to estimate this clock model",
            "Provide the rate if it is fixed"};
    private final String[] columnToolTips2 = {"A group of clock models",
            "<html>Fix mean rate of this group of clock models." +
                    "<br>Select this option to fix the mean substitution rate,<br>" +
                    "rather than try to infer it. If this option is turned off, then<br>" +
                    "either the sequences should have dates or the tree should have<br>" +
                    "sufficient calibration informations specified as priors.<br>" +
                    "In addition, it is only available for multi-clock partitions.</html>",
            "Enter the fixed mean rate here."};
    private static final int MINIMUM_TABLE_HEIGHT = 400;

    JTable dataTable = null;
    DataTableModel dataTableModel = null;
    JScrollPane scrollPane;
    JCheckBox fixedMeanRateCheck = new JCheckBox("Fix mean rate of molecular clock model to: ");
    RealNumberField meanRateField = new RealNumberField(Double.MIN_VALUE, Double.MAX_VALUE);

    JTable clockGroupTable = null;
    ClockGroupTableModel clockGroupTableModel = null;

    BeautiFrame frame = null;
    BeautiOptions options = null;
    boolean settingOptions = false;

    public OldClockModelsPanel(BeautiFrame parent) {

        this.frame = parent;

        dataTableModel = new DataTableModel();
        dataTable = new JTable(dataTableModel) {
            //Implement table header tool tips.
            protected JTableHeader createDefaultTableHeader() {
                return new JTableHeader(columnModel) {
                    public String getToolTipText(MouseEvent e) {
                        Point p = e.getPoint();
                        int index = columnModel.getColumnIndexAtX(p.x);
                        int realIndex = columnModel.getColumn(index).getModelIndex();
                        return columnToolTips[realIndex];
                    }
                };
            }
        };

        initTable(dataTable);

        scrollPane = new JScrollPane(dataTable,
                JScrollPane.VERTICAL_SCROLLBAR_AS_NEEDED,
                JScrollPane.HORIZONTAL_SCROLLBAR_AS_NEEDED);
        scrollPane.setOpaque(false);

        // PanelUtils.setupComponent(clockModelCombo);
        // clockModelCombo.setToolTipText("<html>Select either a strict molecular clock or<br>or a relaxed clock model.</html>");
        // clockModelCombo.addItemListener(comboListener);

        PanelUtils.setupComponent(fixedMeanRateCheck);
        fixedMeanRateCheck.setSelected(false); // default to FixRateType.ESTIMATE
        fixedMeanRateCheck.addItemListener(new ItemListener() {
            public void itemStateChanged(ItemEvent ev) {
                meanRateField.setEnabled(fixedMeanRateCheck.isSelected());
                if (fixedMeanRateCheck.isSelected()) {
                    options.clockModelOptions.fixMeanRate();
                } else {
                    options.clockModelOptions.fixRateOfFirstClockPartition();
                }

                dataTableModel.fireTableDataChanged();
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
//		meanRateField.setEnabled(true);

        JPanel modelPanelParent = new JPanel(new BorderLayout(12, 12));
//        modelPanelParent.setLayout(new BoxLayout(modelPanelParent, BoxLayout.Y_AXIS));
        modelPanelParent.setOpaque(false);
        TitledBorder modelBorder = new TitledBorder("Clock Model : ");
        modelPanelParent.setBorder(modelBorder);

        OptionsPanel panel = new OptionsPanel(12, 12);
        panel.addComponents(fixedMeanRateCheck, meanRateField);


        // The bottom panel is now small enough that this is not necessary
//        JScrollPane scrollPane2 = new JScrollPane(panel);
//        scrollPane2.setOpaque(false);
//        scrollPane2.setPreferredSize(new Dimension(400, 150));

        modelPanelParent.add(scrollPane, BorderLayout.CENTER);
        modelPanelParent.add(panel, BorderLayout.SOUTH);

        //=======  Clock Model Group for Fix Mean function ==========
        clockGroupTableModel = new ClockGroupTableModel();
        clockGroupTable = new JTable(clockGroupTableModel) {
            //Implement table header tool tips.
            protected JTableHeader createDefaultTableHeader() {
                return new JTableHeader(columnModel) {
                    public String getToolTipText(MouseEvent e) {
                        Point p = e.getPoint();
                        int index = columnModel.getColumnIndexAtX(p.x);
                        int realIndex = columnModel.getColumn(index).getModelIndex();
                        return columnToolTips2[realIndex];
                    }
                };
            }
        };
        clockGroupTable.setAutoResizeMode(JTable.AUTO_RESIZE_OFF);
        clockGroupTable.getTableHeader().setReorderingAllowed(false);

        TableColumn col = clockGroupTable.getColumnModel().getColumn(0);
        col.setMinWidth(200);
        col = clockGroupTable.getColumnModel().getColumn(1);
        col.setMinWidth(40);
        col = clockGroupTable.getColumnModel().getColumn(2);
        col.setCellEditor(new RealNumberCellEditor(0, Double.POSITIVE_INFINITY));
        col.setMinWidth(80);
        TableEditorStopper.ensureEditingStopWhenTableLosesFocus(clockGroupTable);

        JScrollPane d_scrollPane = new JScrollPane(clockGroupTable,
                JScrollPane.VERTICAL_SCROLLBAR_AS_NEEDED,
                JScrollPane.HORIZONTAL_SCROLLBAR_AS_NEEDED);
        d_scrollPane.setOpaque(false);
        TitledBorder traitClockBorder = new TitledBorder("Clock Model Group: ");
        d_scrollPane.setBorder(traitClockBorder);

        JSplitPane splitPane = new JSplitPane(JSplitPane.VERTICAL_SPLIT, modelPanelParent, d_scrollPane);
        splitPane.setDividerLocation(MINIMUM_TABLE_HEIGHT);

        setOpaque(false);
        setLayout(new BorderLayout(12, 12));
        setBorder(new BorderUIResource.EmptyBorderUIResource(new Insets(12, 12, 12, 12)));
        add(splitPane, BorderLayout.CENTER);
    }

    private void initTable(JTable dataTable) {
        dataTable.setAutoResizeMode(JTable.AUTO_RESIZE_OFF);
        dataTable.getTableHeader().setReorderingAllowed(false);
//        dataTable.getTableHeader().setDefaultRenderer(
//              new HeaderRenderer(SwingConstants.LEFT, new Insets(0, 4, 0, 4)));

        TableColumn col = dataTable.getColumnModel().getColumn(0);
        col.setCellRenderer(new ClockTableCellRenderer(SwingConstants.LEFT, new Insets(0, 4, 0, 4)));
        col.setMinWidth(200);

        col = dataTable.getColumnModel().getColumn(1);
        ComboBoxRenderer comboBoxRenderer = new ComboBoxRenderer();
        comboBoxRenderer.putClientProperty("JComboBox.isTableCellEditor", Boolean.TRUE);
        col.setCellRenderer(comboBoxRenderer);
        col.setMinWidth(260);

        col = dataTable.getColumnModel().getColumn(2);
        col.setMinWidth(40);

        col = dataTable.getColumnModel().getColumn(3);
        col.setCellRenderer(new ClockTableCellRenderer(SwingConstants.LEFT, new Insets(0, 4, 0, 4)));
        col.setCellEditor(new RealNumberCellEditor(0, Double.POSITIVE_INFINITY));
        col.setMinWidth(80);

        TableEditorStopper.ensureEditingStopWhenTableLosesFocus(dataTable);
    }

    private void modelsChanged() {
        TableColumn col = dataTable.getColumnModel().getColumn(1);
        col.setCellEditor(new DefaultCellEditor(new JComboBox(EnumSet.range(ClockType.STRICT_CLOCK, ClockType.RANDOM_LOCAL_CLOCK).toArray())));

        options.clockModelOptions.updateClockModelGroup();
        clockGroupTableModel.fireTableDataChanged();
    }

    private void fireModelsChanged() {
        options.updatePartitionAllLinks();
        frame.setStatusMessage();
        frame.setDirty();
    }

//    private void updateModelPanelBorder() {
//    	if (options.hasData()) {
//    		modelBorder.setTitle(options.clockModelOptions.getRateOptionClockModel().toString());
//    	} else {
//    		modelBorder.setTitle("Overall clock model(s) parameters");
//    	}
//
//        repaint();
//    }

    public void setOptions(BeautiOptions options) {

        this.options = options;

        settingOptions = true;

        fixedMeanRateCheck.setSelected(options.clockModelOptions.getRateOptionClockModel() == FixRateType.FIX_MEAN);
        fixedMeanRateCheck.setEnabled(!(options.clockModelOptions.getRateOptionClockModel() == FixRateType.TIP_CALIBRATED
                || options.clockModelOptions.getRateOptionClockModel() == FixRateType.NODE_CALIBRATED
                || options.clockModelOptions.getRateOptionClockModel() == FixRateType.RATE_CALIBRATED));
        meanRateField.setValue(options.clockModelOptions.getMeanRelativeRate());

        settingOptions = false;

        int selRow = dataTable.getSelectedRow();
        dataTableModel.fireTableDataChanged();
        if (options.getPartitionClockModels().size() > 0) {
            if (selRow < 0) {
                selRow = 0;
            }
            dataTable.getSelectionModel().setSelectionInterval(selRow, selRow);
        }

//        fireModelsChanged();

        modelsChanged();

        dataTableModel.fireTableDataChanged();
    }

    public void getOptions(BeautiOptions options) {
        if (settingOptions) return;

//        if (fixedMeanRateCheck.isSelected()) {
//        	options.clockModelOptions.fixMeanRate();
//        } else {
//        	options.clockModelOptions.fixRateOfFirstClockPartition();
//        }
        options.clockModelOptions.setMeanRelativeRate(meanRateField.getValue());

//        fireModelsChanged();
    }

    public JComponent getExportableComponent() {
        return dataTable;
    }

    class DataTableModel extends AbstractTableModel {

        private static final long serialVersionUID = -2852144669936634910L;

        //        String[] columnNames = {"Clock Model Name", "Molecular Clock Model"};
        String[] columnNames = {"Name", "Model", "Estimate", "Rate"};

        public DataTableModel() {
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
                    model.setRate((Double) aValue);
                    options.selectParameters();
                    break;
                default:
                    throw new IllegalArgumentException("unknown column, " + col);
            }
            fireModelsChanged();
        }

        public boolean isCellEditable(int row, int col) {
            boolean editable;

            switch (col) {
                case 2:// Check box
                    editable = !fixedMeanRateCheck.isSelected();
                    break;
                case 3:
                    editable = !fixedMeanRateCheck.isSelected() && !((Boolean) getValueAt(row, 2));
                    break;
                default:
                    editable = true;
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

    class ClockGroupTableModel extends AbstractTableModel {

        String[] columnNames = {"Group", "Fix Mean", "Rate"};

        public ClockGroupTableModel() {
        }

        public int getColumnCount() {
            return columnNames.length;
        }

        public int getRowCount() {
            if (options == null) return 0;
            return options.clockModelOptions.clockModelGroupNameList.size();
        }

        public Object getValueAt(int row, int col) {
            switch (col) {
                case 0:
                    return options.clockModelOptions.clockModelGroupNameList.get(row);
                case 1:
                    return options.clockModelOptions.fixMeanList.get(row);
                case 2:
                    return options.clockModelOptions.fixMeanRateList.get(row);
            }
            return null;
        }

        public void setValueAt(Object aValue, int row, int col) {

            switch (col) {
                case 0:
                    String name = ((String) aValue).trim();
                    if (name.length() > 0) {
                        options.clockModelOptions.clockModelGroupNameList.set(row, name);
                    }
                    break;
                case 1:
                    options.clockModelOptions.fixMeanList.set(row, (Boolean) aValue);
                    break;
                case 2:
                    options.clockModelOptions.fixMeanRateList.set(row, (Double) aValue);
                    break;

                default:
                    throw new IllegalArgumentException("unknown column, " + col);
            }

        }

        public boolean isCellEditable(int row, int col) {
            boolean editable;

            switch (col) {
                case 1:// Check box
                    editable = true;//!fixedMeanRateCheck.isSelected();
                    break;
                case 2:
                    editable = /*!fixedMeanRateCheck.isSelected() &&*/ ((Boolean) getValueAt(row, 1));
                    break;
                default:
                    editable = true;
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

}