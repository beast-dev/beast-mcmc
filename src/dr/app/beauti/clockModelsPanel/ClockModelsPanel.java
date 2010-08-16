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
import dr.app.beauti.components.SequenceErrorModelComponentOptions;
import dr.app.beauti.enumTypes.ClockType;
import dr.app.beauti.enumTypes.FixRateType;
import dr.app.beauti.enumTypes.SequenceErrorType;
import dr.app.beauti.options.BeautiOptions;
import dr.app.beauti.options.PartitionClockModel;
import dr.app.beauti.util.PanelUtils;
import org.virion.jam.components.RealNumberField;
import org.virion.jam.framework.Exportable;
import org.virion.jam.panels.OptionsPanel;
import org.virion.jam.table.HeaderRenderer;
import org.virion.jam.table.RealNumberCellEditor;
import org.virion.jam.table.TableEditorStopper;
import org.virion.jam.table.TableRenderer;

import javax.swing.*;
import javax.swing.border.TitledBorder;
import javax.swing.plaf.BorderUIResource;
import javax.swing.table.AbstractTableModel;
import javax.swing.table.TableColumn;
import java.awt.*;
import java.awt.event.ItemEvent;
import java.awt.event.ItemListener;
import java.util.EnumSet;

/**
 * @author Andrew Rambaut
 * @author Alexei Drummond
 * @author Walter Xie
 * @version $Id: ClockModelPanel.java,v 1.17 2006/09/05 13:29:34 rambaut Exp $
 */
public class ClockModelsPanel extends BeautiPanel implements Exportable {

	private static final long serialVersionUID = 2945922234432540027L;

	JTable dataTable = null;
    DataTableModel dataTableModel = null;
    JScrollPane scrollPane;
    JCheckBox fixedMeanRateCheck = new JCheckBox("Fix mean rate of molecular clock model to: ");
    RealNumberField meanRateField = new RealNumberField(Double.MIN_VALUE, Double.MAX_VALUE);

    BeautiFrame frame = null;
    BeautiOptions options = null;
    boolean settingOptions = false;

    JTable discreteTraitTable = null;
    JScrollPane d_scrollPane;
    boolean activateDiscreteTraitsTable = false;

    public ClockModelsPanel(BeautiFrame parent) {

		this.frame = parent;

		dataTableModel = new DataTableModel();
		dataTable = new JTable(dataTableModel);

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
                // todo Rather than validating, wouldn't it be nicer to simply disable the checkbox
                // todo if molecular clocks are linked or there is just one partition?
			    if (!options.clockModelOptions.validateFixMeanRate(fixedMeanRateCheck)) {
			        JOptionPane.showMessageDialog(frame, "It is only necessary to fix mean substitution rate if multiple molecular clock models are being employed.",
		                    "Validation Of Fix Mean Rate",
		                    JOptionPane.WARNING_MESSAGE);
		            fixedMeanRateCheck.setSelected(false);
		            return;
			    }

				meanRateField.setEnabled(fixedMeanRateCheck.isSelected());
				if (fixedMeanRateCheck.isSelected()) {
		        	options.clockModelOptions.fixMeanRate();
		        } else {
		        	options.clockModelOptions.fixRateOfFirstClockPartition();
		        }

				frame.setDirty();
				frame.repaint();
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

		JPanel modelPanelParent = new JPanel(new BorderLayout(12,12));
//        modelPanelParent.setLayout(new BoxLayout(modelPanelParent, BoxLayout.Y_AXIS));
        modelPanelParent.setOpaque(false);
        TitledBorder modelBorder = new TitledBorder("Molecular Clock Model : ");
        modelPanelParent.setBorder(modelBorder);

		OptionsPanel panel = new OptionsPanel(12, 12);
		panel.addComponents(fixedMeanRateCheck, meanRateField);


        // The bottom panel is now small enough that this is not necessary
//        JScrollPane scrollPane2 = new JScrollPane(panel);
//        scrollPane2.setOpaque(false);
//        scrollPane2.setPreferredSize(new Dimension(400, 150));

		modelPanelParent.add(scrollPane, BorderLayout.CENTER);
        modelPanelParent.add(panel, BorderLayout.SOUTH);

        setOpaque(false);
		setLayout(new BorderLayout(12, 12));
		setBorder(new BorderUIResource.EmptyBorderUIResource(new Insets(12, 12, 12, 12)));
		add(modelPanelParent, BorderLayout.CENTER);

        //=======================  Discrete Trait Substitution Model =========================
        discreteTraitTable = new JTable(new DiscreteTraitModelTableModel());

        initTable(discreteTraitTable);
        d_scrollPane = new JScrollPane(discreteTraitTable,
				JScrollPane.VERTICAL_SCROLLBAR_AS_NEEDED,
				JScrollPane.HORIZONTAL_SCROLLBAR_AS_NEEDED);
		d_scrollPane.setOpaque(false);
        d_scrollPane.setPreferredSize(new Dimension(scrollPane.getWidth(), 150));
        TitledBorder traitClockBorder = new TitledBorder("Trait Clock Model : ");
        d_scrollPane.setBorder(traitClockBorder);

    }

    private void initTable(JTable dataTable){
        dataTable.setAutoResizeMode(javax.swing.JTable.AUTO_RESIZE_OFF);
		dataTable.getTableHeader().setReorderingAllowed(false);
		dataTable.getTableHeader().setDefaultRenderer(new HeaderRenderer(SwingConstants.LEFT, new Insets(0, 4, 0, 4)));

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

        if (!options.hasData() && activateDiscreteTraitsTable) {
    		this.remove(d_scrollPane);
            activateDiscreteTraitsTable = false;
        }

        settingOptions = true;

        fixedMeanRateCheck.setSelected(options.clockModelOptions.getRateOptionClockModel() == FixRateType.FIX_MEAN);
        fixedMeanRateCheck.setEnabled(!(options.clockModelOptions.getRateOptionClockModel() == FixRateType.TIP_CALIBRATED
        		|| options.clockModelOptions.getRateOptionClockModel() == FixRateType.NODE_CALIBRATED
        		|| options.clockModelOptions.getRateOptionClockModel() == FixRateType.RATE_CALIBRATED));
        meanRateField.setValue(options.clockModelOptions.getMeanRelativeRate());

        settingOptions = false;

        int selRow = dataTable.getSelectedRow();
        dataTableModel.fireTableDataChanged();
        if (options.getPartitionNonTraitsClockModels().size() > 0) {
            if (selRow < 0) {
                selRow = 0;
            }
            dataTable.getSelectionModel().setSelectionInterval(selRow, selRow);
        }

//        fireModelsChanged();

        modelsChanged();

        dataTableModel.fireTableDataChanged();

        if (options.hasDiscreteIntegerTraitsExcludeSpecies()) {
            if (!activateDiscreteTraitsTable) {
                this.add(d_scrollPane, BorderLayout.SOUTH);
                activateDiscreteTraitsTable = true;
                scrollPane.setPreferredSize(new Dimension(scrollPane.getWidth(), 300));
            }
        } else {
            this.remove(d_scrollPane);
            activateDiscreteTraitsTable = false;
        }
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
            if (options.getPartitionNonTraitsClockModels().size() < 2) {
            	fixedMeanRateCheck.setEnabled(false);
            } else {
            	fixedMeanRateCheck.setEnabled(true);
            }
            return options.getPartitionNonTraitsClockModels().size();
        }

        public Object getValueAt(int row, int col) {
            PartitionClockModel model = options.getPartitionNonTraitsClockModels().get(row);
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
            PartitionClockModel model = options.getPartitionNonTraitsClockModels().get(row);
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

    class DiscreteTraitModelTableModel extends AbstractTableModel {

        String[] columnNames = {"Name", "Model", "Estimate", "Rate"};

        public DiscreteTraitModelTableModel() {
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
//            System.out.println(options.getPartitionTraitsClockModels().size());
            return options.getPartitionTraitsClockModels().size();
        }

        public Object getValueAt(int row, int col) {
            PartitionClockModel model = options.getPartitionTraitsClockModels().get(row);
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
            PartitionClockModel model = options.getPartitionTraitsClockModels().get(row);
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

        }

        public boolean isCellEditable(int row, int col) {

            return true;
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