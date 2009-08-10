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
import dr.app.beauti.components.SequenceErrorType;
import dr.app.beauti.options.BeautiOptions;
import dr.app.beauti.options.ClockType;
import dr.app.beauti.options.FixRateType;
import dr.app.beauti.options.PartitionClockModel;
import dr.app.beauti.options.TreePrior;
import dr.app.beauti.util.PanelUtils;

import org.virion.jam.components.RealNumberField;
import org.virion.jam.framework.Exportable;
import org.virion.jam.panels.OptionsPanel;
import org.virion.jam.table.HeaderRenderer;
import org.virion.jam.table.TableEditorStopper;

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
    
    JComboBox rateOptionCombo = new JComboBox(FixRateType.values());
//    JCheckBox fixedSubstitutionRateCheck = new JCheckBox("Fix mean substitution rate:");
//    JLabel substitutionRateLabel = new JLabel("Mean substitution rate:");
    RealNumberField substitutionRateField = new RealNumberField(Double.MIN_VALUE, Double.MAX_VALUE);
    
    JComboBox errorModelCombo = new JComboBox(SequenceErrorType.values());
    
    SequenceErrorModelComponentOptions comp;

    BeautiFrame frame = null;
    BeautiOptions options = null;
    boolean settingOptions = false;

    public ClockModelsPanel(BeautiFrame parent) {

		this.frame = parent;

		dataTableModel = new DataTableModel();
		dataTable = new JTable(dataTableModel);

		dataTable.getTableHeader().setReorderingAllowed(false);
		dataTable.getTableHeader().setDefaultRenderer(new HeaderRenderer(SwingConstants.LEFT, new Insets(0, 4, 0, 4)));

		TableColumn col = dataTable.getColumnModel().getColumn(1);
		ComboBoxRenderer comboBoxRenderer = new ComboBoxRenderer();
		comboBoxRenderer.putClientProperty("JComboBox.isTableCellEditor", Boolean.TRUE);
		col.setCellRenderer(comboBoxRenderer);

		TableEditorStopper.ensureEditingStopWhenTableLosesFocus(dataTable);

		JScrollPane scrollPane = new JScrollPane(dataTable,
				JScrollPane.VERTICAL_SCROLLBAR_ALWAYS,
				JScrollPane.HORIZONTAL_SCROLLBAR_NEVER);
		scrollPane.setOpaque(false);

		ItemListener comboListener = new ItemListener() {
			public void itemStateChanged(ItemEvent event) {
				fireModelsChanged();
			}
		};

		PanelUtils.setupComponent(errorModelCombo);
		errorModelCombo
				.setToolTipText("<html>Select how to model sequence error or<br>"
						+ "post-mortem DNA damage.</html>");
		errorModelCombo.addItemListener(comboListener);

		// PanelUtils.setupComponent(clockModelCombo);
		// clockModelCombo.setToolTipText("<html>Select either a strict molecular clock or<br>or a relaxed clock model.</html>");
		// clockModelCombo.addItemListener(comboListener);

		PanelUtils.setupComponent(rateOptionCombo);
		rateOptionCombo.addItemListener(new ItemListener() {
			public void itemStateChanged(ItemEvent ev) {
				if (ev.getStateChange() == ItemEvent.SELECTED) {
					substitutionRateField
							.setEnabled((FixRateType) rateOptionCombo
									.getSelectedItem() != FixRateType.ESTIMATE);
					if ((FixRateType) rateOptionCombo.getSelectedItem() == FixRateType.FIX_MEAN) {
						JOptionPane
								.showMessageDialog(
										rateOptionCombo,
										"WARNING: 'Fix mean rate' function is only working for single locus in this release.");
					}
				}
			}
		});
		rateOptionCombo.setToolTipText("<html>Select this option to fix the substitution rate<br>"
						+ "rather than try to infer it. If this option is<br>"
						+ "turned off then either the sequences should have<br>"
						+ "dates or the tree should have sufficient calibration<br>"
						+ "informations specified as priors.</html>");// TODO Alexei

		PanelUtils.setupComponent(substitutionRateField);
		substitutionRateField.setValue(1.0);
		substitutionRateField.addKeyListener(new java.awt.event.KeyAdapter() {
			public void keyTyped(java.awt.event.KeyEvent ev) {
				frame.setDirty();
			}
		});
		substitutionRateField.setToolTipText("<html>Enter the fixed mean rate or 1st partition rate here.</html>");
		substitutionRateField.setEnabled(true);

		JPanel modelPanelParent = new JPanel(new FlowLayout(FlowLayout.CENTER));
        modelPanelParent.setOpaque(false);
        TitledBorder modelBorder = new TitledBorder("Overall model(s) parameters:");
        modelPanelParent.setBorder(modelBorder);
		
		OptionsPanel panel = new OptionsPanel(12, 20);

		panel.addComponentWithLabel("Sequence Error Model:", errorModelCombo);
		// panel.addComponentWithLabel("Molecular Clock Model:", clockModelCombo);
		panel.addComponentWithLabel("Fixed Rate Option:", rateOptionCombo);
		substitutionRateField.setColumns(10);
		panel.addComponentWithLabel("Fixed mean rate / 1st partition rate:", substitutionRateField);

		modelPanelParent.add(panel);
		
		JSplitPane splitPane = new JSplitPane(JSplitPane.HORIZONTAL_SPLIT, scrollPane, modelPanelParent);
		splitPane.setDividerLocation(400);
		splitPane.setContinuousLayout(true);
		splitPane.setBorder(BorderFactory.createEmptyBorder());
		splitPane.setOpaque(false);

		setOpaque(false);
		setLayout(new BorderLayout(0, 0));
		setBorder(new BorderUIResource.EmptyBorderUIResource(new Insets(12, 12, 12, 12)));
		add(splitPane, BorderLayout.CENTER);

		comp = new SequenceErrorModelComponentOptions();
    }

    private void fireDataChanged() {
        options.updatePartitionClockTreeLinks();
        frame.setDirty();
    }

    private void modelsChanged() {
        TableColumn col = dataTable.getColumnModel().getColumn(1);
        col.setCellEditor(new DefaultCellEditor(new JComboBox(EnumSet.range(ClockType.STRICT_CLOCK, ClockType.UNCORRELATED_LOGNORMAL).toArray())));
    }
    
    private void fireModelsChanged() {
        options.updatePartitionClockTreeLinks();
        options.updateFixedRateClockModel();
        frame.setDirty();
    }

    public void setOptions(BeautiOptions options) {

        this.options = options;

        settingOptions = true;
        
//      clockModelCombo.setSelectedItem(options.clockType);
        comp = (SequenceErrorModelComponentOptions) options.getComponentOptions(SequenceErrorModelComponentOptions.class);
        errorModelCombo.setSelectedItem(comp.errorModelType);
      
        rateOptionCombo.setSelectedItem(options.rateOptionClockModel);
        substitutionRateField.setValue(options.meanSubstitutionRate);  
        
        settingOptions = false;
        
        int selRow = dataTable.getSelectedRow();
        dataTableModel.fireTableDataChanged();
        if (options.getPartitionClockModels().size() > 0) {
            if (selRow < 0) {
                selRow = 0;
            }
            dataTable.getSelectionModel().setSelectionInterval(selRow, selRow);
        }
        
        fireModelsChanged();

        modelsChanged();

        dataTableModel.fireTableDataChanged();
    }

    public void getOptions(BeautiOptions options) {
    	if (settingOptions) return;
    	
    	SequenceErrorModelComponentOptions comp = (SequenceErrorModelComponentOptions) options.getComponentOptions(SequenceErrorModelComponentOptions.class);
        comp.errorModelType = (SequenceErrorType) errorModelCombo.getSelectedItem();

        options.rateOptionClockModel = (FixRateType) rateOptionCombo.getSelectedItem();       
        options.meanSubstitutionRate = substitutionRateField.getValue();
       
        fireModelsChanged();    	
    }

    public JComponent getExportableComponent() {
        return dataTable;
    }

    class DataTableModel extends AbstractTableModel {

        private static final long serialVersionUID = -2852144669936634910L;

        String[] columnNames = {"Clock Model Name", "Molecular Clock Model"};

        public DataTableModel() {
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
                case 1:
                    return model.getClockType();
                default:
                    throw new IllegalArgumentException("unknown column, " + col);
            }
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
            }
            fireDataChanged();
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