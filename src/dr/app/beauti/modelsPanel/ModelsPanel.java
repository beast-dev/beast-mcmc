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

package dr.app.beauti.modelsPanel;

import dr.app.beauti.BeautiFrame;
import dr.app.beauti.BeautiPanel;
import dr.app.beauti.components.SequenceErrorModelComponentOptions;
import dr.app.beauti.components.SequenceErrorType;
import dr.app.beauti.options.BeautiOptions;
import dr.app.beauti.options.FixRateType;
import dr.app.beauti.options.PartitionData;
import dr.app.beauti.options.PartitionSubstitutionModel;
import dr.app.beauti.options.TreePrior;
import dr.app.beauti.util.PanelUtils;
import dr.evolution.datatype.DataType;
import org.virion.jam.components.RealNumberField;
import org.virion.jam.framework.Exportable;
import org.virion.jam.panels.ActionPanel;
import org.virion.jam.panels.OptionsPanel;
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
import java.awt.event.ActionEvent;
import java.awt.event.ItemEvent;
import java.awt.event.ItemListener;
import java.util.EnumSet;
import java.util.Enumeration;
import java.util.HashMap;
import java.util.Map;
import java.util.logging.Logger;

/**
 * @author Andrew Rambaut
 * @author Alexei Drummond
 * @version $Id: ModelPanel.java,v 1.17 2006/09/05 13:29:34 rambaut Exp $
 */
public class ModelsPanel extends BeautiPanel implements Exportable {

    public final static boolean DEBUG = false;

    private static final long serialVersionUID = 2778103564318492601L;

    JScrollPane scrollPane = new JScrollPane();
    JTable modelTable = null;
    ModelTableModel modelTableModel = null;
    BeautiOptions options = null;

    JPanel modelPanelParent;
    PartitionSubstitutionModel currentModel = null;

    Map<PartitionSubstitutionModel, PartitionModelPanel> modelPanels = new HashMap<PartitionSubstitutionModel, PartitionModelPanel>();
    TitledBorder modelBorder;

    // Overall model parameters ////////////////////////////////////////////////////////////////////////

    String overallParas = "Overall substitution model(s) parameters:";
    JComboBox rateOptionCombo = new JComboBox(FixRateType.values());
//    JCheckBox fixedSubstitutionRateCheck = new JCheckBox("Fix mean substitution rate:");
//    JLabel substitutionRateLabel = new JLabel("Mean substitution rate:");
    RealNumberField substitutionRateField = new RealNumberField(Double.MIN_VALUE, Double.MAX_VALUE);
    
    
    //    JComboBox clockModelCombo = new JComboBox(ClockType.values());
    ClockModelPanel clockModelPanel = null;

    JComboBox errorModelCombo = new JComboBox(SequenceErrorType.values());

    ////////////////////////////////////////////////////////////////////////////////////////////////////

    BeautiFrame frame = null;
    CreateModelDialog createModelDialog = null;
    boolean settingOptions = false;
    boolean hasAlignment = false;

    SequenceErrorModelComponentOptions comp;

    public ModelsPanel(BeautiFrame parent, Action removeModelAction) {

        super();

        this.frame = parent;
               
        modelTableModel = new ModelTableModel();
        modelTable = new JTable(modelTableModel);

        modelTable.getTableHeader().setReorderingAllowed(false);
        modelTable.getTableHeader().setResizingAllowed(false);
        modelTable.getTableHeader().setDefaultRenderer(
                new HeaderRenderer(SwingConstants.LEFT, new Insets(0, 4, 0, 4)));

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

        scrollPane = new JScrollPane(modelTable,
                JScrollPane.VERTICAL_SCROLLBAR_ALWAYS,
                JScrollPane.HORIZONTAL_SCROLLBAR_NEVER);
        scrollPane.setOpaque(false);

//        ActionPanel actionPanel1 = new ActionPanel(false);
//        actionPanel1.setAddAction(addModelAction);
//        actionPanel1.setRemoveAction(removeModelAction);

        JPanel controlPanel1 = new JPanel(new FlowLayout(FlowLayout.LEFT));
        controlPanel1.setOpaque(false);
//        controlPanel1.add(actionPanel1);

        ItemListener comboListener = new ItemListener() {
            public void itemStateChanged(ItemEvent event) {
                fireModelsChanged();
            }
        };


        PanelUtils.setupComponent(errorModelCombo);
        errorModelCombo.setToolTipText("<html>Select how to model sequence error or<br>" +
                "post-mortem DNA damage.</html>");
        errorModelCombo.addItemListener(comboListener);

//        PanelUtils.setupComponent(clockModelCombo);
//        clockModelCombo.setToolTipText("<html>Select either a strict molecular clock or<br>or a relaxed clock model.</html>");
//        clockModelCombo.addItemListener(comboListener);
        
        PanelUtils.setupComponent(rateOptionCombo);        
        rateOptionCombo.addItemListener(
                new ItemListener() {
                    public void itemStateChanged(ItemEvent ev) {
                    	if (ev.getStateChange() == ItemEvent.SELECTED) {
	                    	substitutionRateField.setEnabled((FixRateType) rateOptionCombo.getSelectedItem() != FixRateType.ESTIMATE);  
	                    	if ((FixRateType) rateOptionCombo.getSelectedItem() == FixRateType.FIX_MEAN) {
	                    		JOptionPane.showMessageDialog(rateOptionCombo, "WARNING: 'Fix mean rate' function is only working for single locus in this release.");
	                    	}
                    	}
                    }
                }
        );
        rateOptionCombo.setToolTipText(
                "<html>Select this option to fix the substitution rate<br>" +
                        "rather than try to infer it. If this option is<br>" +
                        "turned off then either the sequences should have<br>" +
                        "dates or the tree should have sufficient calibration<br>" +
                        "informations specified as priors.</html>");//TODO Alexei
        
        
        PanelUtils.setupComponent(substitutionRateField);
        substitutionRateField.setValue(1.0);
        substitutionRateField.addKeyListener(new java.awt.event.KeyAdapter() {
            public void keyTyped(java.awt.event.KeyEvent ev) {
                frame.setDirty();
            }
        });
        substitutionRateField.setToolTipText("<html>Enter the fixed mean rate or 1st partition rate here.</html>");
        substitutionRateField.setEnabled(true);
        
        setCurrentModel(null);

        OptionsPanel panel = new OptionsPanel(10, 20);
        panel.addSeparator();

        panel.addLabel(overallParas);
        panel.addComponentWithLabel("Sequence Error Model:", errorModelCombo);
//        panel.addComponentWithLabel("Molecular Clock Model:", clockModelCombo);

        panel.addComponentWithLabel("Fixed Rate Option:", rateOptionCombo);
        
        substitutionRateField.setColumns(10);
        panel.addComponentWithLabel("Fixed mean rate / 1st partition rate:", substitutionRateField);
        
        panel.addSeparator();

        clockModelPanel = new ClockModelPanel(frame);

        JSplitPane splitPane1 = new JSplitPane(JSplitPane.HORIZONTAL_SPLIT, clockModelPanel, panel);
        splitPane1.setDividerLocation(500);
        splitPane1.setContinuousLayout(true);
        splitPane1.setBorder(BorderFactory.createEmptyBorder());
        splitPane1.setOpaque(false);

        JPanel panel1 = new JPanel(new BorderLayout(0, 0));
        panel1.setOpaque(false);
        panel1.add(scrollPane, BorderLayout.CENTER);
        panel1.add(controlPanel1, BorderLayout.SOUTH);

        modelPanelParent = new JPanel(new FlowLayout(FlowLayout.CENTER));
        modelPanelParent.setOpaque(false);
        modelBorder = new TitledBorder("Substitution Model");
        modelPanelParent.setBorder(modelBorder);

        JSplitPane splitPane = new JSplitPane(JSplitPane.HORIZONTAL_SPLIT, panel1, modelPanelParent);
        splitPane.setDividerLocation(180);
        splitPane.setContinuousLayout(true);
        splitPane.setBorder(BorderFactory.createEmptyBorder());
        splitPane.setOpaque(false);

        setOpaque(false);
        setBorder(new BorderUIResource.EmptyBorderUIResource(new Insets(12, 12, 12, 12)));
        setLayout(new BorderLayout(0, 0));
        add(splitPane, BorderLayout.NORTH);
        add(splitPane1, BorderLayout.CENTER);

        comp = new SequenceErrorModelComponentOptions();
    }    
       
    public void setOptions(BeautiOptions options) {

        if (DEBUG) {
            Logger.getLogger("dr.app.beauti").info("ModelsPanel.setOptions");
        }

        this.options = options;

        settingOptions = true;

        //setModelOptions(currentModel);

//        clockModelCombo.setSelectedItem(options.clockType);
        comp = (SequenceErrorModelComponentOptions) options.getComponentOptions(SequenceErrorModelComponentOptions.class);
        errorModelCombo.setSelectedItem(comp.errorModelType);

        if (clockModelPanel != null) {
            clockModelPanel.setOptions(options);
        }

        settingOptions = false;

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
        
        rateOptionCombo.setSelectedItem(options.rateOptionClockModel);
        substitutionRateField.setValue(options.meanSubstitutionRate);   
        
        fireModelsChanged();
        
        validate();
        repaint();
    }

    public void getOptions(BeautiOptions options) {

        // This prevents options be overwritten due to listeners calling
        // this function (indirectly through modelChanged()) whilst in the
        // middle of the setOptions() method.
        if (settingOptions) return;

//        options.clockType = (ClockType) clockModelCombo.getSelectedItem();

        SequenceErrorModelComponentOptions comp = (SequenceErrorModelComponentOptions) options.getComponentOptions(SequenceErrorModelComponentOptions.class);
        comp.errorModelType = (SequenceErrorType) errorModelCombo.getSelectedItem();

        options.rateOptionClockModel = (FixRateType) rateOptionCombo.getSelectedItem();       
        options.meanSubstitutionRate = substitutionRateField.getValue();

        if (clockModelPanel != null) {
            clockModelPanel.getOptions(options);
        }
        
        fireModelsChanged();
    }

    private void fireModelsChanged() {
        options.updatePartitionClockTreeLinks();
        options.updateFixedRateClockModel();
        frame.setDirty();
    }

    private void createModel() {
        if (createModelDialog == null) {
            createModelDialog = new CreateModelDialog(frame);
        }

        int result = createModelDialog.showDialog();
        if (result != JOptionPane.CANCEL_OPTION) {
            PartitionSubstitutionModel model = new PartitionSubstitutionModel(options, createModelDialog.getName(), createModelDialog.getDataType());
//            options.addPartitionSubstitutionModel(model);
            modelTableModel.fireTableDataChanged();
            int row = options.getPartitionSubstitutionModels().size() - 1;
            modelTable.getSelectionModel().setSelectionInterval(row, row);
        }

        fireModelsChanged();
    }

//    public void removeSelection() {
//        int selRow = modelTable.getSelectedRow();
//        if (!isUsed(selRow)) {
//            PartitionSubstitutionModel model = options.getPartitionSubstitutionModels().get(selRow);
//            options.getPartitionSubstitutionModels().remove(model);
//        }
//
//        modelTableModel.fireTableDataChanged();
//        int n = options.getPartitionSubstitutionModels().size();
//        if (selRow >= n) {
//            selRow--;
//        }
//        modelTable.getSelectionModel().setSelectionInterval(selRow, selRow);
//        if (n == 0) {
//            setCurrentModel(null);
//        }
//
//        fireModelsChanged();
//    }

    private void selectionChanged() {

        int selRow = modelTable.getSelectedRow();
        if (selRow >= 0) {
            setCurrentModel(options.getPartitionSubstitutionModels().get(selRow));
            frame.modelSelectionChanged(!isUsed(selRow));
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
            default:
                throw new IllegalArgumentException("Unsupported data type");

        }
        modelBorder.setTitle(title + " Substitution Model - " + currentModel.getName());
        repaint();
    }

    private boolean isUsed(int row) {
        PartitionSubstitutionModel model = options.getPartitionSubstitutionModels().get(row);
        for (PartitionData partition : options.dataPartitions) {
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
        String[] columnNames = {"Substitution Model(s)"};

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