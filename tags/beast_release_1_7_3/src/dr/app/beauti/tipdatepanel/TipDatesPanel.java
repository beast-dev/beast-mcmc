/*
 * TipDatesPanel.java
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

package dr.app.beauti.tipdatepanel;

import dr.app.beauti.BeautiFrame;
import dr.app.beauti.BeautiPanel;
import dr.app.beauti.components.tipdatesampling.TipDateSamplingComponentOptions;
import dr.app.beauti.options.BeautiOptions;
import dr.app.beauti.options.ClockModelGroup;
import dr.app.beauti.options.DateGuesser;
import dr.app.beauti.traitspanel.TraitValueDialog;
import dr.app.beauti.types.TipDateSamplingType;
import dr.app.beauti.util.PanelUtils;
import dr.app.gui.table.DateCellEditor;
import dr.app.gui.table.TableEditorStopper;
import dr.app.gui.table.TableSorter;
import dr.evolution.util.*;
import dr.evoxml.util.DateUnitsType;
import jam.framework.Exportable;
import jam.table.TableRenderer;

import javax.swing.*;
import javax.swing.event.ListSelectionEvent;
import javax.swing.event.ListSelectionListener;
import javax.swing.plaf.BorderUIResource;
import javax.swing.table.AbstractTableModel;
import java.awt.*;
import java.awt.event.ActionEvent;
import java.awt.event.ItemEvent;
import java.awt.event.ItemListener;
import java.util.EnumSet;

/**
 * @author Andrew Rambaut
 * @author Alexei Drummond
 * @version $Id: DataPanel.java,v 1.17 2006/09/05 13:29:34 rambaut Exp $
 */
public class TipDatesPanel extends BeautiPanel implements Exportable {

    /**
     *
     */
    private static final long serialVersionUID = 5283922195494563924L;

    JScrollPane scrollPane = new JScrollPane();
    JTable dataTable = null;
    DataTableModel dataTableModel = null;

    SetDatesAction setDatesAction = new SetDatesAction();
    ClearDatesAction clearDatesAction = new ClearDatesAction();
    GuessDatesAction guessDatesAction = new GuessDatesAction();

    JCheckBox usingTipDates = new JCheckBox("Use tip dates");

    JComboBox unitsCombo = new JComboBox(EnumSet.range(DateUnitsType.YEARS, DateUnitsType.DAYS).toArray());
    JComboBox directionCombo = new JComboBox(EnumSet.range(DateUnitsType.FORWARDS, DateUnitsType.BACKWARDS).toArray());

//    JComboBox tipDateSamplingCombo = new JComboBox( TipDateSamplingType.values() );
    JComboBox tipDateSamplingCombo = new JComboBox(new TipDateSamplingType[] { TipDateSamplingType.NO_SAMPLING, TipDateSamplingType.SAMPLE_INDIVIDUALLY });
    JComboBox tipDateTaxonSetCombo = new JComboBox();

    BeautiFrame frame = null;

    BeautiOptions options = null;

    double[] heights = null;

    private GuessDatesDialog guessDatesDialog = null;
    private DateValueDialog dateValueDialog = null;

    public TipDatesPanel(BeautiFrame parent) {

        this.frame = parent;

        dataTableModel = new DataTableModel();
        TableSorter sorter = new TableSorter(dataTableModel);
        dataTable = new JTable(sorter);

        sorter.setTableHeader(dataTable.getTableHeader());

        dataTable.getTableHeader().setReorderingAllowed(false);
//        dataTable.getTableHeader().setDefaultRenderer(
//                new HeaderRenderer(SwingConstants.LEFT, new Insets(0, 4, 0, 4)));

        dataTable.getColumnModel().getColumn(0).setCellRenderer(
                new TableRenderer(SwingConstants.LEFT, new Insets(0, 4, 0, 4)));
        dataTable.getColumnModel().getColumn(0).setPreferredWidth(80);

        dataTable.getColumnModel().getColumn(1).setCellRenderer(
                new TableRenderer(SwingConstants.LEFT, new Insets(0, 4, 0, 4)));
        dataTable.getColumnModel().getColumn(1).setPreferredWidth(80);
        dataTable.getColumnModel().getColumn(1).setCellEditor(
                new DateCellEditor());

        dataTable.getColumnModel().getColumn(2).setCellRenderer(
                new TableRenderer(SwingConstants.LEFT, new Insets(0, 4, 0, 4)));
        dataTable.getColumnModel().getColumn(2).setPreferredWidth(80);

        TableEditorStopper.ensureEditingStopWhenTableLosesFocus(dataTable);

        dataTable.getSelectionModel().addListSelectionListener(new ListSelectionListener() {
            public void valueChanged(ListSelectionEvent evt) {
                selectionChanged();
            }
        });

        scrollPane = new JScrollPane(dataTable,
                JScrollPane.VERTICAL_SCROLLBAR_ALWAYS,
                JScrollPane.HORIZONTAL_SCROLLBAR_ALWAYS);
        scrollPane.setOpaque(false);

        PanelUtils.setupComponent(unitsCombo);
        PanelUtils.setupComponent(directionCombo);

        JToolBar toolBar1 = new JToolBar();
        toolBar1.setFloatable(false);
        toolBar1.setOpaque(false);

        toolBar1.setLayout(new FlowLayout(java.awt.FlowLayout.LEFT, 0, 0));
        JButton button = new JButton(guessDatesAction);
        PanelUtils.setupComponent(button);
        toolBar1.add(button);

        button = new JButton(setDatesAction);
        PanelUtils.setupComponent(button);
        toolBar1.add(button);

        button = new JButton(clearDatesAction);
        PanelUtils.setupComponent(button);
        toolBar1.add(button);

        toolBar1.add(new JToolBar.Separator(new Dimension(12, 12)));
        final JLabel unitsLabel = new JLabel("Dates specified as ");
        toolBar1.add(unitsLabel);
        toolBar1.add(unitsCombo);
        toolBar1.add(directionCombo);

        JPanel panel1 = new JPanel(new BorderLayout(0, 0));
        panel1.setOpaque(false);
        panel1.add(toolBar1, "North");
        panel1.add(scrollPane, "Center");

        JToolBar toolBar2 = new JToolBar();
        toolBar2.setFloatable(false);
        toolBar2.setOpaque(false);

        toolBar2.setLayout(new FlowLayout(java.awt.FlowLayout.LEFT, 0, 0));

        PanelUtils.setupComponent(tipDateSamplingCombo);
        tipDateSamplingCombo.setToolTipText("<html>Select whether to allow sampling<br>" +
                "of all or individual tip dates.</html>");

//        substitutionRateField.setToolTipText("<html>Enter the substitution rate here.</html>");
//        substitutionRateField.setEnabled(true);

        final JLabel tipDateSamplingLabel = new JLabel("Tip date sampling:");
        toolBar2.add(tipDateSamplingLabel);
        toolBar2.add(tipDateSamplingCombo);

        final JLabel tipDateTaxonSetLabel = new JLabel("Apply to taxon set:");
        toolBar2.add(tipDateTaxonSetLabel);
        toolBar2.add(tipDateTaxonSetCombo);


        setOpaque(false);
        setBorder(new BorderUIResource.EmptyBorderUIResource(new java.awt.Insets(12, 12, 12, 12)));
        setLayout(new BorderLayout(0, 0));
        add(usingTipDates, BorderLayout.NORTH);
        add(panel1, BorderLayout.CENTER);
        add(toolBar2, BorderLayout.SOUTH);

        tipDateSamplingCombo.addItemListener(
                new java.awt.event.ItemListener() {
                    public void itemStateChanged(java.awt.event.ItemEvent ev) {
                        boolean samplingOn = tipDateSamplingCombo.getSelectedItem() != TipDateSamplingType.NO_SAMPLING;
                        tipDateTaxonSetLabel.setEnabled(samplingOn);
                        tipDateTaxonSetCombo.setEnabled(samplingOn);
                        fireModelsChanged();
                    }
                }
        );

        clearDatesAction.setEnabled(false);
        guessDatesAction.setEnabled(false);
        directionCombo.setEnabled(false);
        unitsLabel.setEnabled(false);
        unitsCombo.setEnabled(false);
        scrollPane.setEnabled(false);
        dataTable.setEnabled(false);
        tipDateSamplingLabel.setEnabled(false);
        tipDateSamplingCombo.setEnabled(false);
        tipDateTaxonSetLabel.setEnabled(false);
        tipDateTaxonSetCombo.setEnabled(false);

        usingTipDates.addItemListener(new ItemListener() {
            public void itemStateChanged(ItemEvent ev) {
                boolean enabled = usingTipDates.isSelected();
                clearDatesAction.setEnabled(enabled);
                guessDatesAction.setEnabled(enabled);
                unitsLabel.setEnabled(enabled);
                unitsCombo.setEnabled(enabled);
                directionCombo.setEnabled(enabled);
                scrollPane.setEnabled(enabled);
                dataTable.setEnabled(enabled);
                tipDateSamplingCombo.setEnabled(enabled);
                tipDateSamplingLabel.setEnabled(enabled);

                if (options.taxonList != null) timeScaleChanged();
            }
        });

        ItemListener listener = new ItemListener() {
            public void itemStateChanged(ItemEvent ev) {
                timeScaleChanged();
            }
        };
        unitsCombo.addItemListener(listener);
        directionCombo.addItemListener(listener);
    }

    public final void timeScaleChanged() {
        Units.Type units = Units.Type.YEARS;
        switch ((DateUnitsType) unitsCombo.getSelectedItem()) {
            case YEARS:
                units = Units.Type.YEARS;
                break;
            case MONTHS:
                units = Units.Type.MONTHS;
                break;
            case DAYS:
                units = Units.Type.DAYS;
                break;
        }

        boolean backwards = directionCombo.getSelectedItem() == DateUnitsType.BACKWARDS;

        for (int i = 0; i < options.taxonList.getTaxonCount(); i++) {
            Date date = options.taxonList.getTaxon(i).getDate();
            double d = date.getTimeValue();

            Date newDate = createDate(d, units, backwards, 0.0);

            options.taxonList.getTaxon(i).setDate(newDate);
        }

        calculateHeights();

        if (options.clockModelOptions.isTipCalibrated()) { // todo correct?
            for (ClockModelGroup clockModelGroup : options.clockModelOptions.getClockModelGroups()) {
                options.clockModelOptions.tipTimeCalibration(clockModelGroup);
            }
        }

        dataTableModel.fireTableDataChanged();
        frame.setDirty();
    }

    private Date createDate(double timeValue, Units.Type units, boolean backwards, double origin) {
        if (backwards) {
            return Date.createTimeAgoFromOrigin(timeValue, units, origin);
        } else {
            return Date.createTimeSinceOrigin(timeValue, units, origin);
        }
    }

    public void setOptions(BeautiOptions options) {
        this.options = options;

        setupTable();

        unitsCombo.setSelectedItem(options.datesUnits);
        directionCombo.setSelectedItem(options.datesDirection);

        calculateHeights();
        usingTipDates.setSelected(options.clockModelOptions.isTipCalibrated());

        dataTableModel.fireTableDataChanged();

        Object item = tipDateTaxonSetCombo.getSelectedItem();

        tipDateTaxonSetCombo.removeAllItems();
        tipDateTaxonSetCombo.addItem("All taxa");
        for (TaxonList taxa : options.taxonSets) {
            tipDateTaxonSetCombo.addItem(taxa);
        }

        if (item != null) {
            tipDateTaxonSetCombo.setSelectedItem(item);
        }
    }

    private void setupTable() {
        dataTableModel.fireTableDataChanged();
    }

    public void getOptions(BeautiOptions options) {
        options.datesUnits = (DateUnitsType) unitsCombo.getSelectedItem();
        options.datesDirection = (DateUnitsType) directionCombo.getSelectedItem();

        TipDateSamplingComponentOptions comp = (TipDateSamplingComponentOptions) options.getComponentOptions(TipDateSamplingComponentOptions.class);
        comp.tipDateSamplingType = (TipDateSamplingType) tipDateSamplingCombo.getSelectedItem();
        if (tipDateTaxonSetCombo.getSelectedItem() instanceof TaxonList) {
            comp.tipDateSamplingTaxonSet = (TaxonList) tipDateTaxonSetCombo.getSelectedItem();
        } else {
            comp.tipDateSamplingTaxonSet = null;
        }
    }

    private void fireModelsChanged() {
        frame.setDirty();
    }

    public JComponent getExportableComponent() {
        return dataTable;
    }

    public void selectionChanged() {
        // nothing to do
    }

    public void setDates() {
        if (options.taxonList == null) { // validation of check empty taxonList
            return;
        }

        int result;
        do {
            if (dateValueDialog == null) {
                dateValueDialog = new DateValueDialog(frame);
            }

            int[] selRows = dataTable.getSelectedRows();

            if (selRows.length > 0) {
                dateValueDialog.setDescription("Set date values for selected taxa");
            } else {
                dateValueDialog.setDescription("Set date values for all taxa");
            }

            result = dateValueDialog.showDialog();

            if (result == -1 || result == JOptionPane.CANCEL_OPTION) {
                return;
            }

//            currentTrait.guessTrait = true; // ?? no use?
            String value = dateValueDialog.getDateValue();

            java.util.Date origin = new java.util.Date(0);
            double d = Double.parseDouble(value);
            Date date = Date.createTimeSinceOrigin(d, Units.Type.YEARS, origin);

            if (selRows.length > 0) {
                for (int row : selRows) {
                    options.taxonList.getTaxon(row).setAttribute("date", date);
                }
            } else {
                for (Taxon taxon : options.taxonList) {
                    taxon.setAttribute("date", date);
                }
            }

            // adjust the dates to the current timescale...
            timeScaleChanged();

           dataTableModel.fireTableDataChanged();
        } while (result < 0);
    }

    public void clearDates() {
        for (int i = 0; i < options.taxonList.getTaxonCount(); i++) {
            java.util.Date origin = new java.util.Date(0);

            double d = 0.0;

            Date date = Date.createTimeSinceOrigin(d, Units.Type.YEARS, origin);
            options.taxonList.getTaxon(i).setAttribute("date", date);
        }

        // adjust the dates to the current timescale...
        timeScaleChanged();

        dataTableModel.fireTableDataChanged();
    }

    public void guessDates() {

        if (guessDatesDialog == null) {
            guessDatesDialog = new GuessDatesDialog(frame);
        }

        int result = guessDatesDialog.showDialog();

        if (result == -1 || result == JOptionPane.CANCEL_OPTION) {
            return;
        }

        DateGuesser guesser = options.dateGuesser;

        guesser.guessDates = true;
        guessDatesDialog.setupGuesser(guesser);

        String warningMessage = null;

        guesser.guessDates(options.taxonList);

        if (warningMessage != null) {
            JOptionPane.showMessageDialog(this, "Warning: some dates may not be set correctly - \n" + warningMessage,
                    "Error guessing dates",
                    JOptionPane.WARNING_MESSAGE);
        }

        // adjust the dates to the current timescale...
        timeScaleChanged();

        dataTableModel.fireTableDataChanged();
    }

    public class SetDatesAction extends AbstractAction {
        /**
         *
         */
        private static final long serialVersionUID = -7281309694753868635L;

        public SetDatesAction() {
            super("Set Dates");
            setToolTipText("Use this tool to set sampling date values from selected taxa");
        }

        public void actionPerformed(ActionEvent ae) {
            setDates();
        }
    }

    public class ClearDatesAction extends AbstractAction {
        /**
         *
         */
        private static final long serialVersionUID = -7281309694753868635L;

        public ClearDatesAction() {
            super("Clear Dates");
            setToolTipText("Use this tool to remove sampling dates from each taxon");
        }

        public void actionPerformed(ActionEvent ae) {
            clearDates();
        }
    }

    public class GuessDatesAction extends AbstractAction {
        /**
         *
         */
        private static final long serialVersionUID = 8514706149822252033L;

        public GuessDatesAction() {
            super("Guess Dates");
            setToolTipText("Use this tool to guess the sampling dates from the taxon labels");
        }

        public void actionPerformed(ActionEvent ae) {
            guessDates();
        }
    }

    private void calculateHeights() {

        options.maximumTipHeight = 0.0;
        if (options.taxonList == null || options.taxonList.getTaxonCount() == 0) return;

        heights = null;

        dr.evolution.util.Date mostRecent = null;
        for (int i = 0; i < options.taxonList.getTaxonCount(); i++) {
            Date date = options.taxonList.getTaxon(i).getDate();
            if ((date != null) && (mostRecent == null || date.after(mostRecent))) {
                mostRecent = date;
            }
        }

        if (mostRecent != null) {
            heights = new double[options.taxonList.getTaxonCount()];

            TimeScale timeScale = new TimeScale(mostRecent.getUnits(), true, mostRecent.getAbsoluteTimeValue());
            double time0 = timeScale.convertTime(mostRecent.getTimeValue(), mostRecent);

            for (int i = 0; i < options.taxonList.getTaxonCount(); i++) {
                Taxon taxon = options.taxonList.getTaxon(i);
                Date date = taxon.getDate();
                if (date != null) {
                    heights[i] = timeScale.convertTime(date.getTimeValue(), date) - time0;
                    taxon.setAttribute("height", heights[i]);
                    if (heights[i] > options.maximumTipHeight) options.maximumTipHeight = heights[i];
                }
            }
        }

        frame.setStatusMessage();
    }

    class DataTableModel extends AbstractTableModel {

        /**
         *
         */
        private static final long serialVersionUID = -6707994233020715574L;
        String[] columnNames = {"Name", "Date", "Height"};

        public DataTableModel() {
        }

        public int getColumnCount() {
            return columnNames.length;
        }

        public int getRowCount() {
            if (options == null) return 0;
            if (options.taxonList == null) return 0;

            return options.taxonList.getTaxonCount();
        }

        public Object getValueAt(int row, int col) {
            switch (col) {
                case 0:
                    return options.taxonList.getTaxonId(row);
                case 1:
                    Date date = options.taxonList.getTaxon(row).getDate();
                    if (date != null) {
                        return date.getTimeValue();
                    } else {
                        return "-";
                    }
                case 2:
                    if (heights != null) {
                        return heights[row];
                    } else {
                        return "0.0";
                    }
            }
            return null;
        }

        public void setValueAt(Object aValue, int row, int col) {
            if (col == 0) {
                options.taxonList.getTaxon(row).setId(aValue.toString());
            } else if (col == 1) {
                Date date = options.taxonList.getTaxon(row).getDate();
                if (date != null) {
                    double d = (Double) aValue;
                    Date newDate = createDate(d, date.getUnits(), date.isBackwards(), date.getOrigin());
                    options.taxonList.getTaxon(row).setDate(newDate);
                }
            }

            timeScaleChanged();
        }

        public boolean isCellEditable(int row, int col) {
            if (col == 0) return false;
            if (col == 1) {
                Date date = options.taxonList.getTaxon(row).getDate();
                return (date != null);
            }
            return false;
        }

        public String getColumnName(int column) {
            return columnNames[column];
        }

        public Class getColumnClass(int c) {
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
