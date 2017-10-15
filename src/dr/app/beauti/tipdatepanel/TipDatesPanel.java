/*
 * TipDatesPanel.java
 *
 * Copyright (c) 2002-2015 Alexei Drummond, Andrew Rambaut and Marc Suchard
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

package dr.app.beauti.tipdatepanel;

import dr.app.beauti.BeautiFrame;
import dr.app.beauti.BeautiPanel;
import dr.app.beauti.components.tipdatesampling.TipDateSamplingComponentOptions;
import dr.app.beauti.options.*;
import dr.app.beauti.types.TipDateSamplingType;
import dr.app.beauti.util.PanelUtils;
import dr.app.gui.table.DateCellEditor;
import dr.app.gui.table.TableEditorStopper;
import dr.app.gui.table.TableSorter;
import dr.evolution.util.*;
import dr.evolution.util.Date;
import dr.evoxml.util.DateUnitsType;
import dr.util.DataTable;
import jam.framework.Exportable;
import jam.table.HeaderRenderer;
import jam.table.TableRenderer;

import javax.swing.*;
import javax.swing.event.ListSelectionEvent;
import javax.swing.event.ListSelectionListener;
import javax.swing.filechooser.FileNameExtensionFilter;
import javax.swing.plaf.BorderUIResource;
import javax.swing.table.AbstractTableModel;
import java.awt.*;
import java.awt.event.*;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.util.*;
import java.util.List;

/**
 * @author Andrew Rambaut
 * @author Alexei Drummond
 * @author Tommy Lam
 * @version $Id: DataPanel.java,v 1.17 2006/09/05 13:29:34 rambaut Exp $
 */
public class TipDatesPanel extends BeautiPanel implements Exportable {

    /**
     *
     */
    private static final long serialVersionUID = 5283922195494563924L;

    private JScrollPane scrollPane = new JScrollPane();
    private JTable dataTable = null;
    private DataTableModel dataTableModel = null;

    private SetDatesAction setDatesAction = new SetDatesAction();
    private ClearDatesAction clearDatesAction = new ClearDatesAction();
    private GuessDatesAction guessDatesAction = new GuessDatesAction();
    private ImportDatesAction importDatesAction = new ImportDatesAction();

    private SetUncertaintyAction setUncertaintyAction = new SetUncertaintyAction();

    private JCheckBox usingTipDates = new JCheckBox("Use tip dates");
    private JCheckBox specifyOriginDate = new JCheckBox("Specify origin date:");
    private JTextField originDateText = new JTextField(20);
    private JLabel originDateLabel = new JLabel("");

    private JComboBox unitsCombo = new JComboBox(EnumSet.range(DateUnitsType.YEARS, DateUnitsType.DAYS).toArray());
    private JComboBox directionCombo = new JComboBox(EnumSet.range(DateUnitsType.FORWARDS, DateUnitsType.BACKWARDS).toArray());

    //    JComboBox tipDateSamplingCombo = new JComboBox( TipDateSamplingType.values() );
    private JComboBox tipDateSamplingCombo = new JComboBox(new TipDateSamplingType[] {
            TipDateSamplingType.NO_SAMPLING,
            TipDateSamplingType.SAMPLE_INDIVIDUALLY,
//            TipDateSamplingType.SAMPLE_JOINT,
            TipDateSamplingType.SAMPLE_PRECISION
    });
    private JComboBox tipDateTaxonSetCombo = new JComboBox();

    private BeautiFrame frame = null;

    private BeautiOptions options = null;

    private double[] heights = null;

    private GuessDatesDialog guessDatesDialog = null;
    private SetValueDialog dateValueDialog = null;
    private SetValueDialog precisionValueDialog = null;

    public TipDatesPanel(BeautiFrame parent) {

        this.frame = parent;

        dataTableModel = new DataTableModel();
        TableSorter sorter = new TableSorter(dataTableModel);
        dataTable = new JTable(sorter);

        sorter.setTableHeader(dataTable.getTableHeader());

        dataTable.getTableHeader().setReorderingAllowed(false);
//        dataTable.getTableHeader().setDefaultRenderer(
//                new HeaderRenderer(SwingConstants.LEFT, new Insets(0, 4, 0, 4)));

        dataTable.getTableHeader().setReorderingAllowed(false);
        dataTable.getTableHeader().setDefaultRenderer(
                new HeaderRenderer(SwingConstants.LEFT, new Insets(0, 4, 0, 4)));

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
        dataTable.getColumnModel().getColumn(2).setCellEditor(
                new DateCellEditor());

        dataTable.getColumnModel().getColumn(3).setCellRenderer(
                new TableRenderer(SwingConstants.LEFT, new Insets(0, 4, 0, 4)));
        dataTable.getColumnModel().getColumn(3).setPreferredWidth(80);

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
        toolBar1.setBorder(BorderFactory.createEmptyBorder());

        toolBar1.setLayout(new FlowLayout(java.awt.FlowLayout.LEFT, 0, 0));
        JButton button = new JButton(guessDatesAction);
        PanelUtils.setupComponent(button);
        toolBar1.add(button);

        button = new JButton(importDatesAction);
        PanelUtils.setupComponent(button);
        toolBar1.add(button);

        button = new JButton(setDatesAction);
        PanelUtils.setupComponent(button);
        toolBar1.add(button);

        button = new JButton(clearDatesAction);
        PanelUtils.setupComponent(button);
        toolBar1.add(button);

        button = new JButton(setUncertaintyAction);
        PanelUtils.setupComponent(button);
        toolBar1.add(button);

        JToolBar toolBar3 = new JToolBar();
        toolBar3.setFloatable(false);
        toolBar3.setOpaque(false);
        toolBar3.setBorder(BorderFactory.createEmptyBorder());

        toolBar3.setLayout(new FlowLayout(java.awt.FlowLayout.LEFT, 0, 0));

        final JLabel unitsLabel = new JLabel("Dates as:");
        toolBar3.add(unitsLabel);
        toolBar3.add(unitsCombo);
        toolBar3.add(directionCombo);

        toolBar3.add(new JToolBar.Separator(new Dimension(12, 12)));

        toolBar3.add(specifyOriginDate);
        toolBar3.add(originDateText);
        toolBar3.add(originDateLabel);

        JPanel panel2 = new JPanel();
        panel2.setLayout(new BoxLayout(panel2, BoxLayout.PAGE_AXIS));
        panel2.setOpaque(false);
        panel2.add(toolBar1);
        panel2.add(toolBar3);


        JPanel panel1 = new JPanel(new BorderLayout(0, 0));
        panel1.setOpaque(false);
        panel1.add(panel2, "North");
        panel1.add(scrollPane, "Center");

        JToolBar toolBar2 = new JToolBar();
        toolBar2.setFloatable(false);
        toolBar2.setOpaque(false);
        toolBar2.setBorder(BorderFactory.createEmptyBorder());

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

        specifyOriginDate.addItemListener(
                new java.awt.event.ItemListener() {
                    public void itemStateChanged(java.awt.event.ItemEvent ev) {
                        boolean enabled = usingTipDates.isSelected();
                        originDateText.setEnabled(enabled && specifyOriginDate.isSelected());
                        originDateLabel.setEnabled(enabled && specifyOriginDate.isSelected());
                        timeScaleChanged();
                    }
                }
        );

        originDateText.addFocusListener(
                new FocusAdapter() {
                    @Override
                    public void focusLost(FocusEvent focusEvent) {
                        timeScaleChanged();
                    }
                }
        );

        usingTipDates.addItemListener(new ItemListener() {
            @Override
            public void itemStateChanged(ItemEvent e) {
                boolean enabled = usingTipDates.isSelected();
                clearDatesAction.setEnabled(enabled);
                guessDatesAction.setEnabled(enabled);
                importDatesAction.setEnabled(enabled);
                setDatesAction.setEnabled(enabled);
                setUncertaintyAction.setEnabled(enabled);
                unitsLabel.setEnabled(enabled);
                unitsCombo.setEnabled(enabled);
                directionCombo.setEnabled(enabled);
                scrollPane.setEnabled(enabled);
                dataTable.setEnabled(enabled);
                tipDateSamplingCombo.setEnabled(enabled);
                tipDateSamplingLabel.setEnabled(enabled);
                specifyOriginDate.setEnabled(enabled);
                originDateText.setEnabled(enabled && specifyOriginDate.isSelected());
                originDateLabel.setEnabled(enabled && specifyOriginDate.isSelected());

                if (options.taxonList != null) timeScaleChanged();
            }
        });

        // because usingTipDates is listening to itemStateChanged, this will call the above action
        // to set everything disabled initially.
        usingTipDates.setSelected(false);

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

            newDate.setUncertainty(date.getUncertainty());

            options.taxonList.getTaxon(i).setDate(newDate);
        }

        if (specifyOriginDate.isSelected()) {
            String text = originDateText.getText();
            DateGuesser guesser = options.dateGuesser;
            guessDatesDialog.setupGuesser(guesser);

            try {
                options.originDate = guesser.parseDate(text);
            } catch (GuessDatesException e) {
                options.originDate = null;
            }
        } else {
            options.originDate = null;
        }
        if ( options.originDate != null) {
            originDateLabel.setText(" date value: " + Double.toString(options.originDate.getTimeValue()));
        } else {
            originDateLabel.setText(" unable to parse date");
        }

        calculateHeights();

        if (options.clockModelOptions.isTipCalibrated()) { // todo correct?
            for (PartitionTreeModel treeModel : options.getPartitionTreeModels()) {
                treeModel.setTipCalibrations(true);
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

        setupTable();
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

    private void setDates() {
        if (options.taxonList == null) { // validation of check empty taxonList
            return;
        }

        int result;
        do {
            if (dateValueDialog == null) {
                dateValueDialog = new SetValueDialog(frame, "Set Date for Taxa");
            }

            int[] selRows = dataTable.getSelectedRows();

            if (selRows.length == 1) {
                precisionValueDialog.setDescription("Set date value for selected taxon");
            } else if (selRows.length > 1) {
                dateValueDialog.setDescription("Set date values for selected taxa");
            } else {
                dateValueDialog.setDescription("Set date values for all taxa");
            }

            result = dateValueDialog.showDialog();

            if (result == -1 || result == JOptionPane.CANCEL_OPTION) {
                return;
            }

//            currentTrait.guessTrait = true; // ?? no use?
            String value = dateValueDialog.getValue();

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

    private void setPrecisions() {
        if (options.taxonList == null) { // validation of check empty taxonList
            return;
        }

        int result;
        do {
            if (precisionValueDialog == null) {
                precisionValueDialog = new SetValueDialog(frame, "Set Precision for Taxa");
            }

            int[] selRows = dataTable.getSelectedRows();

            if (selRows.length == 1) {
                precisionValueDialog.setDescription("Set precision value for selected taxon");
            } else if (selRows.length > 1) {
                precisionValueDialog.setDescription("Set precision values for selected taxa");
            } else {
                precisionValueDialog.setDescription("Set precision values for all taxa");
            }

            result = precisionValueDialog.showDialog();

            if (result == -1 || result == JOptionPane.CANCEL_OPTION) {
                return;
            }

            double value;
            try {
                value = Double.parseDouble(precisionValueDialog.getValue());
            } catch (NumberFormatException nfe) {
                value = 0.0;
            }
            if (value < 0.0 || Double.isNaN(value) || Double.isInfinite(value)) {
                value = 0.0;
            }

            if (selRows.length > 0) {
                for (int row : selRows) {
                    options.taxonList.getTaxon(row).getDate().setUncertainty(value);
                }
            } else {
                for (Taxon taxon : options.taxonList) {
                    taxon.setAttribute("precision", value);
                }
            }

            dataTableModel.fireTableDataChanged();
        } while (result < 0);
    }

    private void clearDates() {
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

    private void guessDates() {

        if (guessDatesDialog == null) {
            guessDatesDialog = new GuessDatesDialog(frame);
        }

        int[] selRows = dataTable.getSelectedRows();

        if (selRows.length > 0) {
            guessDatesDialog.setDescription("Parse date values for selected taxa");
        } else {
            guessDatesDialog.setDescription("Parse date values for all taxa");
        }


        int result = guessDatesDialog.showDialog();

        if (result == -1 || result == JOptionPane.CANCEL_OPTION) {
            return;
        }

        DateGuesser guesser = options.dateGuesser;

        guesser.guessDates = true;
        guessDatesDialog.setupGuesser(guesser);

        if (selRows.length > 0) {
            Taxa selectedTaxa = new Taxa();

            for (int row : selRows) {
                Taxon taxon = (Taxon) dataTable.getValueAt(row, 0);
                selectedTaxa.addTaxon(taxon);
            }
            guesser.guessDates(selectedTaxa);
        } else {
            guesser.guessDates(options.taxonList);
        }

        // adjust the dates to the current timescale...
        timeScaleChanged();

        dataTableModel.fireTableDataChanged();
    }

    private void importDates() {

        File[] files = frame.selectImportFiles("Import Dates File...", false, new FileNameExtensionFilter[]{
                new FileNameExtensionFilter("Tab-delimited text files", "txt", "tab", "dat")});

        DataTable<String[]> dataTable;

        if (files != null && files.length != 0) {
            try {
                // Load the file as a table
                dataTable = DataTable.Text.parse(new FileReader(files[0]), true, true, true, false);

            } catch (FileNotFoundException fnfe) {
                JOptionPane.showMessageDialog(this, "Unable to open file: File not found",
                        "Unable to open file",
                        JOptionPane.ERROR_MESSAGE);
                return;
            } catch (IOException ioe) {
                JOptionPane.showMessageDialog(this, "Unable to read file: " + ioe.getMessage(),
                        "Unable to read file",
                        JOptionPane.ERROR_MESSAGE);
                return;
            } catch (Exception ex) {
                ex.printStackTrace(System.err);
                JOptionPane.showMessageDialog(this, "Fatal exception: " + ex,
                        "Error reading file",
                        JOptionPane.ERROR_MESSAGE);
                ex.printStackTrace();
                return;
            }
        } else {
            return;
        }

        if (dataTable.getColumnCount() == 0) {
            // expecting at least 2 columns - labels and dates
            JOptionPane.showMessageDialog(frame,
                    "Expecting a tab delimited file with at\n" +
                            "least 2 columns (taxon labels and dates).",
                    "Incompatible values", JOptionPane.ERROR_MESSAGE);
            return;
        }

        String[] columnLabels = dataTable.getColumnLabels();
        String[] taxonNames = dataTable.getRowLabels();

        if (columnLabels.length < 2) {
            // only one column so leave it
            return;
        }

        boolean hasColumnHeadings = options.taxonList.getTaxonIndex(columnLabels[0]) < 0;

        // assume the second column contains the dates
        int dateColumn = 0;

        if (hasColumnHeadings && columnLabels.length > 1) {
            List<Integer> dateColumns = new ArrayList<Integer>();

            // see if there is a column labelled 'dates' or something
            for (int i = 1; i < dataTable.getColumnCount(); i++) {
                if (columnLabels[i].toLowerCase().contains("date")) {
                    dateColumns.add(i - 1);
                }
            }

            if (dateColumns.size() > 0) {
                // if there are multiple date column possibilities, take the first
                // @todo - allow the user to select the column to use
                dateColumn = dateColumns.get(0);
            }
        }

        Map<Taxon, String> taxonDateMap = new HashMap<Taxon, String>();

        String[] values = dataTable.getColumn(dateColumn);

        if (!hasColumnHeadings) {
            final int index = options.taxonList.getTaxonIndex(columnLabels[0]);
            if (index >= 0) {
                taxonDateMap.put(options.taxonList.getTaxon(index), columnLabels[dateColumn + 1]);
            }
        }

        int j = 0;
        for (final String taxonName : taxonNames) {

            final int index = options.taxonList.getTaxonIndex(taxonName);
            if (index >= 0) {
                taxonDateMap.put(options.taxonList.getTaxon(index), values[j]);
            }
            j++;
        }

        if (guessDatesDialog == null) {
            guessDatesDialog = new GuessDatesDialog(frame);
        }

        guessDatesDialog.setDescription("Parse date values from file");

        int result = guessDatesDialog.showDialog(true);

        if (result == -1 || result == JOptionPane.CANCEL_OPTION) {
            return;
        }

        DateGuesser guesser = options.dateGuesser;

        guesser.guessDates = true;
        guessDatesDialog.setupGuesser(guesser);

        guesser.guessDates(options.taxonList, taxonDateMap);

        // adjust the dates to the current timescale...
        timeScaleChanged();

        dataTableModel.fireTableDataChanged();
    }

    public class SetDatesAction extends AbstractAction {
        /**
         *
         */
        private static final long serialVersionUID = -7281309694753868635L;

        SetDatesAction() {
            super("Set Dates");
            setToolTipText("Use this tool to set sampling date values for the selected taxa");
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

        ClearDatesAction() {
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

        GuessDatesAction() {
            super("Parse Dates");
            setToolTipText("Use this tool to parse the sampling dates from the taxon labels");
        }

        public void actionPerformed(ActionEvent ae) {
            guessDates();
        }
    }

    public class ImportDatesAction extends AbstractAction {
        /**
         *
         */
        private static final long serialVersionUID = 8514706149822252033L;

        ImportDatesAction() {
            super("Import Dates");
            setToolTipText("Use this tool to import the sampling dates from a file");
        }

        public void actionPerformed(ActionEvent ae) {
            importDates();
        }
    }

    public class SetUncertaintyAction extends AbstractAction {
        /**
         *
         */
        private static final long serialVersionUID = -7281309694753868639L;

        SetUncertaintyAction() {
            super("Set Uncertainty");
            setToolTipText("Use this tool to set uncertainty in the date for the selected taxa");
        }

        public void actionPerformed(ActionEvent ae) {
            setPrecisions();
        }
    }

    private void calculateHeights() {

        options.maximumTipHeight = 0.0;
        if (options.taxonList == null || options.taxonList.getTaxonCount() == 0) return;

        heights = null;

        dr.evolution.util.Date mostRecent = options.originDate;
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
        String[] columnNames = {"Name", "Date", "Uncertainty", "Height"};

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
            Date date = options.taxonList.getTaxon(row).getDate();

            switch (col) {
                case 0:
                    return options.taxonList.getTaxon(row);
                case 1:
                    if (date != null) {
                        return date.getTimeValue();
                    } else {
                        return "-";
                    }
                case 2:
                    if (date != null) {
                        return date.getUncertainty();
                    } else {
                        return "-";
                    }
                case 3:
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
            } else if (col == 2) {
                Date date = options.taxonList.getTaxon(row).getDate();
                if (date != null) {
                    double d = (Double) aValue;
                    if (d >= 0.0) {
                        date.setUncertainty(d);
                    }
                }
            }

            timeScaleChanged();
        }

        public boolean isCellEditable(int row, int col) {
            if (col == 0) return false;
            if (col == 1 || col == 2) {
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
            StringBuilder sb = new StringBuilder();

            sb.append(getColumnName(0));
            for (int j = 1; j < getColumnCount(); j++) {
                sb.append("\t");
                sb.append(getColumnName(j));
            }
            sb.append("\n");

            for (int i = 0; i < getRowCount(); i++) {
                sb.append(getValueAt(i, 0));
                for (int j = 1; j < getColumnCount(); j++) {
                    sb.append("\t");
                    sb.append(getValueAt(i, j));
                }
                sb.append("\n");
            }

            return sb.toString();
        }
    }
}
