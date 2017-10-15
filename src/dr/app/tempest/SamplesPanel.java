/*
 * SamplesPanel.java
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

package dr.app.tempest;

import dr.app.beauti.options.DateGuesser;
import dr.app.beauti.tipdatepanel.GuessDatesDialog;
import dr.app.beauti.util.PanelUtils;
import dr.app.util.OSType;
import dr.app.util.Utils;
import dr.evolution.util.*;
import dr.app.gui.table.*;
import dr.util.DataTable;
import jam.framework.Exportable;
import jam.table.HeaderRenderer;
import dr.app.gui.table.TableEditorStopper;
import jam.table.TableRenderer;

import javax.swing.*;
import javax.swing.event.ListSelectionEvent;
import javax.swing.event.ListSelectionListener;
import javax.swing.filechooser.FileNameExtensionFilter;
import javax.swing.plaf.BorderUIResource;
import javax.swing.table.AbstractTableModel;
import java.awt.*;
import java.awt.event.ActionEvent;
import java.awt.event.ItemEvent;
import java.awt.event.ItemListener;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Map;

/**
 * @author Andrew Rambaut
 * @version $Id: DataPanel.java,v 1.17 2006/09/05 13:29:34 rambaut Exp $
 */
public class SamplesPanel extends JPanel implements Exportable {

    /**
     *
     */
    private static final long serialVersionUID = 5283922195494563924L;
    JScrollPane scrollPane = new JScrollPane();
    JTable dataTable = null;
    DataTableModel dataTableModel = null;

    ClearDatesAction clearDatesAction = new ClearDatesAction();
    GuessDatesAction guessDatesAction = new GuessDatesAction();
    ImportDatesAction importDatesAction = new ImportDatesAction();

    JComboBox unitsCombo = new JComboBox(new String[]{"Years", "Months", "Days"});
    JComboBox directionCombo = new JComboBox(new String[]{"Since some time in the past", "Before the present"});

    TempestFrame frame = null;

    TaxonList taxonList = null;

    int datesUnits;
    int datesDirection;
    double maximumTipHeight = 0.0;

    DateGuesser guesser = new DateGuesser();

    double[] heights = null;

    GuessDatesDialog guessDatesDialog = null;

    public SamplesPanel(TempestFrame parent, TaxonList taxonList) {

        this.frame = parent;

        frame.setImportAction(importDatesAction);
        dataTableModel = new DataTableModel();
        TableSorter sorter = new TableSorter(dataTableModel);
        dataTable = new JTable(sorter);

        sorter.setTableHeader(dataTable.getTableHeader());

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

        toolBar1.setLayout(new FlowLayout(FlowLayout.LEFT, 0, 0));
        JButton button = new JButton(importDatesAction);
        PanelUtils.setupComponent(button);
        toolBar1.add(button);
         button = new JButton(guessDatesAction);
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

        setOpaque(false);
        setBorder(new BorderUIResource.EmptyBorderUIResource(new Insets(12, 12, 12, 12)));
        setLayout(new BorderLayout(0, 0));

        add(toolBar1, "North");
        add(scrollPane, "Center");

        ItemListener listener = new ItemListener() {
            public void itemStateChanged(ItemEvent ev) {
                timeScaleChanged();
            }
        };
        unitsCombo.addItemListener(listener);
        directionCombo.addItemListener(listener);

        setTaxonList(taxonList);
    }

    public final void timeScaleChanged() {
        Units.Type units = Units.Type.YEARS;
        switch (unitsCombo.getSelectedIndex()) {
            case 0:
                units = Units.Type.YEARS;
                break;
            case 1:
                units = Units.Type.MONTHS;
                break;
            case 2:
                units = Units.Type.DAYS;
                break;
        }

        boolean backwards = directionCombo.getSelectedIndex() == 1;

        for (int i = 0; i < taxonList.getTaxonCount(); i++) {
            Date date = taxonList.getTaxon(i).getDate();
            double d = date.getTimeValue();

            Date newDate = createDate(d, units, backwards, 0.0);

            newDate.setUncertainty(date.getUncertainty());

            taxonList.getTaxon(i).setDate(newDate);
        }

        calculateHeights();

        dataTableModel.fireTableDataChanged();
        frame.timeScaleChanged();
    }

    private Date createDate(double timeValue, Units.Type units, boolean backwards, double origin) {
        if (backwards) {
            return Date.createTimeAgoFromOrigin(timeValue, units, origin);
        } else {
            return Date.createTimeSinceOrigin(timeValue, units, origin);
        }
    }

    private void setTaxonList(TaxonList taxonList) {
        this.taxonList = taxonList;

        setupTable();

        unitsCombo.setSelectedIndex(datesUnits);
        directionCombo.setSelectedIndex(datesDirection);

        calculateHeights();

        dataTableModel.fireTableDataChanged();
    }

    private void setupTable() {
        dataTableModel.fireTableDataChanged();
    }

    public void getOptions() {
        datesUnits = unitsCombo.getSelectedIndex();
        datesDirection = directionCombo.getSelectedIndex();
    }

    public JComponent getExportableComponent() {
        return dataTable;
    }

    public void selectionChanged() {
        // nothing to do
    }

    public void clearDates() {
        for (int i = 0; i < taxonList.getTaxonCount(); i++) {
            java.util.Date origin = new java.util.Date(0);

            double d = 0.0;

            Date date = Date.createTimeSinceOrigin(d, Units.Type.YEARS, origin);
            taxonList.getTaxon(i).setAttribute("date", date);
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

        guesser.guessDates = true;
        guessDatesDialog.setupGuesser(guesser);

        String warningMessage = null;

        guesser.guessDates(taxonList);

        if (warningMessage != null) {
            JOptionPane.showMessageDialog(this, "Warning: some dates may not be set correctly - \n" + warningMessage,
                    "Error guessing dates",
                    JOptionPane.WARNING_MESSAGE);
        }

        // adjust the dates to the current timescale...
        timeScaleChanged();

        dataTableModel.fireTableDataChanged();
    }

    private Map<String, FileDialog> fileDialogs = new HashMap<String, FileDialog>();
    private Map<String, JFileChooser> fileChoosers = new HashMap<String, JFileChooser>();
    /**
     * Use the native file dialog on the Mac because the Swing one is bad. On linux, the native
     * one is bad. No preference on Windows.
     * @param title
     * @return
     */
    public File[] selectImportFiles(final String title, boolean multipleSelection, FileNameExtensionFilter[] fileNameExtensionFilters) {
        if (Boolean.parseBoolean(System.getProperty("use.native.choosers", Boolean.toString(OSType.isMac())))) {
            FileDialog importDialog = fileDialogs.get(title);
            if (importDialog == null) {
                importDialog = new FileDialog(frame, title, FileDialog.LOAD);
                fileDialogs.put(title, importDialog);
            }

            importDialog.setVisible(true);
            if (importDialog.getFile() != null) {
                return new File[] { new File(importDialog.getDirectory(), importDialog.getFile()) };
            }
        } else {
            JFileChooser importChooser = fileChoosers.get(title);
            if (importChooser == null) {
                importChooser = new JFileChooser(Utils.getCWD());

                importChooser.setMultiSelectionEnabled(multipleSelection);
                for (FileNameExtensionFilter fileNameExtensionFilter : fileNameExtensionFilters) {
                    importChooser.setFileFilter(fileNameExtensionFilter);
                }
                importChooser.setDialogTitle(title);

                fileChoosers.put(title, importChooser);
            }

            int returnVal = importChooser.showOpenDialog(this);
            if (returnVal == JFileChooser.APPROVE_OPTION) {
                if (importChooser.isMultiSelectionEnabled()) {
                    return importChooser.getSelectedFiles();
                } else {
                    return new File[] { importChooser.getSelectedFile() };
                }
            }
        }

        return null;
    }

    private void importDates() {

        File[] files = selectImportFiles("Import Dates File...", false, new FileNameExtensionFilter[]{
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

        boolean hasColumnHeadings = taxonList.getTaxonIndex(columnLabels[0]) < 0;

        // assume the second column contains the dates
        int dateColumn = 0;

        if (hasColumnHeadings && columnLabels.length > 1) {
            java.util.List<Integer> dateColumns = new ArrayList<Integer>();

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
            final int index = taxonList.getTaxonIndex(columnLabels[0]);
            if (index >= 0) {
                taxonDateMap.put(taxonList.getTaxon(index), columnLabels[dateColumn + 1]);
            }
        }

        int j = 0;
        for (final String taxonName : taxonNames) {

            final int index = taxonList.getTaxonIndex(taxonName);
            if (index >= 0) {
                taxonDateMap.put(taxonList.getTaxon(index), values[j]);
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

        guesser.guessDates = true;
        guessDatesDialog.setupGuesser(guesser);

        guesser.guessDates(taxonList, taxonDateMap);

        // adjust the dates to the current timescale...
        timeScaleChanged();

        dataTableModel.fireTableDataChanged();
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

    private void calculateHeights() {

        maximumTipHeight = 0.0;
        if (taxonList == null || taxonList.getTaxonCount() == 0) return;

        heights = null;

        Date mostRecent = null;
        for (int i = 0; i < taxonList.getTaxonCount(); i++) {
            Date date = taxonList.getTaxon(i).getDate();
            if ((date != null) && (mostRecent == null || date.after(mostRecent))) {
                mostRecent = date;
            }
        }

        if (mostRecent != null) {
            heights = new double[taxonList.getTaxonCount()];

            TimeScale timeScale = new TimeScale(mostRecent.getUnits(), true, mostRecent.getAbsoluteTimeValue());
            double time0 = timeScale.convertTime(mostRecent.getTimeValue(), mostRecent);

            for (int i = 0; i < taxonList.getTaxonCount(); i++) {
                Date date = taxonList.getTaxon(i).getDate();
                if (date != null) {
                    heights[i] = timeScale.convertTime(date.getTimeValue(), date) - time0;
                    if (heights[i] > maximumTipHeight) maximumTipHeight = heights[i];
                }
            }
        }
    }

    class DataTableModel extends AbstractTableModel {

        /**
         *
         */
        private static final long serialVersionUID = -6707994233020715574L;
        String[] columnNames = {"Name", "Date", "Precision", "Height"};

        public DataTableModel() {
        }

        public int getColumnCount() {
            return columnNames.length;
        }

        public int getRowCount() {
            if (taxonList == null) return 0;

            return taxonList.getTaxonCount();
        }

        public Object getValueAt(int row, int col) {
            Date date = taxonList.getTaxon(row).getDate();
            switch (col) {
                case 0:
                    return taxonList.getTaxonId(row);
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
                taxonList.getTaxon(row).setId(aValue.toString());
            } else if (col == 1) {
                Date date = taxonList.getTaxon(row).getDate();
                if (date != null) {
                    double d = (Double) aValue;
                    Date newDate = createDate(d, date.getUnits(), date.isBackwards(), date.getOrigin());
                    taxonList.getTaxon(row).setDate(newDate);
                }
            } else if (col == 2) {
                Date date = taxonList.getTaxon(row).getDate();
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
            if (col == 0) return true;
            if (col == 1 || col == 2) {
                Date date = taxonList.getTaxon(row).getDate();
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
