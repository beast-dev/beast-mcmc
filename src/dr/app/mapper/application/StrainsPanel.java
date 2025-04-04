/*
 * StrainsPanel.java
 *
 * Copyright © 2002-2024 the BEAST Development Team
 * http://beast.community/about
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
 *
 */

package dr.app.mapper.application;

import dr.app.beauti.options.DateGuesser;
import dr.app.beauti.tipdatepanel.GuessDatesDialog;
import dr.app.beauti.util.PanelUtils;
import dr.app.gui.table.DateCellEditor;
import dr.app.gui.table.TableEditorStopper;
import dr.app.gui.table.TableSorter;
import dr.evolution.util.*;
import jam.framework.Exportable;
import jam.table.HeaderRenderer;
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

/**
 * @author Andrew Rambaut
 */
public class StrainsPanel extends JPanel implements Exportable, MapperDocument.Listener {

    private JScrollPane scrollPane = new JScrollPane();
    private JTable dataTable = null;
    private DataTableModel dataTableModel = null;

    private ClearDatesAction clearDatesAction = new ClearDatesAction();
    private GuessDatesAction guessDatesAction = new GuessDatesAction();

    private JComboBox unitsCombo = new JComboBox(new String[]{"Years", "Months", "Days"});
    private JComboBox directionCombo = new JComboBox(new String[]{"Since some time in the past", "Before the present"});

    private final MapperFrame frame;
    private final MapperDocument document;

    int datesUnits;
    int datesDirection;
    double maximumTipHeight = 0.0;

    DateGuesser guesser = new DateGuesser();

    double[] heights = null;

    GuessDatesDialog guessDatesDialog = null;

    public StrainsPanel(final MapperFrame parent, final MapperDocument document) {

        this.frame = parent;
        this.document = document;

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
        JButton button = new JButton(clearDatesAction);
        PanelUtils.setupComponent(button);
        toolBar1.add(button);
        button = new JButton(guessDatesAction);
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

        java.util.List<Taxon> taxonList = document.getTaxa();

        for (int i = 0; i < taxonList.size(); i++) {
            Date date = taxonList.get(i).getDate();
            double d = date.getTimeValue();

            Date newDate = createDate(d, units, backwards, 0.0);

            taxonList.get(i).setDate(newDate);
        }

        calculateHeights();

        dataTableModel.fireTableDataChanged();
        document.fireTaxaChanged();
    }

    private Date createDate(double timeValue, Units.Type units, boolean backwards, double origin) {
        if (backwards) {
            return Date.createTimeAgoFromOrigin(timeValue, units, origin);
        } else {
            return Date.createTimeSinceOrigin(timeValue, units, origin);
        }
    }

    @Override
    public void taxaChanged() {
        unitsCombo.setSelectedIndex(datesUnits);
        directionCombo.setSelectedIndex(datesDirection);

        calculateHeights();

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
        java.util.List<Taxon> taxonList = document.getTaxa();

        for (int i = 0; i < taxonList.size(); i++) {
            java.util.Date origin = new java.util.Date(0);

            double d = 0.0;

            Date date = Date.createTimeSinceOrigin(d, Units.Type.YEARS, origin);
            taxonList.get(i).setAttribute("date", date);
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

        guesser.guessDates(document.getTaxa());

        if (warningMessage != null) {
            JOptionPane.showMessageDialog(this, "Warning: some dates may not be set correctly - \n" + warningMessage,
                    "Error guessing dates",
                    JOptionPane.WARNING_MESSAGE);
        }

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
            super("Guess Dates");
            setToolTipText("Use this tool to guess the sampling dates from the taxon labels");
        }

        public void actionPerformed(ActionEvent ae) {
            guessDates();
        }
    }

    private void calculateHeights() {
        java.util.List<Taxon> taxonList = document.getTaxa();

        maximumTipHeight = 0.0;
        if (taxonList == null || taxonList.size() == 0) return;

        heights = null;

        Date mostRecent = null;
        for (Taxon taxon : taxonList) {
            Date date = taxon.getDate();
            if ((date != null) && (mostRecent == null || date.after(mostRecent))) {
                mostRecent = date;
            }
        }

        if (mostRecent != null) {
            heights = new double[taxonList.size()];

            TimeScale timeScale = new TimeScale(mostRecent.getUnits(), true, mostRecent.getAbsoluteTimeValue());
            double time0 = timeScale.convertTime(mostRecent.getTimeValue(), mostRecent);

            for (int i = 0; i < taxonList.size(); i++) {
                Date date = taxonList.get(i).getDate();
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
        String[] columnNames = {"Name", "Date", "Height"};

        public DataTableModel() {
        }

        public int getColumnCount() {
            return columnNames.length;
        }

        public int getRowCount() {
            java.util.List<Taxon> taxonList = document.getTaxa();

            if (taxonList == null) return 0;

            return taxonList.size();
        }

        public Object getValueAt(int row, int col) {
            java.util.List<Taxon> taxonList = document.getTaxa();

            switch (col) {
                case 0:
                    return taxonList.get(row);
                case 1:
                    Date date = taxonList.get(row).getDate();
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
            java.util.List<Taxon> taxonList = document.getTaxa();

            if (col == 0) {
                taxonList.get(row).setId(aValue.toString());
            } else if (col == 1) {
                Date date = taxonList.get(row).getDate();
                if (date != null) {
                    double d = (Double) aValue;
                    Date newDate = createDate(d, date.getUnits(), date.isBackwards(), date.getOrigin());
                    taxonList.get(row).setDate(newDate);
                }
            }

            timeScaleChanged();
        }

        public boolean isCellEditable(int row, int col) {
            if (col == 0) return true;
            if (col == 1) {
                Date date = document.getTaxa().get(row).getDate();
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
