/*
 * SamplesPanel.java
 *
 * Copyright (C) 2002-2006 Alexei Drummond and Andrew Rambaut
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

package dr.app.beauti;

import dr.app.beauti.options.BeautiOptions;
import dr.app.beauti.options.DateGuesser;
import dr.evolution.util.Date;
import dr.evolution.util.TimeScale;
import dr.evolution.util.Units;
import dr.gui.table.DateCellEditor;
import dr.gui.table.TableSorter;
import org.virion.jam.components.RealNumberField;
import org.virion.jam.framework.Exportable;
import org.virion.jam.panels.OptionsPanel;
import org.virion.jam.table.HeaderRenderer;
import org.virion.jam.table.TableEditorStopper;
import org.virion.jam.table.TableRenderer;

import javax.swing.*;
import javax.swing.border.EmptyBorder;
import javax.swing.event.ListSelectionEvent;
import javax.swing.event.ListSelectionListener;
import javax.swing.plaf.BorderUIResource;
import javax.swing.table.AbstractTableModel;
import java.awt.*;
import java.awt.event.ActionEvent;
import java.awt.event.ItemEvent;
import java.awt.event.ItemListener;
import java.util.Calendar;
import java.util.GregorianCalendar;

/**
 * @author Andrew Rambaut
 * @author Alexei Drummond
 * @version $Id: DataPanel.java,v 1.17 2006/09/05 13:29:34 rambaut Exp $
 */
public class SamplesPanel extends BeautiPanel implements Exportable {

    /**
     *
     */
    private static final long serialVersionUID = 5283922195494563924L;
    JScrollPane scrollPane = new JScrollPane();
    JTable dataTable = null;
    DataTableModel dataTableModel = null;

    ClearDatesAction clearDatesAction = new ClearDatesAction();
    GuessDatesAction guessDatesAction = new GuessDatesAction();

    JCheckBox usingSamplingDates  = new JCheckBox("Use sampling dates");

    JComboBox unitsCombo = new JComboBox(new String[]{"Years", "Months", "Days"});
    JComboBox directionCombo = new JComboBox(new String[]{"Since some time in the past", "Before the present"});

    BeautiFrame frame = null;

    BeautiOptions options = null;

    double[] heights = null;

    public SamplesPanel(BeautiFrame parent) {

        this.frame = parent;

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

        toolBar1.setLayout(new FlowLayout(java.awt.FlowLayout.LEFT, 0, 0));
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

        JPanel panel1 = new JPanel(new BorderLayout(0, 0));
        panel1.setOpaque(false);
        panel1.add(toolBar1, "North");
        panel1.add(scrollPane, "Center");

        setOpaque(false);
        setBorder(new BorderUIResource.EmptyBorderUIResource(new java.awt.Insets(12, 12, 12, 12)));
        setLayout(new BorderLayout(0, 0));
        add(usingSamplingDates, "North");
        add(panel1, "Center");

        clearDatesAction.setEnabled(false);
        guessDatesAction.setEnabled(false);
        directionCombo.setEnabled(false);
        unitsLabel.setEnabled(false);
        unitsCombo.setEnabled(false);
        scrollPane.setEnabled(false);
        dataTable.setEnabled(false);

        usingSamplingDates.addItemListener(new ItemListener() {
            public void itemStateChanged(ItemEvent ev) {
                boolean enabled = usingSamplingDates.isSelected();
                clearDatesAction.setEnabled(enabled);
                guessDatesAction.setEnabled(enabled);
                unitsLabel.setEnabled(enabled);
                unitsCombo.setEnabled(enabled);
                directionCombo.setEnabled(enabled);
                scrollPane.setEnabled(enabled);
                dataTable.setEnabled(enabled);
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

        for (int i = 0; i < options.taxonList.getTaxonCount(); i++) {
            Date date = options.taxonList.getTaxon(i).getDate();
            double d = date.getTimeValue();

            Date newDate = createDate(d, units, backwards, 0.0);

            options.taxonList.getTaxon(i).setDate(newDate);
        }

        calculateHeights();

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

        unitsCombo.setSelectedIndex(options.datesUnits);
        directionCombo.setSelectedIndex(options.datesDirection);

        calculateHeights();

        dataTableModel.fireTableDataChanged();
    }

    private void setupTable() {
        dataTableModel.fireTableDataChanged();
    }

    public void getOptions(BeautiOptions options) {
        options.datesUnits = unitsCombo.getSelectedIndex();
        options.datesDirection = directionCombo.getSelectedIndex();
    }

    public JComponent getExportableComponent() {
        return dataTable;
    }

    public void selectionChanged() {
        // nothing to do
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

        OptionsPanel optionPanel = new OptionsPanel(12, 12);

        optionPanel.addLabel("The date is given by a numerical field in the taxon label that is:");

        final JRadioButton orderRadio = new JRadioButton("Defined by its order", true);
        final JComboBox orderCombo = new JComboBox(new String[]{"first", "second", "third",
                "fourth", "fourth from last",
                "third from last", "second from last", "last"});

        optionPanel.addComponents(orderRadio, orderCombo);
        optionPanel.addSeparator();

        final JRadioButton prefixRadio = new JRadioButton("Defined by a prefix", false);
        final JTextField prefixText = new JTextField(16);
        prefixText.setEnabled(false);
        optionPanel.addComponents(prefixRadio, prefixText);
        optionPanel.addSeparator();

        final JCheckBox offsetCheck = new JCheckBox("Add the following value to each: ", false);
        final RealNumberField offsetText = new RealNumberField();
        offsetText.setValue(1900);
        offsetText.setColumns(16);
        offsetText.setEnabled(false);
        offsetCheck.addItemListener(new ItemListener() {
            public void itemStateChanged(ItemEvent e) {
                offsetText.setEnabled(offsetCheck.isSelected());
            }
        });
        optionPanel.addComponents(offsetCheck, offsetText);

        final JCheckBox unlessCheck = new JCheckBox("...unless less than:", false);
        final RealNumberField unlessText = new RealNumberField();
        Calendar calendar = GregorianCalendar.getInstance();

        int year = calendar.get(Calendar.YEAR) - 1999;
        unlessText.setValue(year);
        unlessText.setColumns(16);
        unlessText.setEnabled(false);
        optionPanel.addComponents(unlessCheck, unlessText);

        final RealNumberField offset2Text = new RealNumberField();
        offset2Text.setValue(2000);
        offset2Text.setColumns(16);
        offset2Text.setEnabled(false);
        optionPanel.addComponentWithLabel("...in which case add:", offset2Text);

        unlessCheck.addItemListener(new ItemListener() {
            public void itemStateChanged(ItemEvent e) {
                unlessText.setEnabled(unlessCheck.isSelected());
                offset2Text.setEnabled(unlessCheck.isSelected());
            }
        });

        ButtonGroup group = new ButtonGroup();
        group.add(orderRadio);
        group.add(prefixRadio);
        ItemListener listener = new ItemListener() {
            public void itemStateChanged(ItemEvent e) {
                orderCombo.setEnabled(orderRadio.isSelected());
                prefixText.setEnabled(prefixRadio.isSelected());
            }
        };
        orderRadio.addItemListener(listener);
        prefixRadio.addItemListener(listener);

        JOptionPane optionPane = new JOptionPane(optionPanel,
                JOptionPane.QUESTION_MESSAGE,
                JOptionPane.OK_CANCEL_OPTION,
                null,
                null,
                null);
        optionPane.setBorder(new EmptyBorder(12, 12, 12, 12));

        JDialog dialog = optionPane.createDialog(frame, "Guess Dates");
//		dialog.setDefaultCloseOperation(JDialog.DO_NOTHING_ON_CLOSE);
        dialog.setVisible(true);

        if (optionPane.getValue() == null) {
            return;
        }

        int value = (Integer) optionPane.getValue();
        if (value == -1 || value == JOptionPane.CANCEL_OPTION) {
            return;
        }

        DateGuesser guesser = options.dateGuesser;

        guesser.guessDates = true;

        String warningMessage = null;

        for (int i = 0; i < options.taxonList.getTaxonCount(); i++) {
            java.util.Date origin = new java.util.Date(0);

            double d = 0.0;

            try {
                if (orderRadio.isSelected()) {
                    guesser.guessDateFromOrder = true;
                    guesser.order = orderCombo.getSelectedIndex();
                    guesser.fromLast = false;
                    if (guesser.order > 3) {
                        guesser.fromLast = true;
                        guesser.order = 8 - guesser.order - 1;
                    }

                    d = guesser.guessDateFromOrder(options.taxonList.getTaxonId(i), guesser.order, guesser.fromLast);
                } else {
                    guesser.guessDateFromOrder = false;
                    guesser.prefix = prefixText.getText();
                    d = guesser.guessDateFromPrefix(options.taxonList.getTaxonId(i), guesser.prefix);
                }

            } catch (GuessDatesException gfe) {
                warningMessage = gfe.getMessage();
            }

            guesser.offset = 0.0;
            guesser.unlessLessThan = 0.0;
            if (offsetCheck.isSelected()) {
                guesser.offset = offsetText.getValue();
                if (unlessCheck.isSelected()) {
                    guesser.unlessLessThan = unlessText.getValue();
                    guesser.offset2 = offset2Text.getValue();
                    if (d < guesser.unlessLessThan) {
                        d += guesser.offset2;
                    } else {
                        d += guesser.offset;
                    }
                } else {
                    d += guesser.offset;
                }
            }

            Date date = Date.createTimeSinceOrigin(d, Units.Type.YEARS, origin);
            options.taxonList.getTaxon(i).setAttribute("date", date);
        }

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
                Date date = options.taxonList.getTaxon(i).getDate();
                if (date != null) {
                    heights[i] = timeScale.convertTime(date.getTimeValue(), date) - time0;
                    if (heights[i] > options.maximumTipHeight) options.maximumTipHeight = heights[i];
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
            if (col == 0) return true;
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
