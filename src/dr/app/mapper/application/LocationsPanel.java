/*
 * LocationsPanel.java
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

package dr.app.mapper.application;

import dr.app.gui.FileDrop;
import dr.app.gui.table.DateCellEditor;
import dr.app.gui.table.TableEditorStopper;
import dr.app.gui.table.TableSorter;
import dr.app.gui.util.LongTask;
import dr.app.mapper.application.mapper.Layer;
import dr.inference.trace.LogFileTraces;
import dr.inference.trace.TraceException;
import dr.inference.trace.TraceList;
import jam.framework.Exportable;
import jam.panels.ActionPanel;
import jam.table.HeaderRenderer;
import jam.table.TableRenderer;

import javax.swing.*;
import javax.swing.border.Border;
import javax.swing.event.ListSelectionEvent;
import javax.swing.event.ListSelectionListener;
import javax.swing.plaf.BorderUIResource;
import javax.swing.table.AbstractTableModel;
import javax.swing.table.JTableHeader;
import java.awt.*;
import java.awt.List;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.event.MouseEvent;
import java.io.*;
import java.util.*;

/**
 * @author Andrew Rambaut
 * @version $Id: StrainsPanel.java,v 1.17 2006/09/05 13:29:34 rambaut Exp $
 */
public class LocationsPanel extends JPanel implements Exportable, MapperDocument.Listener {

    private JTable traceTable = null;
    private TraceTableModel traceTableModel = null;
    private JSplitPane splitPane1 = null;
    private JPanel topPanel = null;

    private JTable layerTable = null;
    private LayerTableModel layerTableModel = null;

    private JTable dataTable = null;
    private DataTableModel dataTableModel = null;

    private final MapperFrame frame;
    private final MapperDocument document;


    private JScrollPane scrollPane1;

    private JLabel progressLabel;
    private JProgressBar progressBar;

    private final java.util.List<LogFileTraces> traceLists = new ArrayList<LogFileTraces>();
    private final java.util.List<Layer> layers = new ArrayList<Layer>();

    String message = "";
    private int dividerLocation = -1;

    public LocationsPanel(final MapperFrame parent, final MapperDocument document) {

        this.frame = parent;
        this.document = document;

        traceTableModel = new TraceTableModel();
        traceTable = new JTable(traceTableModel);
        TableRenderer renderer = new TableRenderer(SwingConstants.LEFT, new Insets(0, 4, 0, 4));
        traceTable.getColumnModel().getColumn(0).setCellRenderer(renderer);
        traceTable.getColumnModel().getColumn(1).setPreferredWidth(50);
        traceTable.getColumnModel().getColumn(1).setCellRenderer(renderer);
        traceTable.getColumnModel().getColumn(2).setPreferredWidth(50);
        traceTable.getColumnModel().getColumn(2).setCellRenderer(renderer);
        traceTable.setSelectionMode(ListSelectionModel.MULTIPLE_INTERVAL_SELECTION);

        traceTable.getSelectionModel().addListSelectionListener(new ListSelectionListener() {
            public void valueChanged(ListSelectionEvent evt) {
                traceTableSelectionChanged();
            }
        });

        scrollPane1 = new JScrollPane(traceTable, JScrollPane.VERTICAL_SCROLLBAR_ALWAYS,
                JScrollPane.HORIZONTAL_SCROLLBAR_NEVER);

        ActionPanel actionPanel1 = new ActionPanel(false);
        actionPanel1.setAddAction(frame.getImportLocationsAction());
        actionPanel1.setRemoveAction(frame.getDeleteItemAction());
        frame.getDeleteItemAction().setEnabled(false);

        JPanel controlPanel1 = new JPanel(new FlowLayout(FlowLayout.LEFT));
        controlPanel1.add(actionPanel1);

        topPanel = new JPanel(new BorderLayout(0, 0));
        topPanel.setBorder(new BorderUIResource.EmptyBorderUIResource(new java.awt.Insets(0, 0, 6, 0)));
        topPanel.add(new JLabel("Trace Files:"), BorderLayout.NORTH);
        topPanel.add(scrollPane1, BorderLayout.CENTER);
        topPanel.add(controlPanel1, BorderLayout.SOUTH);

        layerTableModel = new LayerTableModel();
        layerTable = new JTable(layerTableModel) {
            //Implement table header tool tips.
            protected JTableHeader createDefaultTableHeader() {
                return new JTableHeader(columnModel) {
//                    public String getToolTipText(MouseEvent e) {
//                        java.awt.Point p = e.getPoint();
//                        int index = columnModel.getColumnIndexAtX(p.x);
//                        int realIndex = columnModel.getColumn(index).getModelIndex();
//                        return columnToolTips[realIndex];
//                    }
                };
            }
        };
        layerTable.getColumnModel().getColumn(0).setPreferredWidth(150);
        layerTable.getColumnModel().getColumn(0).setCellRenderer(renderer);
        layerTable.getColumnModel().getColumn(1).setPreferredWidth(70);
        layerTable.getColumnModel().getColumn(1).setCellRenderer(renderer);
//        layerTable.getColumnModel().getColumn(2).setPreferredWidth(70);
//        layerTable.getColumnModel().getColumn(2).setCellRenderer(renderer);
//        ComboBoxRenderer comboBoxRenderer = new ComboBoxRenderer(TraceFactory.TraceType.values());
//        comboBoxRenderer.putClientProperty("JComboBox.isTableCellEditor", Boolean.TRUE);
//        layerTable.getColumnModel().getColumn(3).setPreferredWidth(20);
//        layerTable.getColumnModel().getColumn(3).setCellRenderer(renderer);
        layerTable.setSelectionMode(ListSelectionModel.MULTIPLE_INTERVAL_SELECTION);

        layerTable.getSelectionModel().addListSelectionListener(new ListSelectionListener() {
            public void valueChanged(ListSelectionEvent evt) {
                layerTableSelectionChanged();
            }
        });

        TableEditorStopper.ensureEditingStopWhenTableLosesFocus(layerTable);

        JScrollPane scrollPane2 = new JScrollPane(layerTable, JScrollPane.VERTICAL_SCROLLBAR_ALWAYS,
                JScrollPane.HORIZONTAL_SCROLLBAR_NEVER);

        JPanel bottomPanel = new JPanel(new BorderLayout(0, 0));
        bottomPanel.add(new JLabel("Layers:"), BorderLayout.NORTH);
        bottomPanel.add(scrollPane2, BorderLayout.CENTER);

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

        JScrollPane scrollPane = new JScrollPane(dataTable,
                JScrollPane.VERTICAL_SCROLLBAR_ALWAYS,
                JScrollPane.HORIZONTAL_SCROLLBAR_ALWAYS);
        scrollPane.setOpaque(false);

        JToolBar toolBar1 = new JToolBar();
        toolBar1.setFloatable(false);
        toolBar1.setOpaque(false);

//        toolBar1.setLayout(new FlowLayout(FlowLayout.LEFT, 0, 0));
//        JButton button = new JButton(clearDatesAction);
//        PanelUtils.setupComponent(button);
//        toolBar1.add(button);
//        button = new JButton(guessDatesAction);
//        PanelUtils.setupComponent(button);
//        toolBar1.add(button);
//        toolBar1.add(new JToolBar.Separator(new Dimension(12, 12)));
//        final JLabel unitsLabel = new JLabel("Dates specified as ");
//        toolBar1.add(unitsLabel);
//        toolBar1.add(unitsCombo);
//        toolBar1.add(directionCombo);

        setOpaque(false);
        setBorder(new BorderUIResource.EmptyBorderUIResource(new Insets(12, 12, 12, 12)));
        setLayout(new BorderLayout(0, 0));

        JPanel leftPanel = new JPanel(new BorderLayout(0, 0));
        leftPanel.setPreferredSize(new Dimension(400, 300));
        splitPane1 = new JSplitPane(JSplitPane.VERTICAL_SPLIT, true, topPanel, bottomPanel);
        splitPane1.setBorder(null);

        JPanel progressPanel = new JPanel(new BorderLayout(0, 0));
        progressLabel = new JLabel("");
        progressBar = new JProgressBar();
        progressPanel.add(progressLabel, BorderLayout.NORTH);
        progressPanel.add(progressBar, BorderLayout.CENTER);
        progressPanel.setBorder(new BorderUIResource.EmptyBorderUIResource(new java.awt.Insets(6, 0, 0, 0)));

        leftPanel.add(splitPane1, BorderLayout.CENTER);
        leftPanel.add(progressPanel, BorderLayout.SOUTH);
        leftPanel.setBorder(new BorderUIResource.EmptyBorderUIResource(new java.awt.Insets(12, 12, 12, 6)));

        JSplitPane splitPane2 = new JSplitPane(JSplitPane.HORIZONTAL_SPLIT, true, leftPanel, scrollPane);
        splitPane2.setBorder(null);
        splitPane2.setDividerLocation(350);

        Color focusColor = UIManager.getColor("Focus.color");
        Border focusBorder = BorderFactory.createMatteBorder(2, 2, 2, 2, focusColor);
        splitPane1.setBorder(BorderFactory.createEmptyBorder(2, 2, 2, 2));
        new FileDrop(null, splitPane1, focusBorder, new FileDrop.Listener() {
            public void filesDropped(java.io.File[] files) {
                frame.importLocationFiles(files);
            }   // end filesDropped
        }); // end FileDrop.Listener

        add(toolBar1, "North");
        add(splitPane2, "Center");

    }

    @Override
    public void taxaChanged() {
        dataTableModel.fireTableDataChanged();
    }

    public void addTraceList(LogFileTraces traceList) {
        int[] selRows = traceTable.getSelectedRows();

        traceLists.add(traceList);

//        setAnalysesEnabled(true);

        traceTableModel.fireTableDataChanged();

        int newRow = traceLists.size() - 1;
        traceTable.getSelectionModel().setSelectionInterval(newRow, newRow);
        if (selRows.length > 1) {
            for (int row : selRows) {
                if (row == traceLists.size() - 1) {
                    row = traceLists.size();
                }
                traceTable.getSelectionModel().addSelectionInterval(row, row);
            }
        }

//        setupDividerLocation();
    }

    private void removeTraceList() {
        int[] selRows = traceTable.getSelectedRows();

        LogFileTraces[] tls = new LogFileTraces[selRows.length];
        int i = 0;
        for (int row : selRows) {
            tls[i] = traceLists.get(row);
            i++;
        }
        for (LogFileTraces tl : tls) {
            traceLists.remove(tl);
        }

        traceTableModel.fireTableDataChanged();
//        layerTableModel.fireTableDataChanged();

        if (traceLists.size() == 0) {
            frame.getDeleteItemAction().setEnabled(false);

//            setAnalysesEnabled(false);

//            layerTableModel.fireTableDataChanged();
        }


        if (traceLists.size() > 0) {
            int row = selRows[0];
            if (row >= traceLists.size()) {
                row = traceLists.size() - 1;
            }
            traceTable.getSelectionModel().addSelectionInterval(row, row);
        }
//        setupDividerLocation();
    }

    public void setBurnIn(int index, int burnIn) {
        LogFileTraces trace = traceLists.get(index);
        trace.setBurnIn(burnIn);
        analyseTraceList(trace);
        updateTraceTables();
    }

    public void updateTraceTables() {
        int[] selectedTraces = traceTable.getSelectedRows();
//        int[] selectedStatistics = layerTable.getSelectedRows();

        traceTableModel.fireTableDataChanged();
//        layerTableModel.fireTableDataChanged();

        traceTable.getSelectionModel().clearSelection();
        for (int row : selectedTraces) {
            traceTable.getSelectionModel().addSelectionInterval(row, row);
        }

//        layerTable.getSelectionModel().clearSelection();
//        for (int row : selectedStatistics) {
//            layerTable.getSelectionModel().addSelectionInterval(row, row);
//        }
    }

    public void traceTableSelectionChanged() {
        int[] selRows = traceTable.getSelectedRows();

        if (selRows.length == 0) {
            frame.getDeleteItemAction().setEnabled(false);
//            setAnalysesEnabled(false);
            return;
        }

//        setAnalysesEnabled(true);

        frame.getDeleteItemAction().setEnabled(true);

//        int[] rows = layerTable.getSelectedRows();
//        layerTableModel.fireTableDataChanged();
//
//        if (rows.length > 0) {
//            for (int row : rows) {
//                layerTable.getSelectionModel().addSelectionInterval(row, row);
//            }
//        } else {
//            layerTable.getSelectionModel().setSelectionInterval(0, 0);
//        }
    }


    public void layerTableSelectionChanged() {

    }

    public void analyseTraceList(TraceList job) {

        if (analyseTask == null) {
            analyseTask = new AnalyseTraceTask();

            javax.swing.Timer timer = new javax.swing.Timer(1000, new ActionListener() {
                public void actionPerformed(ActionEvent evt) {
                    progressBar.setMaximum(analyseTask.getLengthOfTask());
                    progressBar.setValue(analyseTask.getCurrent());
                }
            });

            analyseTask.go();
            timer.start();
        }

        analyseTask.add(job);
    }

    AnalyseTraceTask analyseTask = null;

    class AnalyseTraceTask extends LongTask {

        class AnalysisStack<T> {
            private final java.util.List<T> jobs = new ArrayList<T>();

            public synchronized void add(T job) {
                jobs.add(job);
            }

            public synchronized int getCount() {
                return jobs.size();
            }

            public synchronized T get(int index) {
                return jobs.get(index);
            }

            public synchronized void remove(int index) {
                jobs.remove(index);
            }
        }

        private final AnalysisStack<TraceList> analysisStack = new AnalysisStack<TraceList>();

        public AnalyseTraceTask() {
        }

        public void add(TraceList job) {
            analysisStack.add(job);
            current = 0;
        }

        public int getCurrent() {
            return current;
        }

        public int getLengthOfTask() {
            int count = 0;
            for (int i = 0; i < analysisStack.getCount(); i++) {
                count += analysisStack.get(i).getTraceCount();
            }
            return count;
        }

        public void stop() {
        }

        public boolean done() {
            return false;
        }

        public String getDescription() {
            return "Analysing Trace File...";
        }

        public String getMessage() {
            return null;
        }

        public Object doWork() {

            current = 0;
            boolean textCleared = true;

            do {
                if (analysisStack.getCount() > 0) {
                    Object job = analysisStack.get(0);
                    TraceList tl = (TraceList) job;

                    try {
                        for (int i = 0; i < tl.getTraceCount(); i++) {
                            progressLabel.setText("Analysing " + tl.getName() + ":");
                            textCleared = false;
                            tl.analyseTrace(i);
                            repaint();
                            current += 1;
                        }
                    } catch (final Exception ex) {
                        // do nothing. An exception is sometimes fired when burnin is changed whilst in the
                        // middle of an analysis. This doesn't seem to matter as the analysis is restarted.

                        ex.printStackTrace();
//                        EventQueue.invokeLater (
//								new Runnable () {
//									public void run () {
//										JOptionPane.showMessageDialog(TracerFrame.this, "Fatal exception: " + ex.getMessage(),
//												"Error reading file",
//												JOptionPane.ERROR_MESSAGE);
//									}
//								});
                    }
                    analysisStack.remove(0);
                } else {
                    if (!textCleared) {
                        progressLabel.setText("");
                        textCleared = true;
                    }
                    try {
                        Thread.sleep(500);
                    } catch (InterruptedException ie) {
                        // do nothing
                    }
                }
            } while (true);
        }

        //private int lengthOfTask = 0;
        private int current = 0;
        //private String message;
    }


    public JComponent getExportableComponent() {
        return dataTable;
    }

    public void selectionChanged() {
        // nothing to do
    }

    class TraceTableModel extends AbstractTableModel {
        final String[] columnNames = {"Trace File", "States", "Burn-In"};

        public int getColumnCount() {
            return columnNames.length;
        }

        public int getRowCount() {
            int n = traceLists.size();
            if (n == 0) n++;
            return n;
        }

        public String getColumnName(int col) {
            return columnNames[col];
        }

        public Object getValueAt(int row, int col) {
            TraceList traceList;

            if (traceLists.size() == 0) {
                switch (col) {
                    case 0:
                        return "No files loaded";
                    case 1:
                        return "";
                    case 2:
                        return "";
                }
            } else {
                traceList = traceLists.get(row);
                switch (col) {
                    case 0:
                        return traceList.getName();
                    case 1:
                        return traceList.getMaxState();
                    case 2:
                        return traceList.getBurnIn();
                }
            }

            return null;
        }

        public void setValueAt(Object value, int row, int col) {
            if (col == 2) {
                setBurnIn(row, (Integer) value);
            }
        }

        public Class getColumnClass(int c) {
            return getValueAt(0, c).getClass();
        }

        public boolean isCellEditable(int row, int col) {
            return col == 2 && row < traceLists.size();
        }
    }

    class LayerTableModel extends AbstractTableModel {
        final String[] columnNames = {"Layer", "Type"};

        public int getColumnCount() {
            return columnNames.length;
        }

        public int getRowCount() {
            return layers.size();
        }

        public String getColumnName(int col) {
            return columnNames[col];
        }

        public Object getValueAt(int row, int col) {
            Layer layer = layers.get(row);

            if (col == 0) return layer.getName();
            if (col == 1) return layer.getType().toString();

            return "";
        }

        public boolean isCellEditable(int row, int col) {
            return true;

        }

        public Class getColumnClass(int c) {
            if (getRowCount() == 0) {
                return Object.class;
            }
            return getValueAt(0, c).getClass();
        }
    }

    class DataTableModel extends AbstractTableModel {

        String[] columnNames = {"Serum", "Virus", "Titre", "Table"};

        public DataTableModel() {
        }

        public int getColumnCount() {
            return columnNames.length;
        }

        public int getRowCount() {
//            java.util.List<Taxon> taxonList = document.getTaxa();
//
//            if (taxonList == null) return 0;
//
//            return taxonList.size();
            return 0;
        }

        public Object getValueAt(int row, int col) {
//            java.util.List<Taxon> taxonList = document.getTaxa();
//
//            switch (col) {
//                case 0:
//                    return taxonList.get(row);
//                case 1:
//                    Date date = taxonList.get(row).getDate();
//                    if (date != null) {
//                        return date.getTimeValue();
//                    } else {
//                        return "-";
//                    }
//                case 2:
//                    if (heights != null) {
//                        return heights[row];
//                    } else {
//                        return "0.0";
//                    }
//            }
            return null;
        }

        public void setValueAt(Object aValue, int row, int col) {
//            java.util.List<Taxon> taxonList = document.getTaxa();
//
//            if (col == 0) {
//                taxonList.get(row).setId(aValue.toString());
//            } else if (col == 1) {
//                Date date = taxonList.get(row).getDate();
//                if (date != null) {
//                    double d = (Double) aValue;
//                    Date newDate = createDate(d, date.getUnits(), date.isBackwards(), date.getOrigin());
//                    taxonList.get(row).setDate(newDate);
//                }
//            }
//
//            timeScaleChanged();
        }

        public boolean isCellEditable(int row, int col) {
//            if (col == 0) return true;
//            if (col == 1) {
//                Date date = document.getTaxa().get(row).getDate();
//                return (date != null);
//            }
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
