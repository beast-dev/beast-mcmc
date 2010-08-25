package dr.app.tracer.application;

import com.lowagie.text.Document;
import com.lowagie.text.DocumentException;
import com.lowagie.text.pdf.DefaultFontMapper;
import com.lowagie.text.pdf.PdfContentByte;
import com.lowagie.text.pdf.PdfTemplate;
import com.lowagie.text.pdf.PdfWriter;
import dr.app.gui.FileDrop;
import dr.app.java16compat.FileNameExtensionFilter;
import dr.app.tracer.analysis.*;
import dr.app.tracer.traces.CombinedTraces;
import dr.app.tracer.traces.FilterDialog;
import dr.app.tracer.traces.TracePanel;
import dr.gui.chart.ChartRuntimeException;
import dr.inference.trace.*;
import org.virion.jam.framework.DocumentFrame;
import org.virion.jam.panels.ActionPanel;
import org.virion.jam.table.TableEditorStopper;
import org.virion.jam.table.TableRenderer;
import org.virion.jam.util.LongTask;

import javax.swing.*;
import javax.swing.border.Border;
import javax.swing.event.ListSelectionEvent;
import javax.swing.event.ListSelectionListener;
import javax.swing.plaf.BorderUIResource;
import javax.swing.table.AbstractTableModel;
import java.awt.*;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.event.WindowAdapter;
import java.awt.event.WindowEvent;
import java.awt.geom.Rectangle2D;
import java.io.*;
import java.text.DecimalFormat;
import java.util.*;
import java.util.List;

public class TracerFrame extends DocumentFrame implements TracerFileMenuHandler, AnalysisMenuHandler {

    private TracePanel tracePanel = null;

    private JTable traceTable = null;
    private TraceTableModel traceTableModel = null;
    private JSplitPane splitPane1 = null;
    private JPanel topPanel = null;

    private JTable statisticTable = null;
    private StatisticTableModel statisticTableModel = null;

    private JScrollPane scrollPane1 = null;

    private JLabel progressLabel;
    private JProgressBar progressBar;

    private final java.util.List<LogFileTraces> traceLists = new ArrayList<LogFileTraces>();
    private final java.util.List<FilteredTraceList> currentTraceLists = new ArrayList<FilteredTraceList>();
    private CombinedTraces combinedTraces = null;

    private final java.util.List<String> commonTraceNames = new ArrayList<String>();
    private boolean homogenousTraceFiles = true;

    private final JComboBox filterCombo = new JComboBox(new String[]{"None"});
    private final JLabel filterStatus = new JLabel();
    String message;

    private FilterDialog filterDialog;

    private int dividerLocation = -1;

    private DemographicDialog demographicDialog = null;
    private BayesianSkylineDialog bayesianSkylineDialog = null;
    private GMRFSkyrideDialog gmrfSkyrideDialog = null;
    private TimeDensityDialog timeDensityDialog = null;
    private LineagesThroughTimeDialog lineagesThroughTimeDialog = null;
    private TraitThroughTimeDialog traitThroughTimeDialog = null;
    private NewTemporalAnalysisDialog createTemporalAnalysisDialog = null;

    private BayesFactorsDialog bayesFactorsDialog = null;

    public TracerFrame(String title) {
        super();

        setTitle(title);

        getOpenAction().setEnabled(false);
        getSaveAction().setEnabled(false);
        getSaveAsAction().setEnabled(false);

        getCutAction().setEnabled(false);
        getCopyAction().setEnabled(false);
        getPasteAction().setEnabled(false);
        getDeleteAction().setEnabled(false);
        getSelectAllAction().setEnabled(false);
        getFindAction().setEnabled(false);

        getZoomWindowAction().setEnabled(false);

        AbstractAction importAction = new AbstractAction("Import Trace File...") {
            public void actionPerformed(ActionEvent ae) {
                doImport();
            }
        };
        setImportAction(importAction);
        setExportAction(exportDataAction);

        setAnalysesEnabled(false);
    }

    public void initializeComponents() {

        setSize(new java.awt.Dimension(1000, 700));

        tracePanel = new TracePanel(this);
        tracePanel.setBorder(new BorderUIResource.EmptyBorderUIResource(new java.awt.Insets(12, 6, 12, 12)));

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
        actionPanel1.setAddAction(getImportAction());
        actionPanel1.setRemoveAction(getRemoveTraceAction());
        getRemoveTraceAction().setEnabled(false);

        JPanel controlPanel1 = new JPanel(new FlowLayout(FlowLayout.LEFT));
        controlPanel1.add(actionPanel1);

        topPanel = new JPanel(new BorderLayout(0, 0));
        topPanel.setBorder(new BorderUIResource.EmptyBorderUIResource(new java.awt.Insets(0, 0, 6, 0)));
        topPanel.add(new JLabel("Trace Files:"), BorderLayout.NORTH);
        topPanel.add(scrollPane1, BorderLayout.CENTER);
        topPanel.add(controlPanel1, BorderLayout.SOUTH);

        statisticTableModel = new StatisticTableModel();
        statisticTable = new JTable(statisticTableModel);
        statisticTable.getColumnModel().getColumn(0).setPreferredWidth(150);
        statisticTable.getColumnModel().getColumn(0).setCellRenderer(renderer);
        statisticTable.getColumnModel().getColumn(1).setPreferredWidth(70);
        statisticTable.getColumnModel().getColumn(1).setCellRenderer(renderer);
        statisticTable.getColumnModel().getColumn(2).setPreferredWidth(70);
        statisticTable.getColumnModel().getColumn(2).setCellRenderer(renderer);
//        ComboBoxRenderer comboBoxRenderer = new ComboBoxRenderer(TraceFactory.TraceType.values());
//        comboBoxRenderer.putClientProperty("JComboBox.isTableCellEditor", Boolean.TRUE);
        statisticTable.getColumnModel().getColumn(3).setPreferredWidth(20);
        statisticTable.getColumnModel().getColumn(3).setCellRenderer(renderer);
        statisticTable.setSelectionMode(ListSelectionModel.MULTIPLE_INTERVAL_SELECTION);

        statisticTable.getSelectionModel().addListSelectionListener(new ListSelectionListener() {
            public void valueChanged(ListSelectionEvent evt) {
                statisticTableSelectionChanged();
            }
        });

//        statisticTable.addMouseListener(new MouseAdapter() { // todo bug
//            private void maybeShowPopup(MouseEvent e) {
//                if (traceTable.getSelectedRowCount() < 2) { // TODO generic to CombinedTraces ?
//                    if (e.isPopupTrigger() && statisticTable.isEnabled()) {
//                        Point p = new Point(e.getX(), e.getY());
//                        int col = statisticTable.columnAtPoint(p);
//                        int row = statisticTable.rowAtPoint(p);
//
//                        if (row >= 0 && row < statisticTable.getRowCount() && col == 3) {
//
//                            // create popup menu...
//                            JPopupMenu contextMenu = createContextMenu();
//
//                            // ... and show it
//                            if (contextMenu != null && contextMenu.getComponentCount() > 0) {
//                                contextMenu.show(statisticTable, p.x, p.y);
//                            }
//                        }
//                        statisticTable.setRowSelectionInterval(row, row);
//                        statisticTableSelectionChanged();
//                    }
//                }
//            }
//
//            public void mousePressed(MouseEvent e) {
//                maybeShowPopup(e);
//            }
//
//            public void mouseReleased(MouseEvent e) {
//                maybeShowPopup(e);
//            }
//        });

        TableEditorStopper.ensureEditingStopWhenTableLosesFocus(statisticTable);

        JScrollPane scrollPane2 = new JScrollPane(statisticTable, JScrollPane.VERTICAL_SCROLLBAR_ALWAYS,
                JScrollPane.HORIZONTAL_SCROLLBAR_NEVER);

        JPanel bottomPanel = new JPanel(new BorderLayout(0, 0));
        bottomPanel.setBorder(new BorderUIResource.EmptyBorderUIResource(new java.awt.Insets(6, 0, 0, 0)));
        bottomPanel.add(new JLabel("Traces:"), BorderLayout.NORTH);
        bottomPanel.add(scrollPane2, BorderLayout.CENTER);
        JPanel changeTraceTypePanel = new JPanel(new FlowLayout(FlowLayout.LEFT));
        changeTraceTypePanel.add(new JLabel("Change To:"));
        JButton realButton = new JButton("real(R)");
        JButton integerButton = new JButton("integer(I)");
        JButton categoryButton = new JButton("category(C)");
        ActionListener traceTypeButton = new ActionListener() {
            public void actionPerformed(ActionEvent e) {
                if (((JButton) e.getSource()).getText().equals("integer(I)")) {
                    changeTraceType(TraceFactory.TraceType.INTEGER);
                } else if (((JButton) e.getSource()).getText().equals("category(C)")) {
                    changeTraceType(TraceFactory.TraceType.CATEGORY);
                } else {
                    changeTraceType(TraceFactory.TraceType.CONTINUOUS);
                }
            }
        };
        realButton.addActionListener(traceTypeButton);
        integerButton.addActionListener(traceTypeButton);
        categoryButton.addActionListener(traceTypeButton);
        changeTraceTypePanel.add(integerButton);
        changeTraceTypePanel.add(categoryButton);
        changeTraceTypePanel.add(realButton);
        changeTraceTypePanel.setToolTipText("<html>Traces Type: real(R) is double, integer(I) is integer, category(C) is string.</html>");
        bottomPanel.add(changeTraceTypePanel, BorderLayout.SOUTH);// todo bug
        
//        bottomPanel.add(new JLabel("<html>Traces Type: real(R) is double, integer(I) is integer, " +
//                "category(C) is string. Right click to change trace type in a selected cell.</html>"),
//                BorderLayout.SOUTH);

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

        JSplitPane splitPane2 = new JSplitPane(JSplitPane.HORIZONTAL_SPLIT, true, leftPanel, tracePanel);
        splitPane2.setBorder(null);
        splitPane2.setDividerLocation(350);

        Color focusColor = UIManager.getColor("Focus.color");
        Border focusBorder = BorderFactory.createMatteBorder(2, 2, 2, 2, focusColor);
        splitPane1.setBorder(BorderFactory.createEmptyBorder(2, 2, 2, 2));
        new FileDrop(null, splitPane1, focusBorder, new FileDrop.Listener() {
            public void filesDropped(java.io.File[] files) {
                importFiles(files);
            }   // end filesDropped
        }); // end FileDrop.Listener

        getContentPane().setLayout(new java.awt.BorderLayout(0, 0));
        getContentPane().add(splitPane2, BorderLayout.CENTER);

        splitPane1.setDividerLocation(2000);

        JPanel filterPanel = new JPanel(new FlowLayout(FlowLayout.LEFT));
        filterPanel.setBorder(new BorderUIResource.EmptyBorderUIResource(new java.awt.Insets(6, 8, 12, 12)));
        JButton filterButton = new JButton("Filtered by :");
        filterButton.addActionListener(new ActionListener() {
            public void actionPerformed(ActionEvent e) {
                int rowIndex = statisticTable.getSelectedRow();
                
                if (filterCombo.getSelectedItem().toString().equalsIgnoreCase("None")) {
                    removeAllFilters();
                } else {
                    invokeFilter();
                }

                statisticTableModel.fireTableDataChanged();
                statisticTable.setRowSelectionInterval(rowIndex, rowIndex);
            }
        });
        filterPanel.add(filterButton);
        filterCombo.setPreferredSize(new Dimension(230, 20));
//        filterCombo.addItemListener(new ItemListener () { //todo bug
//            public void itemStateChanged(ItemEvent e) {
//                if (filterCombo.getItemCount() > 0) {
//                String traceName = filterCombo.getSelectedItem().toString();
//                removeAllFilters(traceName);
//                }
//            }
//        });
        filterPanel.add(filterCombo);
        filterPanel.add(filterStatus);

        getContentPane().add(filterPanel, BorderLayout.SOUTH);
    }

    private void invokeFilter() {
        if (filterDialog == null) {
            filterDialog = new FilterDialog(this);
        }

        if (currentTraceLists == null) {
            JOptionPane.showMessageDialog(this, "There is no file being selected !",
                    "Invalid Action",
                    JOptionPane.ERROR_MESSAGE);
        }

        if (hasDiffValues(currentTraceLists)) {
            JOptionPane.showMessageDialog(this, "Cannot filter multi-files containing different values !",
                    "Invalid Action",
                    JOptionPane.ERROR_MESSAGE);
        }

        String traceName = filterCombo.getSelectedItem().toString();

        message = "  " + filterDialog.showDialog(traceName, currentTraceLists, filterStatus.getText());

        filterStatus.setText(message);

    }

    private void removeAllFilters() {
        int n = JOptionPane.showConfirmDialog(this, "Are you removing all filters of selected files?",
                "Filter Configuration", JOptionPane.YES_NO_OPTION);

        if (n == JOptionPane.YES_OPTION) {
            for (FilteredTraceList filteredTraceList : currentTraceLists) {
                filteredTraceList.removeAllFilters();
            }
        }

        filterStatus.setText("");
    }

    private boolean hasDiffValues(List<FilteredTraceList> currentTraceLists) {
        return false;  //Todo 
    }

    private JPopupMenu createContextMenu() {
        JPopupMenu contextMenu = new JPopupMenu();
        JMenuItem menu = new JMenuItem();
        menu.setText(TraceFactory.TraceType.CONTINUOUS + " (" + TraceFactory.TraceType.CONTINUOUS.getBrief() + ")");
        menu.addActionListener(new ActionListener() {
            public void actionPerformed(ActionEvent e) {
                changeTraceType(TraceFactory.TraceType.CONTINUOUS);
            }
        });
        contextMenu.add(menu);

        menu = new JMenuItem();
        menu.setText(TraceFactory.TraceType.INTEGER + " (" + TraceFactory.TraceType.INTEGER.getBrief() + ")");
        menu.addActionListener(new ActionListener() {
            public void actionPerformed(ActionEvent e) {
                changeTraceType(TraceFactory.TraceType.INTEGER);
            }
        });
        contextMenu.add(menu);

        menu = new JMenuItem();
        menu.setText(TraceFactory.TraceType.CATEGORY + " (" + TraceFactory.TraceType.CATEGORY.getBrief() + ")");
        menu.addActionListener(new ActionListener() {
            public void actionPerformed(ActionEvent e) {
                changeTraceType(TraceFactory.TraceType.CATEGORY);
            }
        });
        contextMenu.add(menu);

        return contextMenu;
    }

    private void changeTraceType(TraceFactory.TraceType newType) {
        int rowIndex = statisticTable.getSelectedRow();
        int selRow = traceTable.getSelectedRow();
        int id = traceLists.get(selRow).getTraceIndex(commonTraceNames.get(rowIndex));
//        Trace trace = traceLists.get(selRow).getTrace(id);

        try {
            traceLists.get(selRow).loadTraces(id, newType);
            traceLists.get(selRow).analyseTrace(id);
        } catch (TraceException e) {
            JOptionPane.showMessageDialog(this, e,
                    "Trace Type Exception",
                    JOptionPane.ERROR_MESSAGE);
        } catch (IOException e) {
            System.err.println("selRow = " + selRow + "; new type = " + newType);
            e.printStackTrace();  //To change body of catch statement use File | Settings | File Templates.
        }
        statisticTableModel.fireTableDataChanged();
        statisticTable.setRowSelectionInterval(rowIndex, rowIndex);
    }

    public void setVisible(boolean b) {
        super.setVisible(b);
        setupDividerLocation();
    }

    private void setupDividerLocation() {

        if (dividerLocation == -1 || dividerLocation == splitPane1.getDividerLocation()) {
            int h0 = topPanel.getHeight();
            int h1 = scrollPane1.getViewport().getHeight();
            int h2 = traceTable.getPreferredSize().height;
            dividerLocation = h2 + h0 - h1;

//		   	int h0 = topPanel.getHeight() - scrollPane1.getViewport().getHeight();
// 			dividerLocation = traceTable.getPreferredSize().height + h0;

            if (dividerLocation > 400) dividerLocation = 400;
            splitPane1.setDividerLocation(dividerLocation);
        }
    }

    public void setAnalysesEnabled(boolean enabled) {
        getDemographicAction().setEnabled(enabled);
        getBayesianSkylineAction().setEnabled(enabled);
        getGMRFSkyrideAction().setEnabled(enabled);
        getLineagesThroughTimeAction().setEnabled(enabled);
        getBayesFactorsAction().setEnabled(enabled);
        getCreateTemporalAnalysisAction().setEnabled(enabled);
        getAddDemographicAction().setEnabled(enabled && temporalAnalysisFrame != null);
        getAddBayesianSkylineAction().setEnabled(enabled && temporalAnalysisFrame != null);
        getAddTimeDensityAction().setEnabled(enabled && temporalAnalysisFrame != null);

        getExportAction().setEnabled(enabled);
        getExportDataAction().setEnabled(enabled);
        getExportPDFAction().setEnabled(enabled);
        getCopyAction().setEnabled(true);
    }

    public void addTraceList(LogFileTraces traceList) {

        int[] selRows = traceTable.getSelectedRows();

        traceLists.add(traceList);

        updateCombinedTraces();

        setAnalysesEnabled(true);

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

        setupDividerLocation();
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

        updateCombinedTraces();

        traceTableModel.fireTableDataChanged();
        statisticTableModel.fireTableDataChanged();

        if (traceLists.size() == 0) {
            getRemoveTraceAction().setEnabled(false);

            setAnalysesEnabled(false);

            currentTraceLists.clear();
            statisticTableModel.fireTableDataChanged();
        }


        if (traceLists.size() > 0) {
            int row = selRows[0];
            if (row >= traceLists.size()) {
                row = traceLists.size() - 1;
            }
            traceTable.getSelectionModel().addSelectionInterval(row, row);
        }
        setupDividerLocation();
    }

    public void setBurnIn(int index, int burnIn) {
        LogFileTraces trace = traceLists.get(index);
        trace.setBurnIn(burnIn);
        analyseTraceList(trace);
        updateCombinedTraces();
        updateTraceTables();
    }

    public void updateCombinedTraces() {
        if (traceLists.size() > 1) {
            FilteredTraceList[] traces = new FilteredTraceList[traceLists.size()];
            try {
                traceLists.toArray(traces);
            } catch (ArrayStoreException ase) {
                combinedTraces = null;
//                JOptionPane.showMessageDialog(this, "",
//                        "Trace Type Exception",
//                        JOptionPane.WARNING_MESSAGE);
            }
            try {
                combinedTraces = new CombinedTraces("Combined", traces);

                analyseTraceList(combinedTraces);
            } catch (TraceException te) {
                // do nothing
            }
        } else {
            combinedTraces = null;
        }
    }

    public void updateTraceTables() {
        int[] selectedTraces = traceTable.getSelectedRows();
        int[] selectedStatistics = statisticTable.getSelectedRows();

        traceTableModel.fireTableDataChanged();
        statisticTableModel.fireTableDataChanged();

        traceTable.getSelectionModel().clearSelection();
        for (int row : selectedTraces) {
            traceTable.getSelectionModel().addSelectionInterval(row, row);
        }

        statisticTable.getSelectionModel().clearSelection();
        for (int row : selectedStatistics) {
            statisticTable.getSelectionModel().addSelectionInterval(row, row);
        }
    }

    public void traceTableSelectionChanged() {
        int[] selRows = traceTable.getSelectedRows();

        if (selRows.length == 0) {
            getRemoveTraceAction().setEnabled(false);
            setAnalysesEnabled(false);
            return;
        }

        setAnalysesEnabled(true);

        getRemoveTraceAction().setEnabled(true);

        currentTraceLists.clear();

        for (int row : selRows) {
            if (row == traceLists.size()) {
                // Combined is include in the selection so disable remove
                getRemoveTraceAction().setEnabled(false);
                currentTraceLists.add(combinedTraces);
            }
        }

        // Get the common set of trace names. This is slightly more complicated
        // that it may seem because we want to keep them in order of the first
        // selected trace file (i.e., as a list). So we populate the list with the
        // first trace file, collect the common set, and then retain only those in
        // the set.
        commonTraceNames.clear();
        homogenousTraceFiles = true;
        Set<String> commonSet = new HashSet<String>();
        boolean isFirst = true;
        for (int row : selRows) {
            if (row < traceLists.size()) {
                FilteredTraceList tl = traceLists.get(row);
                Set<String> nameSet = new HashSet<String>();
                for (int i = 0; i < tl.getTraceCount(); i++) {
                    String traceName = tl.getTraceName(i);
                    nameSet.add(traceName);
                    if (isFirst) {
                        // add them in order of the first trace file
                        commonTraceNames.add(traceName);
                    }
                }

                if (isFirst) {
                    commonSet.addAll(nameSet);
                    isFirst = false;
                } else {
                    if (nameSet.size() != commonSet.size()) {
                        homogenousTraceFiles = false;
                    }
                    commonSet.retainAll(nameSet);
                }

                currentTraceLists.add(tl);
            } else if (isFirst) {
                // if the 'Combined' trace is selected but no other trace files, then add all traces
                FilteredTraceList tl = traceLists.get(0);
                Set<String> nameSet = new HashSet<String>();
                for (int i = 0; i < tl.getTraceCount(); i++) {
                    String traceName = tl.getTraceName(i);
                    nameSet.add(traceName);
                    commonTraceNames.add(traceName);
                }
                commonSet.addAll(nameSet);
            }
        }
        commonTraceNames.retainAll(commonSet);

        int[] rows = statisticTable.getSelectedRows();
        statisticTableModel.fireTableDataChanged();

        if (rows.length > 0) {
            for (int row : rows) {
                statisticTable.getSelectionModel().addSelectionInterval(row, row);
            }
        } else {
            statisticTable.getSelectionModel().setSelectionInterval(0, 0);
        }

        getIntersectionOfSelectedTraceLists();

        message = "  " + updateStatusMessage(currentTraceLists);

        filterStatus.setText(message);
    }

    private String updateStatusMessage(List<FilteredTraceList> currentTraceLists) {
        String message = "";
        List<String> traceNameList = new ArrayList<String>();
        List<String> messageList = new ArrayList<String>();

        for (int i = 0; i < currentTraceLists.size(); i++) {
            Filter f = currentTraceLists.get(i).getFilter();

            if (f != null) {
                String tN = f.getTraceName();
                if (!traceNameList.contains(tN)) {
                    traceNameList.add(tN);
                    message = f.getStatusMessage() + " in file(s) " + "\'" + currentTraceLists.get(i).getName() + "\'";
                    messageList.add(message);
                } else {
                    int id = traceNameList.indexOf(tN);
                    message = messageList.get(id) + " and \'" + currentTraceLists.get(i).getName() + "\'";
                    messageList.set(id, message);
                }

                filterCombo.setSelectedItem(tN);  // todo
            }
        }


        message = "";

        for (String s : messageList) {
            message += s + "; ";
        }

        return message;
    }

    private void getIntersectionOfSelectedTraceLists() {
        filterCombo.removeAllItems();
        filterCombo.addItem("None");

//        Map<String, Class> tracesIntersection = new HashMap<String, Class>(); //names have no order
        List<String> tracesIntersection = Collections.synchronizedList(new ArrayList<String>());
        List<Class> tracesIntersectionClass = Collections.synchronizedList(new ArrayList<Class>());
        List<String> incompatibleTrace = Collections.synchronizedList(new ArrayList<String>());
        for (FilteredTraceList tl : currentTraceLists) {
            List<String> currentTrace = new ArrayList<String>();
            for (int i = 0; i < tl.getTraceCount(); i++) {
                String traceName = tl.getTraceName(i);
                currentTrace.add(traceName);
                if (!incompatibleTrace.contains(traceName)) {
                    Class traceType = tl.getTrace(i).getTraceType();
                    if (traceType == null) {
                        incompatibleTrace.add(traceName);
                        break;
                    }

                    if (tracesIntersection.contains(traceName)) {
                        if (traceType != tracesIntersectionClass.get(tracesIntersection.indexOf(traceName))) {
                            tracesIntersectionClass.remove(tracesIntersection.indexOf(traceName));
                            tracesIntersection.remove(traceName);
                            incompatibleTrace.add(traceName);
                            break;
                        }

                    } else if (currentTraceLists.indexOf(tl) == 0) {

                        tracesIntersection.add(traceName);
                        tracesIntersectionClass.add(traceType);
                    }
                }
            } // end i loop

            for (String traceName : tracesIntersection) {
                if (!currentTrace.contains(traceName)) {
                    tracesIntersectionClass.remove(tracesIntersection.indexOf(traceName));
                    tracesIntersection.remove(traceName);
                    incompatibleTrace.add(traceName);
                }
            }

        }

        assert (tracesIntersection.size() == tracesIntersectionClass.size());

        if (!tracesIntersection.isEmpty()) {
            for (String traceName : tracesIntersection) {
                filterCombo.addItem(traceName);
            }
        }
    }

    public void statisticTableSelectionChanged() {

        int[] selRows = statisticTable.getSelectedRows();

        boolean isIncomplete = false;
        for (FilteredTraceList tl : currentTraceLists) {
            if (tl == null || tl.getTraceCount() == 0 || tl.getStateCount() == 0)
                isIncomplete = true;
        }

        java.util.List<String> selectedTraces = new ArrayList<String>();
        for (int selRow : selRows) {
            selectedTraces.add(commonTraceNames.get(selRow));
        }

        if (currentTraceLists.size() == 0 || isIncomplete) {
            tracePanel.setTraces(null, selectedTraces);
        } else {
            TraceList[] tl = new TraceList[currentTraceLists.size()];
            currentTraceLists.toArray(tl);
            try {
                tracePanel.setTraces(tl, selectedTraces);
            } catch (ChartRuntimeException cre) {
                JOptionPane.showMessageDialog(this, "One or more traces contain invalid values and \rare not able to be displayed.",
                        "Problem with tree file",
                        JOptionPane.ERROR_MESSAGE);
            }
        }
    }

    public void analyseTraceList(FilteredTraceList job) {

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

        private final AnalysisStack<FilteredTraceList> analysisStack = new AnalysisStack<FilteredTraceList>();

        public AnalyseTraceTask() {
        }

        public void add(FilteredTraceList job) {
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
                    FilteredTraceList tl = (FilteredTraceList) job;

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

    public final void doExportData() {

        FileDialog dialog = new FileDialog(this,
                "Export Data...",
                FileDialog.SAVE);

        dialog.setVisible(true);
        if (dialog.getFile() != null) {
            File file = new File(dialog.getDirectory(), dialog.getFile());

            try {
                FileWriter writer = new FileWriter(file);
                writer.write(tracePanel.getExportText());
                writer.close();


            } catch (IOException ioe) {
                JOptionPane.showMessageDialog(this, "Unable to write file: " + ioe,
                        "Unable to write file",
                        JOptionPane.ERROR_MESSAGE);
            }
        }

    }

    public final void doExportPDF() {
        FileDialog dialog = new FileDialog(this,
                "Export PDF Image...",
                FileDialog.SAVE);

        dialog.setVisible(true);
        if (dialog.getFile() != null) {
            File file = new File(dialog.getDirectory(), dialog.getFile());

            Rectangle2D bounds = tracePanel.getExportableComponent().getBounds();
            Document document = new Document(new com.lowagie.text.Rectangle((float) bounds.getWidth(), (float) bounds.getHeight()));
            try {
                // step 2
                PdfWriter writer;
                writer = PdfWriter.getInstance(document, new FileOutputStream(file));
// step 3
                document.open();
// step 4
                PdfContentByte cb = writer.getDirectContent();
                PdfTemplate tp = cb.createTemplate((float) bounds.getWidth(), (float) bounds.getHeight());
                Graphics2D g2d = tp.createGraphics((float) bounds.getWidth(), (float) bounds.getHeight(), new DefaultFontMapper());
                tracePanel.getExportableComponent().print(g2d);
                g2d.dispose();
                cb.addTemplate(tp, 0, 0);
            }
            catch (DocumentException de) {
                JOptionPane.showMessageDialog(this, "Error writing PDF file: " + de,
                        "Export PDF Error",
                        JOptionPane.ERROR_MESSAGE);
            }
            catch (FileNotFoundException e) {
                JOptionPane.showMessageDialog(this, "Error writing PDF file: " + e,
                        "Export PDF Error",
                        JOptionPane.ERROR_MESSAGE);
            }
            document.close();
        }
    }


    public final void doImport() {
        final JFileChooser chooser = new JFileChooser(openDefaultDirectory);
        chooser.setMultiSelectionEnabled(true);

        FileNameExtensionFilter filter = new FileNameExtensionFilter("BEAST log (*.log) Files", "log", "txt");
        chooser.setFileFilter(filter);

        final int returnVal = chooser.showOpenDialog(this);
        if (returnVal == JFileChooser.APPROVE_OPTION) {
            File[] files = chooser.getSelectedFiles();
            importFiles(files);
        }
    }

    private void importFiles(File[] files) {
        LogFileTraces[] traces = new LogFileTraces[files.length];

        for (int i = 0; i < files.length; i++) {
            traces[i] = new LogFileTraces(files[i].getName(), files[i]);
        }

        processTraces(traces);
    }

    private File openDefaultDirectory = null;

    private void setDefaultDir(File file) {
        final String s = file.getAbsolutePath();
        String p = s.substring(0, s.length() - file.getName().length());
        openDefaultDirectory = new File(p);
        if (!openDefaultDirectory.isDirectory()) {
            openDefaultDirectory = null;
        }
    }

    protected void processTraces(final LogFileTraces[] tracesArray) {

        final JFrame frame = this;

        // set default dir to directory of last file
        setDefaultDir(tracesArray[tracesArray.length - 1].getFile());

        if (tracesArray.length == 1) {
            try {
                final LogFileTraces traces = tracesArray[0];


                final String fileName = traces.getName();
                final ProgressMonitorInputStream in = new ProgressMonitorInputStream(
                        this,
                        "Reading " + fileName,
                        new FileInputStream(traces.getFile()));
                in.getProgressMonitor().setMillisToDecideToPopup(0);
                in.getProgressMonitor().setMillisToPopup(0);

                final Reader reader = new InputStreamReader(in);

                Thread readThread = new Thread() {
                    public void run() {
                        try {
                            traces.loadTraces(reader, -1, -1, null);

                            EventQueue.invokeLater(
                                    new Runnable() {
                                        public void run() {
                                            analyseTraceList(traces);
                                            addTraceList(traces);
                                        }
                                    });

                        } catch (final TraceException te) {
                            EventQueue.invokeLater(
                                    new Runnable() {
                                        public void run() {
                                            JOptionPane.showMessageDialog(frame, "Problem with trace file: " + te.getMessage(),
                                                    "Problem with tree file",
                                                    JOptionPane.ERROR_MESSAGE);
                                        }
                                    });
                        } catch (final InterruptedIOException iioex) {
                            // The cancel dialog button was pressed - do nothing
                        } catch (final IOException ioex) {
                            EventQueue.invokeLater(
                                    new Runnable() {
                                        public void run() {
                                            JOptionPane.showMessageDialog(frame, "File I/O Error: " + ioex.getMessage(),
                                                    "File I/O Error",
                                                    JOptionPane.ERROR_MESSAGE);
                                        }
                                    });
//                    } catch (final Exception ex) {
//                        EventQueue.invokeLater (
//                                new Runnable () {
//                                    public void run () {
//                                        JOptionPane.showMessageDialog(frame, "Fatal exception: " + ex.getMessage(),
//                                                "Error reading file",
//                                                JOptionPane.ERROR_MESSAGE);
//                                    }
//                                });
                        }

                    }
                };
                readThread.start();

            } catch (FileNotFoundException fnfe) {
                JOptionPane.showMessageDialog(this, "Unable to open file: File not found",
                        "Unable to open file",
                        JOptionPane.ERROR_MESSAGE);
            } catch (IOException ioex) {
                JOptionPane.showMessageDialog(this, "File I/O Error: " + ioex,
                        "File I/O Error",
                        JOptionPane.ERROR_MESSAGE);
            } catch (Exception ex) {
                JOptionPane.showMessageDialog(this, "Fatal exception: " + ex,
                        "Error reading file",
                        JOptionPane.ERROR_MESSAGE);
            }

        } else {
            Thread readThread = new Thread() {
                public void run() {
                    try {
                        for (final LogFileTraces traces : tracesArray) {
                            final Reader reader = new FileReader(traces.getFile());
                            traces.loadTraces(reader, -1, -1, null);

                            EventQueue.invokeLater(
                                    new Runnable() {
                                        public void run() {
                                            analyseTraceList(traces);
                                            addTraceList(traces);
                                        }
                                    });
                        }

                    } catch (final TraceException te) {
                        EventQueue.invokeLater(
                                new Runnable() {
                                    public void run() {
                                        JOptionPane.showMessageDialog(frame, "Problem with trace file: " + te.getMessage(),
                                                "Problem with tree file",
                                                JOptionPane.ERROR_MESSAGE);
                                    }
                                });
                    } catch (final InterruptedIOException iioex) {
                        // The cancel dialog button was pressed - do nothing
                    } catch (final IOException ioex) {
                        EventQueue.invokeLater(
                                new Runnable() {
                                    public void run() {
                                        JOptionPane.showMessageDialog(frame, "File I/O Error: " + ioex.getMessage(),
                                                "File I/O Error",
                                                JOptionPane.ERROR_MESSAGE);
                                    }
                                });
//                    } catch (final Exception ex) {
//                        EventQueue.invokeLater (
//                                new Runnable () {
//                                    public void run () {
//                                        JOptionPane.showMessageDialog(frame, "Fatal exception: " + ex.getMessage(),
//                                                "Error reading file",
//                                                JOptionPane.ERROR_MESSAGE);
//                                    }
//                                });
                    }

                }
            };
            readThread.start();

        }
    }

    protected boolean readFromFile(File file) throws IOException {
        throw new RuntimeException("Cannot read file - use import instead");
    }

    protected boolean writeToFile(File file) {
        throw new RuntimeException("Cannot write file - this is a read-only application");
    }

    public void doCopy() {
        tracePanel.doCopy();
    }

    private TemporalAnalysisFrame temporalAnalysisFrame = null;

    private void doCreateTemporalAnalysis() {
        if (createTemporalAnalysisDialog == null) {
            createTemporalAnalysisDialog = new NewTemporalAnalysisDialog(this);
        }

        if (createTemporalAnalysisDialog.showDialog() == JOptionPane.CANCEL_OPTION) {
            return;
        }

        temporalAnalysisFrame = createTemporalAnalysisDialog.createTemporalAnalysisFrame(this);

        createTemporalAnalysisAction.setEnabled(false);

        addBayesianSkylineAction.setEnabled(true);
        addDemographicAction.setEnabled(true);
        addTimeDensity.setEnabled(true);

        temporalAnalysisFrame.addWindowListener(new WindowAdapter() {

            public void windowClosing(WindowEvent event) {
                temporalAnalysisFrame = null;
                createTemporalAnalysisAction.setEnabled(true);
                addBayesianSkylineAction.setEnabled(false);
                addDemographicAction.setEnabled(false);
                addTimeDensity.setEnabled(false);
            }
        });
    }

    public void doDemographic(boolean add) {
        if (demographicDialog == null) {
            demographicDialog = new DemographicDialog(this);
        }

        if (currentTraceLists.size() != 1) {
            JOptionPane.showMessageDialog(this, "Please select exactly one trace to do\n" +
                    "this analysis on, or select the Combined trace.",
                    "Unable to perform analysis",
                    JOptionPane.INFORMATION_MESSAGE);
        }


        if (add) {
            if (demographicDialog.showDialog(currentTraceLists.get(0), temporalAnalysisFrame) == JOptionPane.CANCEL_OPTION) {
                return;
            }

            demographicDialog.addToTemporalAnalysis(currentTraceLists.get(0), temporalAnalysisFrame);
        } else {
            if (demographicDialog.showDialog(currentTraceLists.get(0), null) == JOptionPane.CANCEL_OPTION) {
                return;
            }

            demographicDialog.createDemographicFrame(currentTraceLists.get(0), this);
        }
    }

    public void doBayesianSkyline(boolean add) {
        if (bayesianSkylineDialog == null) {
            bayesianSkylineDialog = new BayesianSkylineDialog(this);
        }

        if (currentTraceLists.size() != 1) {
            JOptionPane.showMessageDialog(this, "Please select exactly one trace to do\n" +
                    "this analysis on, (but not the Combined trace).",
                    "Unable to perform analysis",
                    JOptionPane.INFORMATION_MESSAGE);
        }

        if (add) {
            if (bayesianSkylineDialog.showDialog(currentTraceLists.get(0), temporalAnalysisFrame) == JOptionPane.CANCEL_OPTION) {
                return;
            }

            bayesianSkylineDialog.addToTemporalAnalysis(currentTraceLists.get(0), temporalAnalysisFrame);
        } else {
            if (bayesianSkylineDialog.showDialog(currentTraceLists.get(0), null) == JOptionPane.CANCEL_OPTION) {
                return;
            }

            bayesianSkylineDialog.createBayesianSkylineFrame(currentTraceLists.get(0), this);
        }
    }


    public void doGMRFSkyride(boolean add) {
        if (gmrfSkyrideDialog == null) {
            gmrfSkyrideDialog = new GMRFSkyrideDialog(this);
        }

        if (currentTraceLists.size() != 1) {
            JOptionPane.showMessageDialog(this, "Please select exactly one trace to do\n" +
                    "this analysis on, (but not the Combined trace).",
                    "Unable to perform analysis",
                    JOptionPane.INFORMATION_MESSAGE);
        }

        if (add) {
            if (gmrfSkyrideDialog.showDialog(currentTraceLists.get(0), temporalAnalysisFrame) == JOptionPane.CANCEL_OPTION) {
                return;
            }

            gmrfSkyrideDialog.addToTemporalAnalysis(currentTraceLists.get(0), temporalAnalysisFrame);
        } else {
            if (gmrfSkyrideDialog.showDialog(currentTraceLists.get(0), null) == JOptionPane.CANCEL_OPTION) {
                return;
            }

            gmrfSkyrideDialog.createGMRFSkyrideFrame(currentTraceLists.get(0), this);
        }
    }


    public void doLineagesThroughTime(boolean add) {
        if (lineagesThroughTimeDialog == null) {
            lineagesThroughTimeDialog = new LineagesThroughTimeDialog(this);
        }

        if (currentTraceLists.size() != 1) {
            JOptionPane.showMessageDialog(this, "Please select exactly one trace to do\n" +
                    "this analysis on, (but not the Combined trace).",
                    "Unable to perform analysis",
                    JOptionPane.INFORMATION_MESSAGE);
        }

        if (add) {
            if (lineagesThroughTimeDialog.showDialog(currentTraceLists.get(0), temporalAnalysisFrame) == JOptionPane.CANCEL_OPTION) {
                return;
            }

            lineagesThroughTimeDialog.addToTemporalAnalysis(currentTraceLists.get(0), temporalAnalysisFrame);
        } else {
            if (lineagesThroughTimeDialog.showDialog(currentTraceLists.get(0), null) == JOptionPane.CANCEL_OPTION) {
                return;
            }

            lineagesThroughTimeDialog.createLineagesThroughTimeFrame(currentTraceLists.get(0), this);
        }
    }

    public void doTraitThroughTime(boolean add) {
        if (traitThroughTimeDialog == null) {
            traitThroughTimeDialog = new TraitThroughTimeDialog(this);
        }

        if (currentTraceLists.size() != 1) {
            JOptionPane.showMessageDialog(this, "Please select exactly one trace to do\n" +
                    "this analysis on, (but not the Combined trace).",
                    "Unable to perform analysis",
                    JOptionPane.INFORMATION_MESSAGE);
        }

        if (add) {
            if (traitThroughTimeDialog.showDialog(currentTraceLists.get(0), temporalAnalysisFrame) == JOptionPane.CANCEL_OPTION) {
                return;
            }

            traitThroughTimeDialog.addToTemporalAnalysis(currentTraceLists.get(0), temporalAnalysisFrame);
        } else {
            if (traitThroughTimeDialog.showDialog(currentTraceLists.get(0), null) == JOptionPane.CANCEL_OPTION) {
                return;
            }

            traitThroughTimeDialog.createTraitThroughTimeFrame(currentTraceLists.get(0), this);
        }
    }

    private void doAddTimeDensity() {
        if (timeDensityDialog == null) {
            timeDensityDialog = new TimeDensityDialog(this);
        }

        if (currentTraceLists.size() != 1) {
            JOptionPane.showMessageDialog(this, "Please select exactly one trace to do\n" +
                    "this analysis on, (or the Combined trace).",
                    "Unable to perform analysis",
                    JOptionPane.INFORMATION_MESSAGE);
        }

        if (timeDensityDialog.showDialog(currentTraceLists.get(0), temporalAnalysisFrame) == JOptionPane.CANCEL_OPTION) {
            return;
        }

        timeDensityDialog.addToTemporalAnalysis(currentTraceLists.get(0), temporalAnalysisFrame);
    }

    private void doCalculateBayesFactors() {
        if (bayesFactorsDialog == null) {
            bayesFactorsDialog = new BayesFactorsDialog(this);
        }

        if (bayesFactorsDialog.showDialog(currentTraceLists) == JOptionPane.CANCEL_OPTION) {
            return;
        }

        bayesFactorsDialog.createBayesFactorsFrame(currentTraceLists, this);

    }

    public JComponent getExportableComponent() {

        return tracePanel.getExportableComponent();
    }

    class TraceTableModel extends AbstractTableModel {
        final String[] columnNames = {"Tree File", "States", "Burn-In"};

        public int getColumnCount() {
            return columnNames.length;
        }

        public int getRowCount() {
            int n = traceLists.size();
            if (n == 0 || combinedTraces != null) n++;
            return n;
        }

        public String getColumnName(int col) {
            return columnNames[col];
        }

        public Object getValueAt(int row, int col) {
            FilteredTraceList traceList;

            if (traceLists.size() == 0) {
                switch (col) {
                    case 0:
                        return "No files loaded";
                    case 1:
                        return "";
                    case 2:
                        return "";
                }
            } else if (row == traceLists.size()) {
                traceList = combinedTraces;
                switch (col) {
                    case 0:
                        return traceList.getName();
                    case 1:
                        return traceList.getMaxState();
                    case 2:
                        return "-";
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

    class StatisticTableModel extends AbstractTableModel {
        final String[] columnNames = {"Statistic", "Mean", "ESS", "Type"};

        private final DecimalFormat formatter = new DecimalFormat("0.###E0");
        private final DecimalFormat formatter2 = new DecimalFormat("####0.###");

        public int getColumnCount() {
            return columnNames.length;
        }

        public int getRowCount() {
            if (currentTraceLists.size() == 0 || currentTraceLists.get(0) == null) return 0;
            return commonTraceNames.size();
        }

        public String getColumnName(int col) {
            return columnNames[col];
        }

        public Object getValueAt(int row, int col) {
            String traceName = commonTraceNames.get(row);

            if (col == 0) return traceName;

            if (!homogenousTraceFiles) {
                return "n/a";
            }

            TraceDistribution td = currentTraceLists.get(0).getDistributionStatistics(row);
            if (td == null) return "-";
            if (col == 3) return td.getTraceType().getBrief();

            double value = 0.0;
            boolean warning = false;
            boolean extremeWarning = false;
            switch (col) {
                case 1:
                    value = td.getMean();
                    break;
                case 2:
                    if (!td.isValid()) return "-";
                    value = td.getESS();
                    if( Double.isNaN(value) || value < 1 ) {
                        // assume not applicable; should be tested in the computation
                       return "-";
                    }
                    if (value < 200.0) warning = true;
                    if (value < 100.0) extremeWarning = true;
                    value = Math.round(value);
                    break;
            }

            String string;
            if (Math.abs(value) < 0.1 || Math.abs(value) >= 100000.0) {
                string = formatter.format(value);
            } else string = formatter2.format(value);

            if (warning) {
                return "<html><font color=\"" + (extremeWarning ? "#EE0000" : "#EEAA00") + "\">" + string + "</font></html> ";
            }

            return string;
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

    public Action getExportDataAction() {
        return exportDataAction;
    }

    public Action getExportPDFAction() {
        return exportPDFAction;
    }


    public Action getRemoveTraceAction() {
        return removeTraceAction;
    }

    public Action getDemographicAction() {
        return demographicAction;
    }

    public Action getBayesianSkylineAction() {
        return bayesianSkylineAction;
    }

    public Action getGMRFSkyrideAction() {
        return gmrfSkyrideAction;
    }

    public Action getLineagesThroughTimeAction() {
        return lineagesThroughTimeAction;
    }

    public Action getTraitThroughTimeAction() {
        return traitThroughTimeAction;
    }

    public Action getCreateTemporalAnalysisAction() {
        return createTemporalAnalysisAction;
    }

    public Action getAddDemographicAction() {
        return addDemographicAction;
    }

    public Action getAddBayesianSkylineAction() {
        return addBayesianSkylineAction;
    }

    public Action getAddTimeDensityAction() {
        return addTimeDensity;
    }

    public Action getBayesFactorsAction() {
        return bayesFactorsAction;
    }

    private final AbstractAction demographicAction = new AbstractAction(AnalysisMenuFactory.DEMOGRAPHIC_RECONSTRUCTION) {
        public void actionPerformed(ActionEvent ae) {
            doDemographic(false);
        }
    };

    private final AbstractAction bayesianSkylineAction = new AbstractAction(AnalysisMenuFactory.BAYESIAN_SKYLINE_RECONSTRUCTION) {
        public void actionPerformed(ActionEvent ae) {
            doBayesianSkyline(false);
        }
    };

    private final AbstractAction gmrfSkyrideAction = new AbstractAction(AnalysisMenuFactory.GMRF_SKYRIDE_RECONSTRUCTION) {
        public void actionPerformed(ActionEvent ae) {
            doGMRFSkyride(false);
        }
    };

    private final AbstractAction lineagesThroughTimeAction = new AbstractAction(AnalysisMenuFactory.LINEAGES_THROUGH_TIME) {
        public void actionPerformed(ActionEvent ae) {
            doLineagesThroughTime(false);
        }
    };

    private final AbstractAction traitThroughTimeAction = new AbstractAction(AnalysisMenuFactory.TRAIT_THROUGH_TIME) {
        public void actionPerformed(ActionEvent ae) {
            doTraitThroughTime(false);
        }
    };

    private final AbstractAction createTemporalAnalysisAction = new AbstractAction(AnalysisMenuFactory.CREATE_TEMPORAL_ANALYSIS) {
        public void actionPerformed(ActionEvent ae) {
            doCreateTemporalAnalysis();
        }
    };

    private final AbstractAction addDemographicAction = new AbstractAction(AnalysisMenuFactory.ADD_DEMOGRAPHIC_RECONSTRUCTION) {
        public void actionPerformed(ActionEvent ae) {
            doDemographic(true);
        }
    };

    private final AbstractAction addBayesianSkylineAction = new AbstractAction(AnalysisMenuFactory.ADD_BAYESIAN_SKYLINE_RECONSTRUCTION) {
        public void actionPerformed(ActionEvent ae) {
            doBayesianSkyline(true);
        }
    };

    private final AbstractAction addTimeDensity = new AbstractAction(AnalysisMenuFactory.ADD_TIME_DENSITY) {
        public void actionPerformed(ActionEvent ae) {
            doAddTimeDensity();
        }
    };

    private final AbstractAction bayesFactorsAction = new AbstractAction(AnalysisMenuFactory.CALCULATE_BAYES_FACTORS) {
        public void actionPerformed(ActionEvent ae) {
            doCalculateBayesFactors();
        }
    };

    private final AbstractAction removeTraceAction = new AbstractAction() {
        public void actionPerformed(ActionEvent ae) {
            removeTraceList();
        }
    };

    private final AbstractAction exportDataAction = new AbstractAction("Export Data...") {
        public void actionPerformed(ActionEvent ae) {
            doExportData();
        }
    };

    private final AbstractAction exportPDFAction = new AbstractAction("Export PDF...") {
        public void actionPerformed(ActionEvent ae) {
            doExportPDF();
        }
    };

}