/*
 * MapperFrame.java
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

/*
 * BeautiFrame.java
 *
 * (c) 2002-2005 BEAST Development Core Team
 *
 * This package may be distributed under the
 * Lesser Gnu Public Licence (LGPL)
 */
package dr.app.mapper.application;

import dr.app.mapper.application.menus.MapperFileMenuHandler;
import dr.evolution.tree.FlexibleTree;
import dr.evolution.tree.Tree;
import dr.inference.trace.LogFileTraces;
import dr.inference.trace.TraceException;
import dr.util.DataTable;
import jam.framework.DocumentFrame;
import jam.framework.Exportable;

import javax.swing.*;
import javax.swing.filechooser.FileNameExtensionFilter;
import javax.swing.plaf.BorderUIResource;
import java.awt.*;
import java.awt.datatransfer.Clipboard;
import java.awt.datatransfer.StringSelection;
import java.awt.event.ActionEvent;
import java.io.*;
import java.util.ArrayList;
import java.util.List;

/**
 * @author Andrew Rambaut
 * @version $Id: BeautiFrame.java,v 1.22 2006/09/09 16:07:06 rambaut Exp $
 */
public class MapperFrame extends DocumentFrame implements MapperFileMenuHandler {

    private final JTabbedPane tabbedPane = new JTabbedPane();
    private final JLabel statusLabel = new JLabel("No data loaded");

    private StrainsPanel strainsPanel;
    private MeasurementsPanel measurementsPanel;
    private LocationsPanel locationsPanel;
//    private AnalysisPanel analysisPanel;

    MapperDocument document = new MapperDocument();

    List<Tree> trees = new ArrayList<Tree>();

    public MapperFrame(String title) {
        super();

        setTitle(title);

        getExportDataAction().setEnabled(false);
        getExportPDFAction().setEnabled(false);

        // Prevent the application to close in requestClose()
        // after a user cancel or a failure in beast file generation
        setDefaultCloseOperation(DO_NOTHING_ON_CLOSE);

        getOpenAction().setEnabled(true);
        getSaveAction().setEnabled(false);
        getSaveAsAction().setEnabled(false);

        getFindAction().setEnabled(false);

        getCutAction().setEnabled(false);
        getPasteAction().setEnabled(false);
        getDeleteAction().setEnabled(false);
        getSelectAllAction().setEnabled(false);

        getCopyAction().setEnabled(false);

        getZoomWindowAction().setEnabled(false);
    }

    public void initializeComponents() {

        strainsPanel = new StrainsPanel(this, document);
        measurementsPanel = new MeasurementsPanel(this, document);
        locationsPanel = new LocationsPanel(this, document);
//        analysisPanel = new AnalysisPanel(this, trees.get(0));

        tabbedPane.addTab("Measurements", measurementsPanel);
        tabbedPane.addTab("Strains", strainsPanel);
        tabbedPane.addTab("Locations", locationsPanel);
//        tabbedPane.addTab("Analysis", analysisPanel);

        JPanel panel = new JPanel(new BorderLayout(6, 6));
        panel.setBorder(new BorderUIResource.EmptyBorderUIResource(new Insets(12, 12, 12, 12)));
        panel.add(tabbedPane, BorderLayout.CENTER);

        panel.add(statusLabel, BorderLayout.SOUTH);

        getContentPane().setLayout(new BorderLayout(0, 0));
        getContentPane().add(panel, BorderLayout.CENTER);

        setSize(new Dimension(1024, 768));

        setStatusMessage();
    }

    public void timeScaleChanged() {
//        analysisPanel.timeScaleChanged();
        setStatusMessage();
    }

    public void measurementsChanged() {
        //To change body of created methods use File | Settings | File Templates.
    }


    protected boolean readFromFile(File file) throws IOException {
        Reader reader = new FileReader(file);

        BufferedReader bufferedReader = new BufferedReader(reader);
        String line = bufferedReader.readLine();
        while (line != null && line.length() == 0) {
            line = bufferedReader.readLine();
        }

        boolean isNexus = (line != null && line.toUpperCase().contains("#NEXUS"));

        reader = new FileReader(file);

        Tree tree = null;
        try {
//            if (isNexus) {
//                NexusImporter importer = new NexusImporter(reader);
//                tree = importer.importTree(taxa);
//            } else {
//                NewickImporter importer = new NewickImporter(reader);
//                tree = importer.importTree(taxa);
//            }

//        } catch (Importer.ImportException ime) {
//            JOptionPane.showMessageDialog(this, "Error parsing imported file: " + ime,
//                    "Error reading file",
//                    JOptionPane.ERROR_MESSAGE);
//            ime.printStackTrace();
//            return false;
//        } catch (IOException ioex) {
//            JOptionPane.showMessageDialog(this, "File I/O Error: " + ioex,
//                    "File I/O Error",
//                    JOptionPane.ERROR_MESSAGE);
//            ioex.printStackTrace();
//            return false;
        } catch (Exception ex) {
            JOptionPane.showMessageDialog(this, "Fatal exception: " + ex,
                    "Error reading file",
                    JOptionPane.ERROR_MESSAGE);
            ex.printStackTrace();
            return false;
        }


        if (tree == null) {
            JOptionPane.showMessageDialog(this, "The file is not in a suitable format or contains no trees.",
                    "Error reading file",
                    JOptionPane.ERROR_MESSAGE);
            return false;
        }

        FlexibleTree binaryTree = new FlexibleTree(tree, true);
        binaryTree.resolveTree();
        trees.add(binaryTree);
//        if (taxa == null) {
//            taxa = binaryTree;
//        }

        getExportDataAction().setEnabled(true);

        return true;
    }

    public final void doImport() {
        final JFileChooser chooser = new JFileChooser();
        chooser.setMultiSelectionEnabled(true);

        FileNameExtensionFilter filter = new FileNameExtensionFilter("Tab or Comma delimited tables", "csv", "txt");
        chooser.setFileFilter(filter);

        final int returnVal = chooser.showOpenDialog(this);
        if (returnVal == JFileChooser.APPROVE_OPTION) {
            File[] files = chooser.getSelectedFiles();
            importMeasurementFiles(files);
        }
    }

    void importMeasurementFiles(File[] files) {
        List<DataTable<String[]>> dataTables = new ArrayList<DataTable<String[]>>();

        for (int i = 0; i < files.length; i++) {
            String fileName = files[i].getName();

            DataTable<String[]> dataTable;
            try {
                dataTable = DataTable.Text.parse(new FileReader(files[i]), true, false);
                dataTables.add(dataTable);
            } catch (IOException e) {
                JOptionPane.showMessageDialog(this, "Unable to read measurements from file, " + fileName,
                        "Error reading file",
                        JOptionPane.ERROR_MESSAGE);
            }
        }

        for (DataTable<String[]> dataTable : dataTables) {
            document.addTable(dataTable);
        }
    }


    public final void doImportLocations() {
        final JFileChooser chooser = new JFileChooser();
        chooser.setMultiSelectionEnabled(true);

        FileNameExtensionFilter filter = new FileNameExtensionFilter("BEAST log (*.log) Files", "log", "txt");
        chooser.setFileFilter(filter);

        final int returnVal = chooser.showOpenDialog(this);
        if (returnVal == JFileChooser.APPROVE_OPTION) {
            File[] files = chooser.getSelectedFiles();
            importLocationFiles(files);
        }
    }

    void importLocationFiles(File[] files) {
        LogFileTraces[] traces = new LogFileTraces[files.length];

        for (int i = 0; i < files.length; i++) {
            traces[i] = new LogFileTraces(files[i].getName(), files[i]);
        }

        processTraces(traces);
    }

    protected void processTraces(final LogFileTraces[] tracesArray) {

        final JFrame frame = this;

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

//                final Reader reader = new InputStreamReader(in);

                Thread readThread = new Thread() {
                    public void run() {
                        try {
                            traces.loadTraces(in);

                            EventQueue.invokeLater(
                                    new Runnable() {
                                        public void run() {
//                                            analyseTraceList(traces);
//                                            addTraceList(traces);
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
//                            final Reader reader = new FileReader(traces.getFile());
                            traces.loadTraces();

                            EventQueue.invokeLater(
                                    new Runnable() {
                                        public void run() {
//                                            analyseTraceList(traces);
//                                            addTraceList(traces);
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

    protected boolean writeToFile(File file) throws IOException {
        return false;  //To change body of implemented methods use File | Settings | File Templates.
    }

    protected void doExportData() {
        FileDialog dialog = new FileDialog(this,
                "Export Data File...",
                FileDialog.SAVE);

        dialog.setVisible(true);
        if (dialog.getFile() != null) {
            File file = new File(dialog.getDirectory(), dialog.getFile());

            Writer writer = null;
            try {
                writer = new PrintWriter(file);
//                analysisPanel.writeDataFile(writer);
                writer.close();
            } catch (IOException ioe) {
                JOptionPane.showMessageDialog(this, "Error writing data file: " + ioe.getMessage(),
                        "Export Error",
                        JOptionPane.ERROR_MESSAGE);
            }

        }
    }

    private void setStatusMessage() {
//        Tree tree = treesPanel.getTree();
//        if (tree != null) {
//            String message = "";
//            message += "Tree loaded, " + tree.getTaxonCount() + " taxa";
//
//            TemporalRooting tr = treesPanel.getTemporalRooting();
//            if (tr.isContemporaneous()) {
//                message += ", contemporaneous tips";
//            } else {
//                NumberFormatter nf = new NumberFormatter(3);
//                message += ", dated tips with range " + nf.format(tr.getDateRange());
//            }
//            statusLabel.setText(message);
//        }
    }

    public JComponent getExportableComponent() {

        JComponent exportable = null;
        Component comp = tabbedPane.getSelectedComponent();

        if (comp instanceof Exportable) {
            exportable = ((Exportable) comp).getExportableComponent();
        } else if (comp instanceof JComponent) {
            exportable = (JComponent) comp;
        }

        return exportable;
    }

    @Override
    public void doCopy() {
        StringWriter writer = new StringWriter();
        PrintWriter pwriter = new PrintWriter(writer);

//        for (String tip : treesPanel.getSelectedTips()) {
//            pwriter.println(tip);
//        }

        Clipboard clipboard = Toolkit.getDefaultToolkit().getSystemClipboard();
        StringSelection selection = new StringSelection(writer.toString());
        clipboard.setContents(selection, selection);
    }

    @Override
    public Action getImportAction() {
        return getImportMeasurementsAction();
    }

    @Override
    public Action getImportMeasurementsAction() {
        return importMeasurementsAction;
    }

    @Override
    public Action getImportLocationsAction() {
        return importLocationsAction;
    }

    @Override
    public Action getImportTreesAction() {
        return importTreesAction;
    }

    public Action getDeleteItemAction() {
        return getDeleteAction();
    }

    @Override
    public Action getExportAction() {
        return getExportDataAction();
    }

    @Override
    public Action getExportDataAction() {
        return exportDataAction;
    }

    @Override
    public Action getExportPDFAction() {
        return exportPDFAction;
    }

    protected AbstractAction importMeasurementsAction = new AbstractAction("Import Measurements...") {
        public void actionPerformed(ActionEvent ae) {
            doImport();
        }
    };

    protected AbstractAction importLocationsAction = new AbstractAction("Import Locations...") {
        public void actionPerformed(ActionEvent ae) {
//            doImportLocations();
        }
    };

    protected AbstractAction importTreesAction = new AbstractAction("Import Trees...") {
        public void actionPerformed(ActionEvent ae) {
//            doImportTrees();
        }
    };

    protected AbstractAction exportDataAction = new AbstractAction("Export Data...") {
        public void actionPerformed(ActionEvent ae) {
            doExportData();
        }
    };

    protected AbstractAction exportPDFAction = new AbstractAction("Export PDF...") {
        public void actionPerformed(ActionEvent ae) {
//            doExportPDF();
        }
    };

}