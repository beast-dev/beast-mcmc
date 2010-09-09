package pyromania.sourcespanel;

import jam.framework.Exportable;
import jam.panels.ActionPanel;
import pyromania.database.Pyromania;
import pyromania.PyromaniaFrame;
import pyromania.database.Source;
import pyromania.gui.DeleteActionResponder;
import pyromania.gui.MultiLineTableCellRenderer;
import pyromania.gui.table.TableEditorStopper;

import javax.swing.*;
import javax.swing.event.ListSelectionEvent;
import javax.swing.event.ListSelectionListener;
import javax.swing.plaf.BorderUIResource;
import javax.swing.table.AbstractTableModel;
import javax.swing.table.TableColumn;
import java.awt.*;
import java.awt.datatransfer.*;
import java.awt.event.*;
import java.io.File;
import java.io.IOException;
import java.util.*;

/**
 * @author Andrew Rambaut
 * @version $Id$
 */
public class SourcesPanel extends JPanel implements DeleteActionResponder, Exportable {
    private JScrollPane scrollPane = new JScrollPane();
    private JTable dataTable = null;
    private DataTableModel dataTableModel = null;

    private PyromaniaFrame frame = null;

    private SourcesSettingsDialog sourcesSettingsDialog = null;
    private final Pyromania document;

    private java.util.List<Source> sources;

    public SourcesPanel(final PyromaniaFrame parent, final Pyromania document, final Action addDataAction) {

        this.frame = parent;
        this.document = document;

        dataTableModel = new DataTableModel();
        dataTable = new JTable(dataTableModel);

        dataTable.getTableHeader().setReorderingAllowed(false);
//        dataTable.getTableHeader().setDefaultRenderer(
//                new HeaderRenderer(SwingConstants.LEFT, new Insets(0, 4, 0, 4)));

        TableColumn col = dataTable.getColumnModel().getColumn(0);
        col.setCellRenderer(new MultiLineTableCellRenderer());

        dataTable.setRowHeight(dataTable.getRowHeight() * 2);

        dataTable.setDragEnabled(false);
        dataTable.setTransferHandler(new FSTransfer());

        TableEditorStopper.ensureEditingStopWhenTableLosesFocus(dataTable);

        dataTable.getSelectionModel().addListSelectionListener(new ListSelectionListener() {
            public void valueChanged(ListSelectionEvent evt) {
                selectionChanged();
            }
        });

        dataTable.addMouseListener(new MouseAdapter() {
            public void mouseClicked(MouseEvent e) {
                if (e.getClickCount() == 2) {
                    editSelection();
                }
            }
        });

        scrollPane = new JScrollPane(dataTable,
                JScrollPane.VERTICAL_SCROLLBAR_ALWAYS,
                JScrollPane.HORIZONTAL_SCROLLBAR_ALWAYS);
        scrollPane.setOpaque(false);

//        JToolBar toolBar1 = new JToolBar();
//        toolBar1.setFloatable(false);
//        toolBar1.setOpaque(false);
//        toolBar1.setLayout(new FlowLayout(FlowLayout.LEFT, 0, 0));

//        JButton button = new JButton(unlinkModelsAction);
//        unlinkModelsAction.setEnabled(false);
//        PanelUtils.setupComponent(button);
//        toolBar1.add(button);


        ActionPanel actionPanel1 = new ActionPanel(true);
        actionPanel1.setAddAction(addDataAction);
        actionPanel1.setRemoveAction(removeAction);

        removeAction.setEnabled(false);

        JPanel controlPanel1 = new JPanel(new FlowLayout(FlowLayout.LEFT));
        controlPanel1.setOpaque(false);
        controlPanel1.add(actionPanel1);

        setOpaque(false);
        setBorder(new BorderUIResource.EmptyBorderUIResource(new Insets(12, 12, 12, 12)));
        setLayout(new BorderLayout(0, 0));
//        add(toolBar1, BorderLayout.NORTH);
        add(scrollPane, BorderLayout.CENTER);
        add(controlPanel1, BorderLayout.SOUTH);

//        document.addListener(new SpreadDocument.Listener() {
//            public void dataChanged() {
//                dataTableModel.fireTableDataChanged();
//            }
//
//            public void settingsChanged() { }
//        });

        sources = document.getSources();
    }

    public void selectionChanged() {
        int[] selRows = dataTable.getSelectedRows();
        boolean hasSelection = (selRows != null && selRows.length != 0);
//        removeAction.setEnabled(hasSelection);

//        frame.setRemoveActionEnabled(this, hasSelection);
    }

    public JComponent getExportableComponent() {
        return dataTable;
    }

    public void editSelection() {
        int selRow = dataTable.getSelectedRow();
        if (selRow >= 0) {
            Source source = document.getSources().get(selRow);
            editSettings(source);
        }
    }

    private void editSettings(Source source) {
        if (sourcesSettingsDialog == null) {
            sourcesSettingsDialog = new SourcesSettingsDialog(frame);
        }

        int result = sourcesSettingsDialog.showDialog(source);

        if (result != JOptionPane.CANCEL_OPTION) {
            sourcesSettingsDialog.getReadsFile(); // force update of builder settings
            document.fireDataChanged();
        }
    }

    public void delete() {
        int[] selRows = dataTable.getSelectedRows();
        Set<Source> dataToRemove = new HashSet<Source>();
        for (int row : selRows) {
            dataToRemove.add(document.getSources().get(row));
        }

        // TODO: would probably be a good idea to check if the user wants to remove the last file
        document.getSources().removeAll(dataToRemove);
        document.fireDataChanged();

        if (document.getSources().size() == 0) {
            // all data partitions removed so reset the taxa
//            frame.setStatusMessage("No data loaded");
        }
    }

    public void selectAll() {
        dataTable.selectAll();
    }

    class DataTableModel extends AbstractTableModel {

        private static final long serialVersionUID = -6707994233020715574L;
        String[] columnNames = {"Reads"};

        public DataTableModel() {
        }

        public int getColumnCount() {
            return columnNames.length;
        }

        public int getRowCount() {
            return sources.size();
        }

        public Object getValueAt(int row, int col) {
            Source source = sources.get(row);
            switch (col) {
                case 0:
                    return source;
                default:
                    throw new IllegalArgumentException("unknown column, " + col);
            }
        }

        public boolean isCellEditable(int row, int col) {
            return false;
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

    public class FSTransfer extends TransferHandler {
        public boolean importData(JComponent comp, Transferable t) {
            // Make sure we have the right starting points
            if (!(comp instanceof JTable)) {
                return false;
            }
            if (!t.isDataFlavorSupported(DataFlavor.javaFileListFlavor)) {
                return false;
            }

            try {
                java.util.List data = (java.util.List) t.getTransferData(DataFlavor.javaFileListFlavor);
                Iterator i = data.iterator();
                while (i.hasNext()) {
                    File f = (File) i.next();
                    frame.importDataFile(f);
                }
                return true;
            } catch (UnsupportedFlavorException ufe) {
                System.err.println("Ack! we should not be here.\nBad Flavor.");
            } catch (IOException ioe) {
                System.out.println("Something failed during import:\n" + ioe);
            }
            return false;
        }

        public boolean canImport(JComponent comp, DataFlavor[] transferFlavors) {
            if (comp instanceof JTable) {
                for (int i = 0; i < transferFlavors.length; i++) {
                    if (!transferFlavors[i].equals(DataFlavor.javaFileListFlavor)) {
                        return false;
                    }
                }
                return true;
            }
            return false;
        }
    }

    public Action getDeleteAction() {
        return removeAction;
    }

    private Action removeAction = new RemoveDataAction();
    private class RemoveDataAction extends AbstractAction {
        public RemoveDataAction() {
            super("Remove");
            setToolTipText("Use this button to remove selected reads data from the table");
        }

        public void actionPerformed(ActionEvent ae) {
            delete();
        }
    }


}