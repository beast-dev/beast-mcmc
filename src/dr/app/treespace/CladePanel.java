package dr.app.treespace;

import dr.app.gui.table.MultiLineTableCellRenderer;
import dr.app.gui.table.TableEditorStopper;
import jam.panels.OptionsPanel;

import javax.swing.*;
import javax.swing.event.ListSelectionEvent;
import javax.swing.event.ListSelectionListener;
import javax.swing.table.AbstractTableModel;
import javax.swing.table.TableColumn;
import java.awt.*;
import java.awt.event.MouseAdapter;
import java.awt.event.MouseEvent;

/**
 * @author Andrew Rambaut
 * @version $Id$
 */
public class CladePanel extends JPanel {

    private static final long serialVersionUID = -3710586474593827540L;

    private final TreeSpaceFrame frame;
    private final TreeSpaceDocument document;
    private final CladeTableModel cladeTableModel;
    private final JTable cladeTable;

    public CladePanel(final TreeSpaceFrame parent, final TreeSpaceDocument document) {

        this.frame = parent;
        this.document = document;

        cladeTableModel = new CladeTableModel();
        cladeTable = new JTable(cladeTableModel);

        cladeTable.getTableHeader().setReorderingAllowed(false);

        TableEditorStopper.ensureEditingStopWhenTableLosesFocus(cladeTable);

        cladeTable.getSelectionModel().addListSelectionListener(new ListSelectionListener() {
            public void valueChanged(ListSelectionEvent evt) {
                selectionChanged();
            }
        });

        cladeTable.addMouseListener(new MouseAdapter() {
            public void mouseClicked(MouseEvent e) {
                if (e.getClickCount() == 2) {
                    editSelection();
                }
            }
        });

        JScrollPane scrollPane = new JScrollPane(cladeTable,
                JScrollPane.VERTICAL_SCROLLBAR_ALWAYS,
                JScrollPane.HORIZONTAL_SCROLLBAR_ALWAYS);
        scrollPane.setOpaque(false);

        setLayout(new BorderLayout());

        JPanel panel = new JPanel(new BorderLayout());

        JSplitPane splitPane = new JSplitPane(JSplitPane.HORIZONTAL_SPLIT, scrollPane, panel);
        splitPane.setDividerLocation(240);

        add(splitPane, BorderLayout.CENTER);

        document.addListener(new TreeSpaceDocument.Listener() {
            public void dataChanged() {
                cladeTableModel.fireTableDataChanged();
            }

            public void settingsChanged() {
            }
        });

    }

    private void selectionChanged() {
    }

    private void editSelection() {
    }

    class CladeTableModel extends AbstractTableModel {

        String[] columnNames = {"Clade", "Frequency"};

        public CladeTableModel() {
        }

        public int getColumnCount() {
            return columnNames.length;
        }

        public int getRowCount() {
            return document.getCladeSystem().getClades().size();
        }

        public Object getValueAt(int row, int col) {
            CladeSystem.Clade clade = document.getCladeSystem().getClades().get(row);
            switch (col) {
                case 0:
                    return row + 1;
                case 1:
                    return clade.getCredibility();
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

}