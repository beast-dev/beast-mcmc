package dr.app.tracer.traces;

import dr.gui.chart.JChartPanel;
import org.virion.jam.framework.Exportable;

import javax.swing.*;
import javax.swing.table.DefaultTableModel;
import javax.swing.table.JTableHeader;
import java.awt.*;

/**
 * A panel that displays correlation plots of 2 traces in integer and category type
 *
 * @author Andrew Rambaut
 * @author Alexei Drummond
 * @author Walter Xie
 * @version $Id: TablePanel.java,v 1.1.1.2 2006/04/25 23:00:09 rambaut Exp $
 */
public class TablePanel extends JChartPanel implements Exportable {

    JTable dataTable = null;
    JScrollPane scrollPane;
    DataTableModel dataTableModel;
    ListModel lm;
    Object[] rowNames = new Object[0];
    Object[] columnNames = new Object[0];
    double[][] data;

    public TablePanel(String title, String xAxisTitle, String yAxisTitle) {
        super(null, title, xAxisTitle, yAxisTitle);

        lm = new RowHeaderModel();
        dataTableModel = new DataTableModel(); //(int rowCount, int columnCount)
        dataTable = new JTable(dataTableModel);
        dataTable.setAutoResizeMode(JTable.AUTO_RESIZE_OFF);

        scrollPane = new JScrollPane(dataTable,
                JScrollPane.VERTICAL_SCROLLBAR_ALWAYS,
                JScrollPane.HORIZONTAL_SCROLLBAR_ALWAYS);
        scrollPane.setOpaque(false);
        add(scrollPane, "Table");

    }

    public void setTable(Object[] rowNames, Object[] columnNames, double[][] data) {
        this.rowNames = rowNames;
        this.columnNames = columnNames;
        this.data = data;

        lm = new RowHeaderModel();

        dataTableModel.fireTableDataChanged();

        JList rowHeader = new JList(lm);
        rowHeader.setFixedCellWidth(50);
        rowHeader.setFixedCellHeight(dataTable.getRowHeight() + dataTable.getRowMargin());
        //                           + table.getIntercellSpacing().height);
        rowHeader.setCellRenderer(new RowHeaderRenderer(dataTable));

        scrollPane.setRowHeaderView(rowHeader);
//        removeAll();

        validate();
        repaint();
    }

    public JComponent getExportableComponent() {
        return this;
    }

    class RowHeaderModel extends AbstractListModel {
        public int getSize() {
            return rowNames.length;
        }

        public Object getElementAt(int index) {
            return rowNames[index];
        }
    }

    class RowHeaderRenderer extends JLabel implements ListCellRenderer {

        RowHeaderRenderer(JTable table) {
            JTableHeader header = table.getTableHeader();
            setOpaque(true);
            setBorder(UIManager.getBorder("TableHeader.cellBorder"));
            setHorizontalAlignment(CENTER);
            setForeground(header.getForeground());
            setBackground(header.getBackground());
            setFont(header.getFont());
        }

        public Component getListCellRendererComponent(JList list, Object value,
                                                      int index, boolean isSelected, boolean cellHasFocus) {
            setText((value == null) ? "" : value.toString());
            return this;
        }
    }


    class DataTableModel extends DefaultTableModel {

//            public DataTableModel(int rowCount, int columnCount) {
//                super(rowCount, columnCount);
//            }

        public DataTableModel() {
        }

        public int getColumnCount() {
            return columnNames.length;
        }

        public int getRowCount() {
            return lm.getSize();
        }

        public Object getValueAt(int row, int col) {
            return data[row][col];
        }

        public boolean isCellEditable(int row, int col) {
            return false;
        }

        public String getColumnName(int column) {
            return columnNames[column].toString();
        }

//            public Class getColumnClass(int c) {
//                if (getRowCount() == 0) {
//                    return Object.class;
//                }
//                return getValueAt(0, c).getClass();
//            }

//            public String toString() {
//                StringBuffer buffer = new StringBuffer();
//
//                buffer.append(getColumnName(0));
//                for (int j = 1; j < getColumnCount(); j++) {
//                    buffer.append("\t");
//                    buffer.append(getColumnName(j));
//                }
//                buffer.append("\n");
//
//                for (int i = 0; i < getRowCount(); i++) {
//                    buffer.append(getValueAt(i, 0));
//                    for (int j = 1; j < getColumnCount(); j++) {
//                        buffer.append("\t");
//                        buffer.append(getValueAt(i, j));
//                    }
//                    buffer.append("\n");
//                }
//
//                return buffer.toString();
//            }
    }

}

