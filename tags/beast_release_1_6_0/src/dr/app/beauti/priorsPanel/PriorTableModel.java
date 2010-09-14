package dr.app.beauti.priorsPanel;

import dr.app.beauti.options.Parameter;

import javax.swing.table.AbstractTableModel;

/**
 * @author Alexei Drummond
 * @author Andrew Rambaut
 */
class PriorTableModel extends AbstractTableModel {

    private static final long serialVersionUID = -8864178122484971872L;

    String[] columnNames = {"Parameter", "Prior", "Bound", "Description"};
    private PriorsPanel priorsPanel;

    public PriorTableModel(PriorsPanel priorsPanel) {
        this.priorsPanel = priorsPanel;
    }

    public int getColumnCount() {
        return columnNames.length;
    }

    public int getRowCount() {
        return priorsPanel.parameters.size();
    }

    public Object getValueAt(int row, int col) {
        Parameter param = priorsPanel.parameters.get(row);
        switch (col) {
            case 0:
                return param.getName();
            case 1:
                return param.priorType.getPriorString(param);
            case 2:
                return param.priorType.getPriorBoundString(param);
            case 3:
                return param.getDescription();
        }
        return null;
    }

    public String getColumnName(int column) {
        return columnNames[column];
    }

    public Class getColumnClass(int c) {
        return getValueAt(0, c).getClass();
    }

    public boolean isCellEditable(int row, int col) {
        return col == 1;
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
