package dr.inference.model;

import dr.inference.loggers.LogColumn;
import dr.inference.loggers.NumberColumn;

/**
 * Created by max on 4/6/16.
 */
public class AdaptableSizeFastMatrixParameter extends FastMatrixParameter {
    public AdaptableSizeFastMatrixParameter(String id, int rowDimension, int colDimension, int maxRow, int maxCol) {
        super(id, maxRow, maxCol);
        if(maxRow < rowDimension){
            throw new RuntimeException("Row Dimension: " + rowDimension + ", is greater than Max Row Dimension: " + maxRow + " in " + getParameterName());
        }
        if(maxCol < colDimension){
            throw new RuntimeException("Column Dimension: " + colDimension + ", is greater than Max Column Dimension: " + columnDimension + " in " + getParameterName());
        }
        this.rowDimension = rowDimension;
        this.columnDimension = colDimension;

    }

    void setRowDimension(int rowDimension){
        if(rowDimension > this.rowDimension)
            fireParameterChangedEvent(-1, ChangeType.ADDED);
        else
            fireParameterChangedEvent(-1, ChangeType.REMOVED);
        this.rowDimension = rowDimension;
    }

    void setColumnDimension(int columnDimension){
        if(columnDimension > this.columnDimension)
            fireParameterChangedEvent(-1, ChangeType.ADDED);
        else
            fireParameterChangedEvent(-1, ChangeType.REMOVED);
        this.columnDimension = columnDimension;
    }

    public int getRowDimension(){
        return rowDimension;
    }

    public int getColumnDimension(){
        return columnDimension;
    }

    public int getDimension(){
        return rowDimension * columnDimension;
    }

    public void setParameterValue(int index, double value){
        int row = index % rowDimension;
        int col = index / rowDimension;

        super.setParameterValueQuietly(row, col, value);
        fireParameterChangedEvent(index, null);
    }

    public void setParameterValueQuietly(int index, double value){
        int row = index % rowDimension;
        int col = index / rowDimension;

        super.setParameterValueQuietly(row, col, value);
    }

    public double getParameterValue(int index){
        int row = index % rowDimension;
        int col = index / rowDimension;

        return super.getParameterValue(row, col);
    }

    public double[] getParameterValues(){
        double[] answer = new double[getDimension()];
        for (int i = 0; i < getDimension(); i++) {
            answer[i]= getParameterValue(i);
        }
        return answer;
    }

    int rowDimension;
    int columnDimension;
    int maxRow;
    int maxCol;

    public LogColumn[] getColumns(){
        LogColumn[] bigMatrixColumn = new ASFMPColumn[1];
        bigMatrixColumn[0] = new ASFMPColumn(getParameterName());
        return bigMatrixColumn;
    }

    private class ASFMPColumn extends NumberColumn {

        public ASFMPColumn(String label) {
            super(label);
        }

        protected String getFormattedValue(){
            String fullMatrix = "{";
            for (int i = 0; i <getRowDimension() ; i++) {
                fullMatrix += " { ";
                for (int j = 0; j < getColumnDimension(); j++) {
                    fullMatrix += formatValue(getParameterValue(i,j));
                    if(j != getColumnDimension() - 1){
                        fullMatrix += ", ";
                    }
                }
                if(i == getRowDimension()-1)
                    fullMatrix += " } ";
                else
                    fullMatrix += " },";
            }
            fullMatrix += "}";
            return fullMatrix;
        }

        @Override
        public double getDoubleValue() {
            return 0;
        }
    }
}

