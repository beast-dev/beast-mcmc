package dr.inference.model;

import dr.inference.loggers.LogColumn;
import dr.inference.loggers.NumberColumn;

/**
 * Created by max on 4/6/16.
 */
public class AdaptableSizeFastMatrixParameter extends FastMatrixParameter {
    public AdaptableSizeFastMatrixParameter(String id, int rowDimension, int colDimension, int maxRow, int maxCol, double startingValue, boolean lowerTriangle) {
        super(id, maxRow, maxCol, startingValue);
        if(maxRow < rowDimension){
            throw new RuntimeException("Row Dimension: " + rowDimension + ", is greater than Max Row Dimension: " + maxRow + " in " + getParameterName());
        }
        if(maxCol < colDimension){
            throw new RuntimeException("Column Dimension: " + colDimension + ", is greater than Max Column Dimension: " + columnDimension + " in " + getParameterName());
        }
        this.rowDimension = rowDimension;
        this.columnDimension = colDimension;
        this.maxRow = maxRow;
        this.maxCol = maxCol;
        this.storedColumnDimension = rowDimension;
        this.storedColumnDimension = columnDimension;
        this.lowerTriangle = lowerTriangle;

    }

    public void setRowDimension(int rowDimension){
        int oldRowDimension = this.rowDimension;
        this.rowDimension = rowDimension;
        if(rowDimension > maxRow){
            throw new RuntimeException("Row Dimension Larger Than Maximum");
        }
        if(rowDimension > oldRowDimension)
            fireParameterChangedEvent(-1, ChangeType.ADDED);
        else
            fireParameterChangedEvent(-1, ChangeType.REMOVED);
    }

    public void setColumnDimension(int columnDimension){
        int oldColumnDimension = this.columnDimension;
        if(columnDimension > maxCol){
            throw new RuntimeException("Column Dimension Larger Than Maximum");
        }
        this.columnDimension = columnDimension;
        if(columnDimension > oldColumnDimension)
            fireParameterChangedEvent(-1, ChangeType.ADDED);
        else
            fireParameterChangedEvent(-1, ChangeType.REMOVED);

    }

    @Override
    protected void storeValues(){
        super.storeValues();
        storedRowDimension = rowDimension;
        storedColumnDimension = columnDimension;
    }

    protected void restoreValues(){
        super.restoreValues();
        rowDimension = storedRowDimension;
        columnDimension = storedColumnDimension;
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

    public int getMaxRowDimension(){return maxRow;}

    public int getMaxColumnDimension(){return maxCol;}

    public void setParameterValue(int index, double value) {
        int row = index % rowDimension;
        int col = index / rowDimension;

        if (row >= col || !lowerTriangle){
            super.setParameterValueQuietly(row, col, value);
            fireParameterChangedEvent(index, ChangeType.VALUE_CHANGED);
        }
    }

    public void setParameterValue(int row, int column, double value){
        if(row >= column || !lowerTriangle){
            super.setParameterValueQuietly(row, column, value);
            fireParameterChangedEvent(getRowDimension() * column + row, ChangeType.VALUE_CHANGED);
        }

    }

    public void setParameterValueQuietly(int index, double value){
        int row = index % rowDimension;
        int col = index / rowDimension;

        if(row >= col || !lowerTriangle) {
            super.setParameterValueQuietly(row, col, value);
        }
    }

    public void setParameterValueQuietly(int row, int column, double value){
        if(row >= column || !lowerTriangle){
            super.setParameterValueQuietly(row, column, value);
            }

    }

    public double getParameterValue(int index){
        int row = index % rowDimension;
        int col = index / rowDimension;

        if(row >= col || !lowerTriangle) {
                return super.getParameterValue(row, col);
            }
            else return 0;
    }

    public double getParameterValue(int row, int col){
        if(row >= col || !lowerTriangle) {
            return super.getParameterValue(row, col);
        }
        else return 0;
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
    int storedRowDimension;
    int storedColumnDimension;
    int maxRow;
    int maxCol;
    boolean lowerTriangle;

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

