package dr.inference.model;

/**
 * Created by max on 4/6/16.
 */
public class AdaptableSizeFastMatrixParameter extends FastMatrixParameter {
    public AdaptableSizeFastMatrixParameter(String id, int rowDimension, int colDimension, int maxCol, int maxRow) {
        super(id, maxRow, maxCol);
        this.rowDimension = rowDimension;
        this.columnDimension = colDimension;

    }

    void setRowDimension(int rowDimension){
        this.rowDimension = rowDimension;
    }

    void setColumnDimension(int columnDimension){
        this.columnDimension = columnDimension;
    }

    public int getRowDimension(){
        return rowDimension;
    }

    public int getColumnDimension(){
        return columnDimension;
    }

    public int getDimension(){
        return rowDimension*columnDimension;
    }

    public void setParameterValue(int index, double value){
        int row = index / rowDimension;
        int col = index % rowDimension;

        super.setParameterValue(row, col, value);
    }

    public void setParameterValueQuietly(int index, double value){
        int row = index / rowDimension;
        int col = index % rowDimension;

        super.setParameterValueQuietly(row, col, value);
    }

    public double getParameterValue(int index){
        int row = index / rowDimension;
        int col = index % rowDimension;

        return super.getParameterValue(row, col);
    }

    int rowDimension;
    int columnDimension;
    int maxRow;
    int maxCol;
}

