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

    int rowDimension;
    int columnDimension;
    int maxRow;
    int maxCol;
}

