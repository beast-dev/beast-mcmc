package dr.util;

/**
 * Tabular data provider
 *
 *  A very modest start. will evolve further according to needs.
 */
public abstract class TabularData {
    public abstract int nColumns();

    public abstract String columnName(int nColumn);

    public abstract int nRows();

    public abstract Object data(int nRow, int nColumn);

    public int getColumn(String name) {
        for(int n = 0; n < nColumns(); ++n) {
            if( columnName(n).equals(name) ) {
                return n;
            }
        }
        return -1;
    }
}

