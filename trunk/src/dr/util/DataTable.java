package dr.util;

import java.io.*;
import java.util.*;

/**
 * @author Andrew Rambaut
 * @version $Id$
 */
public interface DataTable<T> {
    int getColumnCount();

    int getRowCount();

    String[] getColumnLabels();
    String[] getRowLabels();

    T getColumn(int columnIndex);
    T getRow(int rowIndex);

    T[] getData();

    public class Double implements DataTable<double[]> {

        private Double(Reader source) throws IOException {
            BufferedReader reader = new BufferedReader(source);

            String line = reader.readLine();
            if (line == null) {
                throw new IllegalArgumentException("Empty file");
            }

            List<String> columnLabels = new ArrayList<String>();

            StringTokenizer tokenizer = new StringTokenizer(line, "\t");

            // potentially the first token is the name of the row labels
            String name = tokenizer.nextToken();

            while (tokenizer.hasMoreTokens()) {
                String label = tokenizer.nextToken();
                columnLabels.add(label);
            }

            List<String> rowLabels = new ArrayList<String>();
            List<double[]> rows = new ArrayList<double[]>();

            int rowIndex = 1;

            line = reader.readLine();
            while (line != null) {
                tokenizer = new StringTokenizer(line, "\t");

                String label = tokenizer.nextToken();
                rowLabels.add(label);

                double[] row = new double[this.columnLabels.length];

                int columnIndex = 0;
                while (tokenizer.hasMoreTokens()) {
                    String value = tokenizer.nextToken();
                    try {
                        row[columnIndex] = java.lang.Double.valueOf(value);
                    } catch (NumberFormatException nfe) {
                        throw new IllegalArgumentException("Non numerical value at row " + (rowIndex + 1) +
                                ", column " + (columnIndex + 1));
                    }

                    columnIndex ++;
                }
                if (columnIndex != this.columnLabels.length) {
                    throw new IllegalArgumentException("Wrong number of values on row " + (rowIndex + 1) +
                            ", expecting " + this.columnLabels.length + " but actually " + columnIndex);
                }

                rows.add(row);

                line = reader.readLine();
                rowIndex++;
            }

            this.rowLabels = rowLabels.toArray(this.rowLabels);
            data = rows.toArray(data);

            this.columnLabels = columnLabels.toArray(this.columnLabels);
        }

        public int getColumnCount() {
            return columnLabels.length;
        }

        public int getRowCount() {
            return rowLabels.length;
        }

        public String[] getColumnLabels() {
            return columnLabels;
        }

        public String[] getRowLabels() {
            return rowLabels;
        }

        public double[] getColumn(final int columnIndex) {
            double[] column = new double[data.length];
            for (int i = 0; i < data.length; i++) {
                column[i] = data[i][columnIndex];
            }
            return column;
        }

        public double[] getRow(final int rowIndex) {
            return data[rowIndex];
        }

        public double[][] getData() {
            return data;
        }

        private String[] columnLabels;
        private String[] rowLabels;

        private double[][] data;

        public static DataTable<double []> parse(Reader source) throws IOException {
            return new DataTable.Double(source);
        }
    }



}
