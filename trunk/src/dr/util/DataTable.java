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

        private Double(Reader source, boolean hasColumnLabels, boolean hasRowLabels) throws IOException {
            BufferedReader reader = new BufferedReader(source);

            String line = reader.readLine();
            if (line == null) {
                throw new IllegalArgumentException("Empty file");
            }

            int columnCount = -1;

            if (hasColumnLabels) {
                List<String> columnLabels = new ArrayList<String>();

                StringTokenizer tokenizer = new StringTokenizer(line, "\t");

                if (hasRowLabels && !line.startsWith("\t")) {
                    // potentially the first token is the name of the row labels
                    String name = tokenizer.nextToken();
                }

                while (tokenizer.hasMoreTokens()) {
                    String label = tokenizer.nextToken();
                    columnLabels.add(label);
                }

                this.columnLabels = new String[columnLabels.size()];
                columnLabels.toArray(this.columnLabels);

                columnCount = columnLabels.size();

                line = reader.readLine();
            }

            List<String> rowLabels = new ArrayList<String>();
            List<double[]> rows = new ArrayList<double[]>();

            int rowIndex = 1;

            while (line != null) {
                StringTokenizer tokenizer = new StringTokenizer(line, "\t");

                if (columnCount == -1) {
                    columnCount = tokenizer.countTokens();
                    if (hasRowLabels) {
                        columnCount --;
                    }
                }

                if (hasRowLabels) {
                    String label = tokenizer.nextToken();
                    rowLabels.add(label);
                }

                double[] row = new double[columnCount];

                int columnIndex = 0;
                while (tokenizer.hasMoreTokens()) {
                    String value = tokenizer.nextToken();
                    try {
                        row[columnIndex] = java.lang.Double.valueOf(value);
                    } catch (NumberFormatException nfe) {
                        row[columnIndex] = java.lang.Double.NaN;
//                        throw new IllegalArgumentException("Non numerical value at row " + (rowIndex + 1) +
//                                ", column " + (columnIndex + 1));
                    }

                    columnIndex ++;
                }
                if (columnIndex != columnCount) {
                    throw new IllegalArgumentException("Wrong number of values on row " + (rowIndex + 1) +
                            ", expecting " + columnCount + " but actually " + columnIndex);
                }

                rows.add(row);

                line = reader.readLine();
                rowIndex++;
            }

            if (hasRowLabels) {
                this.rowLabels = new String[rowLabels.size()];
                rowLabels.toArray(this.rowLabels);
            }

            data = new double[rows.size()][];
            rows.toArray(data);

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
            return new DataTable.Double(source, true, true);
        }

        public static DataTable<double []> parse(Reader source, boolean columnLabels, boolean rowLabels) throws IOException {
            return new DataTable.Double(source, columnLabels, rowLabels);
        }
    }

    class Text implements DataTable<String[]> {

        private Text(Reader source, boolean hasColumnLabels, boolean hasRowLabels) throws IOException {
            BufferedReader reader = new BufferedReader(source);

            String line = reader.readLine();
            if (line == null) {
                throw new IllegalArgumentException("Empty file");
            }

            int columnCount = -1;

            if (hasColumnLabels) {
                List<String> columnLabels = new ArrayList<String>();

                StringTokenizer tokenizer = new StringTokenizer(line, "\t");

                if (hasRowLabels && !line.startsWith("\t")) {
                    // potentially the first token is the name of the row labels
                    String name = tokenizer.nextToken();
                }

                while (tokenizer.hasMoreTokens()) {
                    String label = tokenizer.nextToken();
                    columnLabels.add(label);
                }

                this.columnLabels = new String[columnLabels.size()];
                columnLabels.toArray(this.columnLabels);

                columnCount = columnLabels.size();

                line = reader.readLine();
            }

            List<String> rowLabels = new ArrayList<String>();
            List<String[]> rows = new ArrayList<String[]>();
            int rowIndex = 1;

            while (line != null) {
                StringTokenizer tokenizer = new StringTokenizer(line, "\t");

                if (columnCount == -1) {
                    columnCount = tokenizer.countTokens();
                    if (hasRowLabels) {
                        columnCount --;
                    }
                }

                if (hasRowLabels) {
                    String label = tokenizer.nextToken();
                    rowLabels.add(label);
                }

                String[] row = new String [columnCount];

                int columnIndex = 0;
                while (tokenizer.hasMoreTokens()) {
                    row[columnIndex] = tokenizer.nextToken().trim();

                    columnIndex ++;
                }
                if (columnIndex != columnCount) {
                    throw new IllegalArgumentException("Wrong number of values on row " + (rowIndex + 1) +
                            ", expecting " + columnCount + " but actually " + columnIndex);
                }

                rows.add(row);

                line = reader.readLine();
                rowIndex++;
            }

            if (hasRowLabels) {
                this.rowLabels = new String[rowLabels.size()];
                rowLabels.toArray(this.rowLabels);
            }

            data = new String[rows.size()][];
            rows.toArray(data);

        }

        public int getColumnCount() {
            return data[0].length;
        }

        public int getRowCount() {
            return data.length;
        }

        public String[] getColumnLabels() {
            return columnLabels;
        }

        public String[] getRowLabels() {
            return rowLabels;
        }

        public String[] getColumn(final int columnIndex) {
            String[] column = new String[data.length];
            for (int i = 0; i < data.length; i++) {
                column[i] = data[i][columnIndex];
            }
            return column;
        }

        public String[] getRow(final int rowIndex) {
            return data[rowIndex];
        }

        public String[][] getData() {
            return data;
        }

        private String[] columnLabels;
        private String[] rowLabels;

        private String[][] data;

        public static DataTable<String []> parse(Reader source) throws IOException {
            return new DataTable.Text(source, true, true);
        }

        public static DataTable<String []> parse(Reader source, boolean columnLabels, boolean rowLabels) throws IOException {
            return new DataTable.Text(source, columnLabels, rowLabels);
        }
    }


}
