/*
 * DataTable.java
 *
 * Copyright (c) 2002-2015 Alexei Drummond, Andrew Rambaut and Marc Suchard
 *
 * This file is part of BEAST.
 * See the NOTICE file distributed with this work for additional
 * information regarding copyright ownership and licensing.
 *
 * BEAST is free software; you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as
 * published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.
 *
 *  BEAST is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with BEAST; if not, write to the
 * Free Software Foundation, Inc., 51 Franklin St, Fifth Floor,
 * Boston, MA  02110-1301  USA
 */

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

        private Double(Reader source, boolean hasColumnLabels, boolean hasRowLabels, boolean isCSV) throws IOException {
            this(source, hasColumnLabels, hasRowLabels, false, isCSV);
        }

        private Double(Reader source, boolean hasColumnLabels, boolean hasRowLabels, boolean includeFirstColumnLabel, boolean isCSV) throws IOException {
            BufferedReader reader = new BufferedReader(source);

            String line = reader.readLine();
            if (line == null) {
                throw new IllegalArgumentException("Empty file");
            }

            int columnCount = -1;

            String delim = (isCSV ? "," : "\t");

            if (hasColumnLabels) {
                List<String> columnLabels = new ArrayList<String>();

                StringTokenizer tokenizer = new StringTokenizer(line, delim);

                if (hasRowLabels && !line.startsWith(delim)) {
                    // potentially the first token is the name of the row labels
                    String name = tokenizer.nextToken().trim();
                    if (includeFirstColumnLabel) {
                        columnLabels.add(name);
                    }
                }

                while (tokenizer.hasMoreTokens()) {
                    String label = tokenizer.nextToken().trim();
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
                StringTokenizer tokenizer = new StringTokenizer(line, delim);

                if (columnCount == -1) {
                    columnCount = tokenizer.countTokens();
                    if (hasRowLabels) {
                        columnCount --;
                    }
                }

                if (hasRowLabels) {
                    String label = tokenizer.nextToken().trim();
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
                if (columnIndex != columnCount - (includeFirstColumnLabel ? 1 : 0)) {
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
            return new DataTable.Double(source, true, true, false);
        }

        public static DataTable<double []> parse(Reader source, boolean isCSV) throws IOException {
            return new DataTable.Double(source, true, true, isCSV);
        }

        public static DataTable<double []> parse(Reader source, boolean columnLabels, boolean rowLabels) throws IOException {
            return new DataTable.Double(source, columnLabels, rowLabels, false);
        }

        public static DataTable<double []> parse(Reader source, boolean columnLabels, boolean rowLabels, boolean isCSV) throws IOException {
            return new DataTable.Double(source, columnLabels, rowLabels, isCSV);
        }
    }

    class Text implements DataTable<String[]> {

        private Text(Reader source, boolean hasColumnLabels, boolean hasRowLabels, boolean isCSV) throws IOException {
            this(source, hasColumnLabels, hasRowLabels, false, isCSV);
        }

        private Text(Reader source, boolean hasColumnLabels, boolean hasRowLabels, boolean includeFirstColumnLabel, boolean isCSV) throws IOException {
            BufferedReader reader = new BufferedReader(source);

            String delim = (isCSV ? "," : "\t");

            String line = reader.readLine();
            if (line == null) {
                throw new IllegalArgumentException("Empty file");
            }

            int columnCount = -1;

            if (hasColumnLabels) {
                List<String> columnLabels = new ArrayList<String>();

                StringTokenizer tokenizer = new StringTokenizer(line, delim);

                if (hasRowLabels && !line.startsWith(delim)) {
                    // potentially the first token is the name of the row labels
                    String name = tokenizer.nextToken().trim();
                    if (includeFirstColumnLabel) {
                        columnLabels.add(name);
                    }
                }

                while (tokenizer.hasMoreTokens()) {
                    String label = tokenizer.nextToken().trim();
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
                StringTokenizer tokenizer = new StringTokenizer(line, delim);

                if (columnCount == -1) {
                    columnCount = tokenizer.countTokens();
                    if (hasRowLabels) {
                        columnCount --;
                    }
                }

                if (hasRowLabels) {
                    String label = tokenizer.nextToken().trim();
                    rowLabels.add(label);
                }

                String[] row = new String [columnCount];

                int columnIndex = 0;
                while (tokenizer.hasMoreTokens()) {
                    row[columnIndex] = tokenizer.nextToken().trim();

                    columnIndex ++;
                }
                if (columnIndex != columnCount - (includeFirstColumnLabel ? 1 : 0)) {
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
            return new DataTable.Text(source, true, true, false);
        }

        public static DataTable<String []> parse(Reader source, boolean isCSV) throws IOException {
            return new DataTable.Text(source, true, true, isCSV);
        }

        public static DataTable<String []> parse(Reader source, boolean columnLabels, boolean rowLabels) throws IOException {
            return new DataTable.Text(source, columnLabels, rowLabels, false);
        }

        public static DataTable<String []> parse(Reader source, boolean columnLabels, boolean rowLabels, boolean isCSV) throws IOException {
            return new DataTable.Text(source, columnLabels, rowLabels, isCSV);
        }

        public static DataTable<String []> parse(Reader source, boolean columnLabels, boolean rowLabels, boolean firstColumnLabel, boolean isCSV) throws IOException {
            return new DataTable.Text(source, columnLabels, rowLabels, firstColumnLabel, isCSV);
        }
    }


}
