/*
 * CheckpointableMCMCOperator.java
 *
 * Copyright © 2002-2024 the BEAST Development Team
 * http://beast.community/about
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
 *
 */

package dr.inference.operators;

import java.io.PrintStream;

/**
 * Optional extension point for operators with internal state beyond the generic
 * accept/reject and adaptable-parameter fields saved by the main checkpoint.
 */
public interface CheckpointableMCMCOperator {

    String getCheckpointStateType();

    void writeCheckpointState(CheckpointStateWriter writer);

    void readCheckpointState(CheckpointStateReader reader);

    class CheckpointStateWriter {
        private final PrintStream out;

        public CheckpointStateWriter(PrintStream out) {
            this.out = out;
        }

        public void writeInt(int value) {
            out.print("\t");
            out.print(value);
        }

        public void writeVector(double[] values) {
            if (values == null) {
                throw new IllegalArgumentException("Cannot write null checkpoint vector");
            }
            out.print("\t");
            out.print(values.length);
            for (double value : values) {
                out.print("\t");
                out.print(value);
            }
        }

        public void writeMatrix(double[][] values) {
            if (values == null) {
                throw new IllegalArgumentException("Cannot write null checkpoint matrix");
            }
            out.print("\t");
            out.print(values.length);
            int columns = 0;
            if (values.length > 0) {
                if (values[0] == null) {
                    throw new IllegalArgumentException("Cannot write ragged checkpoint matrix");
                }
                columns = values[0].length;
            }
            out.print("\t");
            out.print(columns);
            for (double[] row : values) {
                if (row == null || row.length != columns) {
                    throw new IllegalArgumentException("Cannot write ragged checkpoint matrix");
                }
                for (double value : row) {
                    out.print("\t");
                    out.print(value);
                }
            }
        }
    }

    class CheckpointStateReader {
        public interface DoubleParser {
            double parseDouble(String string);
        }

        private final String[] fields;
        private final DoubleParser parser;
        private final String stateType;
        private final String operatorName;
        private int index;

        public CheckpointStateReader(String[] fields, int index, String stateType, String operatorName,
                                     DoubleParser parser) {
            this.fields = fields;
            this.index = index;
            this.stateType = stateType;
            this.operatorName = operatorName;
            this.parser = parser;
        }

        public int readInt() {
            try {
                return Integer.parseInt(nextField("integer"));
            } catch (NumberFormatException nfe) {
                throw malformed("integer", nfe);
            }
        }

        public double[] readVector() {
            int length = readLength("vector");
            double[] values = new double[length];
            for (int i = 0; i < length; i++) {
                values[i] = readDouble("vector value");
            }
            return values;
        }

        public double[][] readMatrix() {
            int rows = readLength("matrix row count");
            int columns = readLength("matrix column count");
            double[][] values = new double[rows][columns];
            for (int i = 0; i < rows; i++) {
                for (int j = 0; j < columns; j++) {
                    values[i][j] = readDouble("matrix value");
                }
            }
            return values;
        }

        private int readLength(String fieldName) {
            int length;
            try {
                length = Integer.parseInt(nextField(fieldName));
            } catch (NumberFormatException nfe) {
                throw malformed(fieldName, nfe);
            }
            if (length < 0) {
                throw malformed(fieldName);
            }
            return length;
        }

        private double readDouble(String fieldName) {
            try {
                return parser.parseDouble(nextField(fieldName));
            } catch (RuntimeException rte) {
                throw malformed(fieldName, rte);
            }
        }

        private String nextField(String fieldName) {
            if (index >= fields.length) {
                throw malformed(fieldName);
            }
            return fields[index++];
        }

        private RuntimeException malformed(String fieldName) {
            return new RuntimeException("Malformed " + stateType + " operator state for checkpointable operator " +
                    operatorName + ": missing or invalid " + fieldName);
        }

        private RuntimeException malformed(String fieldName, RuntimeException cause) {
            return new RuntimeException("Malformed " + stateType + " operator state for checkpointable operator " +
                    operatorName + ": missing or invalid " + fieldName, cause);
        }
    }
}
