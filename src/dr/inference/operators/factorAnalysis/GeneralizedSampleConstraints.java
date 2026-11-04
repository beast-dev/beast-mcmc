/*
 * GeneralizedSampleConstraints.java
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

package dr.inference.operators.factorAnalysis;

import dr.xml.*;

import java.util.ArrayList;


/**
 * @author Gabriel Hassler
 * @author Marc A. Suchard
 */

public class GeneralizedSampleConstraints implements LoadingsSamplerConstraints {
    private final ArrayList<int[]> columnIndices;
    private final ArrayList<Integer> uniqueLengths;
    private final int[] arrayIndices;


    GeneralizedSampleConstraints(ArrayList<int[]> columnIndices) {
        this.columnIndices = columnIndices;

        int nCols = columnIndices.size();
        this.uniqueLengths = new ArrayList<>();
        this.arrayIndices = new int[nCols];
        for (int i = 0; i < nCols; i++) {
            int dim = columnIndices.get(i).length;
            int arrayIndex = -1;
            for (int j = 0; j < uniqueLengths.size(); j++) {
                if (uniqueLengths.get(j) == dim) {
                    arrayIndex = j;
                    break;
                }
            }

            if (arrayIndex == -1) {
                arrayIndex = uniqueLengths.size();
                uniqueLengths.add(dim);
            }

            arrayIndices[i] = arrayIndex;
        }
    }


    @Override
    public int getColumnDim(int colIndex, int nRows) {
        return columnIndices.get(colIndex).length; //TODO: return int[] of actual indices, not just assume 0:(n- 1)
    }

    @Override
    public int getArrayIndex(int colIndex, int nRows) {
        return arrayIndices[colIndex];
    }

    @Override
    public void allocateStorage(ArrayList<double[][]> precisionArray, ArrayList<double[]> midMeanArray, ArrayList<double[]> meanArray, int nRows) {
        for (int i = 0; i < uniqueLengths.size(); i++) {
            int dim = uniqueLengths.get(i);
            precisionArray.add(new double[dim][dim]);
            midMeanArray.add(new double[dim]);
            meanArray.add(new double[dim]);
        }
    }

    public static final AbstractXMLObjectParser PARSER = new AbstractXMLObjectParser() {

        private static final String SAMPLE_COLUMNS = "sampleColumns";
        private static final String INDICES = "indices";
        private static final String TRAITS = "traits";
        private static final String ROWS = "rows";

        @Override
        public Object parseXMLObject(XMLObject xo) throws XMLParseException {
            ArrayList<int[]> columnIndices = new ArrayList<>();
            int nextTrait = 1;
            for (int i = 0; i < xo.getChildCount(); i++) {
                XMLObject cxo = (XMLObject) xo.getChild(i);
                int[] traits = cxo.getIntegerArrayAttribute(TRAITS);
                int[] rows = cxo.getIntegerArrayAttribute(ROWS);

                for (int trait : traits) {
                    if (trait != nextTrait) {
                        throw new XMLParseException("Currently only implemented for sequential '" + TRAITS +
                                "' values.");
                    }

                    for (int j = 0; j < rows.length; j++) {
                        rows[j] -= 1;
                    }
                    columnIndices.add(rows);
                    nextTrait++;
                }

            }
            return new GeneralizedSampleConstraints(columnIndices);
        }

        @Override
        public XMLSyntaxRule[] getSyntaxRules() {
            return new XMLSyntaxRule[]{
                    new ElementRule(INDICES, new XMLSyntaxRule[]{
                            AttributeRule.newIntegerArrayRule(TRAITS, false),
                            AttributeRule.newIntegerArrayRule(ROWS, false)
                    }, 1, Integer.MAX_VALUE)
            };
        }

        @Override
        public String getParserDescription() {
            return "Sample from only certain elements of the loadings matrix";
        }

        @Override
        public Class getReturnType() {
            return LoadingsSamplerConstraints.class;
        }

        @Override
        public String getParserName() {
            return SAMPLE_COLUMNS;
        }
    };
}
