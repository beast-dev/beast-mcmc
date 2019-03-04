/*
 * StandarizeTraits.java
 *
 * Copyright (c) 2002-2017 Alexei Drummond, Andrew Rambaut and Marc Suchard
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

package dr.evomodel.continuous;

import dr.inference.model.MatrixParameterInterface;
import dr.xml.*;

import java.util.List;
import java.util.logging.Logger;

/**
 * Created by Marc A. Suchard on 2/3/17.
 */
public class StandardizeTraits {

    private static final String STANDARDIZE_TRAITS = "standardizeTraits";

    private final MatrixParameterInterface matrix;
    private final boolean[] missing;

    private final double targetSd;

    public StandardizeTraits(MatrixParameterInterface matrix) {
        this(matrix, null, 1.0);
    }

    public StandardizeTraits(MatrixParameterInterface matrix, List<Integer> missingIndices, double targetSd) {

        this.matrix = matrix;
        this.missing = new boolean[matrix.getDimension()];

        if (missingIndices != null) {
            for (int m : missingIndices) {
                missing[m] = true;
            }
        }

        this.targetSd = targetSd;
    }

    public String doStandardization(boolean byColumn) {

        StringBuilder sb = new StringBuilder();
        sb.append("Trait standardization report:\n");

        final int dim = (byColumn ? matrix.getColumnDimension() : matrix.getRowDimension());

        for (int index = 0; index < dim; ++index) {
            MeanVariance mv = getStatistics(matrix, index, byColumn);

            sb.append("\tBEFORE\n");
            sb.append("\tindex: ").append(index + 1).append("\n");
            sb.append("\tmean : ").append(mv.mean).append("\n");
            sb.append("\tvar  : ").append(mv.variance).append("\n");
            sb.append("\tcnt  : ").append(mv.count).append("\n");

            updateValues(matrix, mv, index, byColumn);

            mv = getStatistics(matrix, index, byColumn);
            sb.append("\tAFTER\n");
            sb.append("\tindex: ").append(index + 1).append("\n");
            sb.append("\tmean : ").append(mv.mean).append("\n");
            sb.append("\tvar  : ").append(mv.variance).append("\n");
            sb.append("\tcnt  : ").append(mv.count).append("\n\n");

        }

        return sb.toString();
    }

//    public void doStandardization() throws Exception {
//        List<Taxon> taxonList = taxa.asList();
//
//
//        List<double[]> values = new ArrayList<double[]>();
//
//        for (Taxon taxon : taxonList) {
//            String attribute = (String) taxon.getAttribute(traitName);
//            if (attribute == null) {
//                throw new Exception("Taxon " + taxon.getId() + " does not contain trait " + traitName);
//            }
//            values.add(convert(attribute));
//        }
//
//        final int dim = values.get(0).length;
//
//        for (int col = 0; col < dim; ++col) {
//            standardize(values, col);
//        }
//
//        int index = 0;
//        for (Taxon taxon : taxonList) {
//            taxon.setAttribute(traitName, convert(values.get(index)));
//            ++index;
//        }
////        System.exit(-1);
//    }
//
//    private double[] convert(String string) {
//        StringTokenizer st = new StringTokenizer(string);
//        double[] values = new double[st.countTokens()];
//
//        for (int i = 0; i < st.countTokens(); ++i) {
//            String str = st.nextToken();
//
//            double value = Double.NaN;
//            if (str.equals("NA") || str.equals("?")) {
//                // Do nothing
//            } else {
//                value = Double.valueOf(str);
//            }
//            values[i] = value;
//        }
//
//        return values;
//    }
//
//    private String convert(double[] values) {
//        StringBuilder sb = new StringBuilder();
//        for (double x : values) {
//            if (!Double.isNaN(x)) {
//                sb.append(x);
//            } else {
//                sb.append("NA");
//            }
//            sb.append(" ");
//        }
//        return sb.toString().trim();
//    }
//
//    private void standardize(List<double[]> values, int column) {
//
//        MeanVariance mv = getStatistics(values, column);
//        final double sd = Math.sqrt(mv.variance);
//
//
//        System.err.println("Col:  " + column);
//        System.err.println("Mean: " + mv.mean);
//        System.err.println("Var:  " + mv.variance);
//        System.err.println("Count:" + mv.count);
//        System.err.println("");
//
//        for (double[] row : values) {
//            double x = row[column];
//            if (!Double.isNaN(x)) {
//                row[column] = (x - mv.mean) / sd;
//            }
//        }
//    }

    private void updateValues(MatrixParameterInterface matrix, final MeanVariance mv, int major, boolean byColumn) {

        final int dim = (byColumn ? matrix.getRowDimension() : matrix.getColumnDimension());
        final int rowDim = matrix.getRowDimension();
        final double sd = Math.sqrt(mv.variance);

        int offset = byColumn ? dim * major : major;

        for (int index = 0; index < dim; ++index) {

            final int row = byColumn ? index : major;
            final int col = byColumn ? major : index;

            double x = matrix.getParameterValue(row, col);
            if (!Double.isNaN(x) && !missing[offset]) {
                x = (x - mv.mean) / sd * targetSd;
                matrix.setParameterValueQuietly(row, col, x);
            }

            offset += byColumn ? 1 : rowDim;

        }

        matrix.fireParameterChangedEvent();
    }


//    private MeanVariance getStatistics(List<double[]> values, int column) {
//
//        double s = 0.0;
//        double ss = 0.0;
//        int c = 0;
//        for (double[] row : values) {
//            double x = row[column];
//            if (!Double.isNaN(x)) {
//                s += x;
//                ss += x * x;
//                ++c;
//            }
//        }
//
//        MeanVariance mv = new MeanVariance();
//        mv.mean = s / c;
//        mv.variance =  ss / c - mv.mean * mv.mean;
//        mv.count = c;
//
//        return mv;
//    }

    private MeanVariance getStatistics(MatrixParameterInterface matrix, int major, boolean byColumn) {

        double s = 0.0;
        double ss = 0.0;
        int c = 0;

        final int dim = (byColumn ? matrix.getRowDimension() : matrix.getColumnDimension());
        final int rowDim = matrix.getRowDimension();

        int offset = byColumn ? dim * major : major;
        for (int index = 0; index < dim; ++index) {
            double x = byColumn ? matrix.getParameterValue(index, major) : matrix.getParameterValue(major, index);
            if (!Double.isNaN(x) && !missing[offset]) {
                s += x;
                ss += x * x;
                ++c;
            }

            offset += byColumn ? 1 : rowDim;
        }

        MeanVariance mv = new MeanVariance();
        mv.mean = s / c;
        mv.variance = ss / c - mv.mean * mv.mean;
        mv.count = c;

        return mv;
    }

    private class MeanVariance {
        double mean;
        double variance;
        int count;
    }

    public static XMLObjectParser PARSER = new AbstractXMLObjectParser() {

        @Override
        public Object parseXMLObject(XMLObject xo) {

            MatrixParameterInterface matrix = (MatrixParameterInterface) xo.getChild(MatrixParameterInterface.class);

            StandardizeTraits st = new StandardizeTraits(matrix);

            String report = st.doStandardization(false);

            Logger.getLogger("dr.evomodel.continuous").info(report);

            return st;
        }

        /**
         * @return an array of syntax rules required by this element.
         * Order is not important.
         */
        @Override
        public XMLSyntaxRule[] getSyntaxRules() {
            return rules;
        }

        @Override
        public String getParserDescription() {
            return null;
        }

        @Override
        public Class getReturnType() {
            return StandardizeTraits.class;
        }

        /**
         * @return Parser name, which is identical to name of xml element parsed by it.
         */
        @Override
        public String getParserName() {
            return STANDARDIZE_TRAITS;
        }

        private final XMLSyntaxRule[] rules = new XMLSyntaxRule[] {
//                AttributeRule.newStringRule(TRAIT_NAME),
//                new ElementRule(Taxa.class),
                new ElementRule(MatrixParameterInterface.class),
        };
    };
}