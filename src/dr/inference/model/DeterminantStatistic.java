/*
 * DeterminantStatistic.java
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

package dr.inference.model;

import dr.xml.*;
import org.ejml.data.DenseMatrix64F;
import org.ejml.ops.CommonOps;

/**
 * @author Gabriel Hassler
 * @author Marc A. Suchard
 */

public class DeterminantStatistic extends Statistic.Abstract implements VariableListener {

    private final MatrixParameterInterface matrix;
    private final int matrixDim;
    private boolean detKnown = false;
    private double det;

    public DeterminantStatistic(String name, MatrixParameterInterface matrix) {
        super(name);

        this.matrix = matrix;
        this.matrixDim = matrix.getRowDimension();
        matrix.addParameterListener(this);

    }

    @Override
    public int getDimension() {
        return 1;
    }

    @Override
    public double getStatisticValue(int dim) {
        if (!detKnown) {
            double[] values = matrix.getParameterValues();
            DenseMatrix64F M = DenseMatrix64F.wrap(matrixDim, matrixDim, values);
            det = CommonOps.det(M);
            detKnown = true;
        }

        return det;
    }

    @Override
    public void variableChangedEvent(Variable variable, int index, Variable.ChangeType type) {
        detKnown = false;
    }

    public static final AbstractXMLObjectParser PARSER = new AbstractXMLObjectParser() {

        private static final String DETERMINANT_STATISTIC = "determinant";


        @Override
        public Object parseXMLObject(XMLObject xo) throws XMLParseException {
            MatrixParameterInterface matrix = (MatrixParameterInterface) xo.getChild(MatrixParameterInterface.class);
            if (matrix.getColumnDimension() != matrix.getRowDimension()) {
                throw new XMLParseException("can only calculate determinant for square matrices");
            }
            final String name;
            if (xo.hasId()) {
                name = xo.getId();
            } else if (matrix.getId() != null) {
                name = "determinant." + matrix.getId();
            } else {
                name = "determinant";
            }
            return new DeterminantStatistic(name, matrix);
        }

        @Override
        public XMLSyntaxRule[] getSyntaxRules() {
            return new XMLSyntaxRule[]{
                    new ElementRule(MatrixParameterInterface.class)
            };
        }

        @Override
        public String getParserDescription() {
            return "Statistic that computes the determinant of a matrix";
        }

        @Override
        public Class getReturnType() {
            return DeterminantStatistic.class;
        }

        @Override
        public String getParserName() {
            return DETERMINANT_STATISTIC;
        }
    };


}
