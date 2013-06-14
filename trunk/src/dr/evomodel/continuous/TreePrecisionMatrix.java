/*
 * TreePrecisionMatrix.java
 *
 * Copyright (c) 2002-2013 Alexei Drummond, Andrew Rambaut and Marc Suchard
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

import dr.inference.model.MatrixParameter;
import dr.inference.model.Model;
import dr.inference.model.ModelListener;
import dr.math.matrixAlgebra.Matrix;
import dr.xml.*;

/**
 * @author Marc Suchard
 */
public class TreePrecisionMatrix extends MatrixParameter implements ModelListener {

    private final FullyConjugateMultivariateTraitLikelihood traitModel;
    private final boolean conditionOnRoot;
    private double[][] precisionMatrix = null;
    private final int dim;

    public TreePrecisionMatrix(FullyConjugateMultivariateTraitLikelihood traitModel, boolean conditionOnRoot) {
        super(MATRIX_PARAMETER);

        traitModel.addModelListener(this);
        this.traitModel = traitModel;

        this.conditionOnRoot = conditionOnRoot;
        dim = traitModel.getTreeModel().getExternalNodeCount() * traitModel.getDimTrait();
    }

    public String toString() {
        if (precisionMatrix == null) {
            getPrecisionMatrix();
        }
        StringBuilder buffer = new StringBuilder();
        buffer.append(new Matrix(precisionMatrix));
        return buffer.toString();
    }

    public double getParameterValue(int index) {
        return getParameterValue(index / getColumnDimension(), index % getColumnDimension());
    }

    public double[] getAttributeValue() {
        double[] stats = new double[dim * dim];
        int index = 0;
        for (int i = 0; i < dim; i++) {
            for (int j = 0; j < dim; j++) {
                stats[index] = getParameterValue(i, j);
                index++;
            }
        }
        return stats;
    }

    private void getPrecisionMatrix() {
        precisionMatrix = MultivariateTraitUtils.computeTreeTraitPrecision(traitModel, conditionOnRoot);
    }

    public double getParameterValue(int row, int col) {
        if (precisionMatrix == null) {
            getPrecisionMatrix();
        }
        return precisionMatrix[row][col];
    }

    public double[][] getParameterAsMatrix() {
        if (precisionMatrix == null) {
            getPrecisionMatrix();
        }
        return precisionMatrix;
    }

    public int getColumnDimension() {
        return dim;
    }

    public int getRowDimension() {
        return dim;
    }

    public void modelChangedEvent(Model model, Object object, int index) {
        precisionMatrix = null;
    }

    public void modelRestored(Model model) {
        precisionMatrix = null;
    }

    public static final String PARSE_OBJECT_NAME = "treeTraitPrecisionMatrix";
    public static final String CONDITION = "conditionOnRoot";


    public static dr.xml.XMLObjectParser PARSER = new dr.xml.AbstractXMLObjectParser() {

        public String getParserName() {
            return PARSE_OBJECT_NAME;
        }

        public Object parseXMLObject(XMLObject xo) throws XMLParseException {

            FullyConjugateMultivariateTraitLikelihood traitModel = (FullyConjugateMultivariateTraitLikelihood)
                    xo.getChild(FullyConjugateMultivariateTraitLikelihood.class);

            boolean condition = xo.getAttribute(CONDITION, false);

            return new TreePrecisionMatrix(traitModel, condition);
        }

        //************************************************************************
        // AbstractXMLObjectParser implementation
        //************************************************************************

        public XMLSyntaxRule[] getSyntaxRules() {
            return rules;
        }

        private final XMLSyntaxRule[] rules = {
                new ElementRule(FullyConjugateMultivariateTraitLikelihood.class),
                AttributeRule.newBooleanRule(CONDITION, true),
        };

        public String getParserDescription() {
            return "This element returns the precision matrix of a multivariate diffusion along a tree.";
        }

        public Class getReturnType() {
            return TreePrecisionMatrix.class;
        }
    };

}
