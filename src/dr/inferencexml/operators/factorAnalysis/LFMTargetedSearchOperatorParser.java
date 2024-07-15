/*
 * LFMTargetedSearchOperatorParser.java
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

package dr.inferencexml.operators.factorAnalysis;

import dr.inference.model.MatrixParameterInterface;
import dr.inference.operators.factorAnalysis.LFMTargetedSearchOperator;
import dr.xml.*;

import java.util.*;

public class LFMTargetedSearchOperatorParser extends AbstractXMLObjectParser {
    public static final String LFM_TARGETED_SEARCH_OPERATOR = "LFMTargetedSearchOperator";
    public static final String SPARSE_MATRIX = "sparseMatrix";
    public static final String LOADINGS_MATRIX = "loadingsMatrix";
    public static final String FACTORS_MATRIX = "factorsMatrix";
    public static final String SPARSE_TARGET_MATRICES = "sparseTargetMatrices";
    public static final String LOADINGS_TARGET_MATRICES = "loadingsTargetMatrices";
    public static final String FACTORS_TARGET_MATRICES = "factorsTargetMatrices";
    public static final String CUTOFFS = "cutoffs";
    public static final String WEIGHT = "weight";


    @Override
    public Object parseXMLObject(XMLObject xo) throws XMLParseException {
        MatrixParameterInterface sparseMatrix = (MatrixParameterInterface) xo.getChild(SPARSE_MATRIX).getChild(MatrixParameterInterface.class);
        MatrixParameterInterface loadingsMatrix = (MatrixParameterInterface) xo.getChild(LOADINGS_MATRIX).getChild(MatrixParameterInterface.class);
        MatrixParameterInterface factorsMatrix = (MatrixParameterInterface) xo.getChild(FACTORS_MATRIX).getChild(MatrixParameterInterface.class);
        MatrixParameterInterface cutoffs = (MatrixParameterInterface) xo.getChild(CUTOFFS).getChild(MatrixParameterInterface.class);

        XMLObject cxo = xo.getChild(SPARSE_TARGET_MATRICES);
        ArrayList<MatrixParameterInterface> sparseTargetList = new ArrayList<MatrixParameterInterface>();
        for (int i = 0; i < cxo.getChildCount(); i++) {
            sparseTargetList.add((MatrixParameterInterface) cxo.getChild(i));
        }
        XMLObject dxo = xo.getChild(LOADINGS_TARGET_MATRICES);
        ArrayList<MatrixParameterInterface> loadingsTargetList = new ArrayList<MatrixParameterInterface>();
        for (int i = 0; i < dxo.getChildCount(); i++) {
            loadingsTargetList.add((MatrixParameterInterface) dxo.getChild(i));
        }
        XMLObject exo = xo.getChild(FACTORS_TARGET_MATRICES);
        ArrayList<MatrixParameterInterface> factorsTargetList = new ArrayList<MatrixParameterInterface>();
        for (int i = 0; i < exo.getChildCount(); i++) {
            factorsTargetList.add((MatrixParameterInterface) exo.getChild(i));
        }
        double weight = xo.getDoubleAttribute(WEIGHT);

        return new LFMTargetedSearchOperator(weight, sparseMatrix, sparseTargetList, factorsMatrix, factorsTargetList, loadingsMatrix, loadingsTargetList, cutoffs);
    }

    @Override
    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }
    XMLSyntaxRule[] rules = new XMLSyntaxRule[]{
            new ElementRule(SPARSE_MATRIX, new XMLSyntaxRule[]{new ElementRule(MatrixParameterInterface.class)}),
            new ElementRule(SPARSE_TARGET_MATRICES, new XMLSyntaxRule[]{new ElementRule(MatrixParameterInterface.class, 1, Integer.MAX_VALUE)}),
            new ElementRule(FACTORS_MATRIX, new XMLSyntaxRule[]{new ElementRule(MatrixParameterInterface.class)}),
            new ElementRule(FACTORS_TARGET_MATRICES, new XMLSyntaxRule[]{new ElementRule(MatrixParameterInterface.class, 1, Integer.MAX_VALUE)}),
            new ElementRule(LOADINGS_MATRIX, new XMLSyntaxRule[]{new ElementRule(MatrixParameterInterface.class)}),
            new ElementRule(LOADINGS_TARGET_MATRICES, new XMLSyntaxRule[]{new ElementRule(MatrixParameterInterface.class, 1, Integer.MAX_VALUE)}),
            new ElementRule(CUTOFFS, new XMLSyntaxRule[]{new ElementRule(MatrixParameterInterface.class)}),

    };

    @Override
    public String getParserDescription() {
        return "Targeted search for Sparse Latent Factor Model";
    }

    @Override
    public Class getReturnType() {
        return LFMTargetedSearchOperator.class;
    }

    @Override
    public String getParserName() {
        return LFM_TARGETED_SEARCH_OPERATOR;
    }
}
