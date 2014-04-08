/*
 * LatentFactorModelParser.java
 *
 * Copyright (c) 2002-2014 Alexei Drummond, Andrew Rambaut and Marc Suchard
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

package dr.evomodelxml.continuous;

import dr.evomodel.continuous.LatentFactorModel;
import dr.evomodel.tree.TreeModel;
import dr.evomodelxml.treelikelihood.TreeTraitParserUtilities;
import dr.inference.model.*;
import dr.math.matrixAlgebra.Matrix;
import dr.math.matrixAlgebra.Vector;
import dr.xml.*;

import java.util.List;

/**
 * @author Max Tolkoff
 * @author Marc Suchard
 */

public class LatentFactorModelParser extends AbstractXMLObjectParser {
    public final static String LATENT_FACTOR_MODEL = "latentFactorModel";
    public final static String NUMBER_OF_FACTORS = "factorNumber";
    public final static String FACTORS = "factors";
    public final static String DATA = "data";
    public final static String LOADINGS = "loadings";
    public static final String ROW_PRECISION = "rowPrecision";
    public static final String COLUMN_PRECISION = "columnPrecision";


    public String getParserName() {
        return LATENT_FACTOR_MODEL;
    }

    public Object parseXMLObject(XMLObject xo) throws XMLParseException {
//        Parameter latent  = null;
        CompoundParameter factors = (CompoundParameter) xo.getChild(FACTORS).getChild(CompoundParameter.class);

        TreeTraitParserUtilities utilities = new TreeTraitParserUtilities();
//        String traitName = TreeTraitParserUtilities.DEFAULT_TRAIT_NAME;
        String traitName = (String) xo.getAttribute(TreeTraitParserUtilities.TRAIT_NAME);

        TreeModel treeModel = (TreeModel) xo.getChild(TreeModel.class);
        System.err.println("TN: " + traitName);

        TreeTraitParserUtilities.TraitsAndMissingIndices returnValue =
                utilities.parseTraitsFromTaxonAttributes(xo, traitName, treeModel, true);
        CompoundParameter dataParameter = returnValue.traitParameter;
        List<Integer> missingIndices = returnValue.missingIndices;
        traitName = returnValue.traitName;
//
//
//
//
//        MatrixParameter data = (MatrixParameter) xo.getChild(DATA).getChild(MatrixParameter.class);



//        int colDim=treeModel.getTaxonCount();
//        int rowDim=dataParameter.getDimension()/treeModel.getTaxonCount();
//        Parameter[] dataTemp=new Parameter[colDim];
//        for(int i=0; i<colDim; i++)
//        {
//            dataTemp[i] = new Parameter.Default(rowDim);
//            for(int j=0; j<rowDim; j++)
//            {
//                dataTemp[i].setParameterValue(j, dataParameter.getParameterValue(i*rowDim+j));
//            }
//
//        }
//        MatrixParameter dataMatrix=new MatrixParameter(null, dataTemp);
//        System.err.print(new Matrix(dataMatrix.getParameterAsMatrix()));
//        System.err.print(dataMatrix.getRowDimension());
        LowerTriangularMatrixParameter loadings = (LowerTriangularMatrixParameter) xo.getChild(LOADINGS).getChild(MatrixParameter.class);
        DiagonalMatrix rowPrecision = (DiagonalMatrix) xo.getChild(ROW_PRECISION).getChild(MatrixParameter.class);
        DiagonalMatrix colPrecision = (DiagonalMatrix) xo.getChild(COLUMN_PRECISION).getChild(MatrixParameter.class);
        int numFactors = xo.getAttribute(NUMBER_OF_FACTORS, 4);

        return new LatentFactorModel(dataParameter, factors, loadings, rowPrecision, colPrecision, numFactors);
    }

    private static final XMLSyntaxRule[] rules = {
            AttributeRule.newIntegerRule(NUMBER_OF_FACTORS),
            new ElementRule(TreeModel.class),
            AttributeRule.newStringRule(TreeTraitParserUtilities.TRAIT_NAME),
            new ElementRule(TreeTraitParserUtilities.TRAIT_PARAMETER, new XMLSyntaxRule[]{
                    new ElementRule(Parameter.class)
            }),
            new ElementRule(FACTORS, new XMLSyntaxRule[]{
                    new ElementRule(CompoundParameter.class),
            }),
            new ElementRule(LOADINGS, new XMLSyntaxRule[]{
                    new ElementRule(LowerTriangularMatrixParameter.class)
            }),
            new ElementRule(ROW_PRECISION, new XMLSyntaxRule[]{
                    new ElementRule(DiagonalMatrix.class)
            }),
            new ElementRule(COLUMN_PRECISION, new XMLSyntaxRule[]{
                    new ElementRule(DiagonalMatrix.class)
            }),
    };

//    <latentFactorModel>
//      <factors>
//         <parameter idref="factors"/>
//      </factors>
//    </latentFactorModel>


    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }

    @Override
    public String getParserDescription() {
        return "Sets up a latent factor model, with starting guesses for the loadings and factor matrices as well as the data for the factor analysis";
    }

    @Override
    public Class getReturnType() {
        return LatentFactorModel.class;
    }


}
