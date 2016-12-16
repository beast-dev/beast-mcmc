/*
 * LatentFactorModelParser.java
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

package dr.inferencexml.model;

import dr.inference.model.*;
import dr.xml.*;

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
    public static final String SCALE_DATA="scaleData";
    public static final String CONTINUOUS="continuous";
    public static final String COMPUTE_RESIDUALS_FOR_DISCRETE="computeResidualsForDiscrete";
    public static final String RECOMPUTE_RESIDUALS="recomputeResiduals";
    public static final String RECOMPUTE_FACTORS="recomputeFactors";
    public static final String RECOMPUTE_LOADINGS="recomputeLoadings";
    public static final String MISSING_INDICATOR = "missingIndicator";


    public String getParserName() {
        return LATENT_FACTOR_MODEL;
    }

    public Object parseXMLObject(XMLObject xo) throws XMLParseException {

        MatrixParameterInterface factors;
        if (xo.getChild(FACTORS).getChild(FastMatrixParameter.class) == null)
        {
            CompoundParameter factorsTemp = (CompoundParameter) xo.getChild(FACTORS).getChild(CompoundParameter.class);
            factors = MatrixParameter.recast(factorsTemp.getParameterName(), factorsTemp);
        }
        else {
            factors = (MatrixParameterInterface) xo.getChild(FACTORS).getChild(MatrixParameterInterface.class);
        }
        MatrixParameterInterface dataParameter = (MatrixParameterInterface) xo.getChild(DATA).getChild(MatrixParameterInterface.class);
        MatrixParameterInterface loadings = (MatrixParameterInterface) xo.getChild(LOADINGS).getChild(MatrixParameterInterface.class);
        Parameter missingIndicator = null;
        if(xo.hasChildNamed(MISSING_INDICATOR))
            missingIndicator = (Parameter) xo.getChild(MISSING_INDICATOR).getChild(Parameter.class);
        DiagonalMatrix rowPrecision = (DiagonalMatrix) xo.getChild(ROW_PRECISION).getChild(MatrixParameter.class);
        DiagonalMatrix colPrecision = (DiagonalMatrix) xo.getChild(COLUMN_PRECISION).getChild(MatrixParameter.class);
        boolean newModel= xo.getAttribute(COMPUTE_RESIDUALS_FOR_DISCRETE, true);
        boolean computeResiduals= xo.getAttribute(RECOMPUTE_RESIDUALS, true);
        boolean computeFactors=xo.getAttribute(RECOMPUTE_FACTORS, true);
        boolean computeLoadings=xo.getAttribute(RECOMPUTE_LOADINGS, true);
        Parameter continuous=null;
        if(xo.getChild(CONTINUOUS)!=null)
            continuous=(Parameter) xo.getChild(CONTINUOUS).getChild(Parameter.class);
        else
            continuous=new Parameter.Default(colPrecision.getRowDimension(), 1.0);
        boolean scaleData=xo.getAttribute(SCALE_DATA, true);
 //       int numFactors = xo.getAttribute(NUMBER_OF_FACTORS, 4);
        Parameter temp=null;
//        for(int i=0; i<loadings.getColumnDimension(); i++)
//        {
//            if(loadings.getParameterValue(i,i)<0)
//            {
//               loadings.setParameterValue(i, i, temp.getParameterValue(i));
//            }
//        }


        return new LatentFactorModel(dataParameter, factors, loadings, rowPrecision, colPrecision, missingIndicator, scaleData, continuous, newModel, computeResiduals, computeFactors, computeLoadings);
    }

    private static final XMLSyntaxRule[] rules = {
            AttributeRule.newIntegerRule(NUMBER_OF_FACTORS),
            AttributeRule.newBooleanRule(SCALE_DATA, true),
            AttributeRule.newBooleanRule(COMPUTE_RESIDUALS_FOR_DISCRETE, true),
            AttributeRule.newBooleanRule(RECOMPUTE_FACTORS, true),
            AttributeRule.newBooleanRule(RECOMPUTE_RESIDUALS, true),
            AttributeRule.newBooleanRule(RECOMPUTE_LOADINGS,true),
            new ElementRule(DATA, new XMLSyntaxRule[]{
                    new ElementRule(MatrixParameterInterface.class),
            }),
            new ElementRule(FACTORS, new XMLSyntaxRule[]{
                    new ElementRule(CompoundParameter.class)
            }),
            new ElementRule(LOADINGS, new XMLSyntaxRule[]{
                    new ElementRule(MatrixParameterInterface.class)
            }),
            new ElementRule(ROW_PRECISION, new XMLSyntaxRule[]{
                    new ElementRule(DiagonalMatrix.class)
            }),
            new ElementRule(COLUMN_PRECISION, new XMLSyntaxRule[]{
                    new ElementRule(DiagonalMatrix.class)
            }),
            new ElementRule(CONTINUOUS, new XMLSyntaxRule[]{
                    new ElementRule(Parameter.class)
            }, true),
            new ElementRule(MISSING_INDICATOR, new XMLSyntaxRule[]{
                    new ElementRule(Parameter.class)
            }, true),
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
