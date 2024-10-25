/*
 * ScaleData.java
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

package dr.inferencexml.model;

import dr.inference.model.CompoundParameter;
import dr.inference.model.MatrixParameter;
import dr.inference.model.MatrixParameterInterface;
import dr.inference.model.Parameter;
import dr.xml.*;

/**
 * Created by maxryandolinskytolkoff on 12/8/16.
 */
public class ScaleData extends AbstractXMLObjectParser{

    public final static String SCALE_DATA = "scaleData";
    public final static String CONTINUOUS = "continuous";

    @Override
    public String getParserName() {
        return SCALE_DATA;
    }

    @Override
    public Object parseXMLObject(XMLObject xo) throws XMLParseException {
        MatrixParameterInterface scaling = null;
        CompoundParameter temp = null;
        if(xo.getChild(MatrixParameterInterface.class)!=null){
            scaling = (MatrixParameterInterface) xo.getChild(MatrixParameterInterface.class);
        }
        else{
            temp = (CompoundParameter) xo.getChild(CompoundParameter.class);
            scaling = MatrixParameter.recast(temp.getParameterName(), temp);
        }
        Parameter continuous;
        if(xo.hasChildNamed(CONTINUOUS))
            continuous = (Parameter) xo.getChild(CONTINUOUS).getChild(Parameter.class);
        else
            continuous = new Parameter.Default(scaling.getColumnDimension(), 1);

        System.out.println(continuous.getParameterValue(11));

        double[][] aData = scaling.getParameterAsMatrix();
        double[] meanList = new double[scaling.getRowDimension()];
        double[] varList = new double[scaling.getRowDimension()];
        double[] count = new double[scaling.getRowDimension()];
        for (int i = 0; i < scaling.getColumnDimension(); i++) {
            for (int j = 0; j < scaling.getRowDimension(); j++) {
                if (scaling.getParameterValue(j, i) != 0) {
                    meanList[j] += scaling.getParameterValue(j, i);
                    count[j]++;
                }
            }
        }
        for (int i = 0; i < scaling.getRowDimension(); i++) {
            if (continuous.getParameterValue(i) == 1)
                meanList[i] = meanList[i] / count[i];
             else
                meanList[i] = 0;
        }

        double[][] answerTemp = new double[scaling.getRowDimension()][scaling.getColumnDimension()];
        for (int i = 0; i < scaling.getColumnDimension(); i++) {
            for (int j = 0; j < scaling.getRowDimension(); j++) {
                if (aData[j][i] != 0) {
                    answerTemp[j][i] = aData[j][i] - meanList[j];
                }
            }
        }
//        System.out.println(new Matrix(answerTemp));

        for (int i = 0; i < scaling.getColumnDimension(); i++) {
            for (int j = 0; j < scaling.getRowDimension(); j++) {
                varList[j] += answerTemp[j][i] * answerTemp[j][i];
            }
        }

        for (int i = 0; i < scaling.getRowDimension(); i++) {
            if (continuous.getParameterValue(i) == 1) {
                varList[i] = varList[i] / (count[i] - 1);
                varList[i] = StrictMath.sqrt(varList[i]);
            } else {
                varList[i] = 1;
            }
        }

//        System.out.println(data.getColumnDimension());
//        System.out.println(data.getRowDimension());

        for (int i = 0; i < scaling.getColumnDimension(); i++) {
            for (int j = 0; j < scaling.getRowDimension(); j++) {
                scaling.setParameterValue(j, i, answerTemp[j][i] / varList[j]);
            }
        }





        return null;
    }

    @Override
    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }
    private XMLSyntaxRule[] rules = new XMLSyntaxRule[]{
            new OrRule(new ElementRule(MatrixParameterInterface.class), new ElementRule(CompoundParameter.class)),
     new ElementRule(CONTINUOUS, new XMLSyntaxRule[]{
            new ElementRule(Parameter.class)
    }, true)};

    @Override
    public String getParserDescription() {
        return "Standardizes the rows of a data matrix";
    }

    @Override
    public Class getReturnType() {
        return null;
    }
}
