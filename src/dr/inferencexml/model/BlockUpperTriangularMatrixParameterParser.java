/*
 * BlockUpperTriangularMatrixParameterParser.java
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

import dr.inference.model.BlockUpperTriangularMatrixParameter;
import dr.inference.model.Parameter;
import dr.xml.*;

/**
@author Max Tolkoff
 */
public class BlockUpperTriangularMatrixParameterParser extends AbstractXMLObjectParser {
    private static final String BLOCK_UPPER_TRIANGULAR_MATRIX="blockUpperTriangularMatrixParameter";
    private static final String COLUMN_DIMENSION="columnDimension";
    private static final String TRANSPOSE="transpose";
    private static final String DIAGONAL_RESTRICTION="diagonalRestriction";

    @Override
    public Object parseXMLObject(XMLObject xo) throws XMLParseException {

        final String name = xo.hasId() ? xo.getId() : null;
        final boolean transpose= xo.getAttribute(TRANSPOSE, false);
//        int rowDim=xo.getChildCount();
//        int colDim;
        final boolean diagonalRestriction=xo.getAttribute(DIAGONAL_RESTRICTION, false);
        Parameter temp=null;
//        if(xo.hasAttribute(COLUMN_DIMENSION)) {
//            colDim = xo.getAttribute(COLUMN_DIMENSION, 1);
//        }
//        else
//        {
//            temp=(Parameter) xo.getChild(xo.getChildCount()-1);
//            colDim=temp.getDimension();
//        }

        Parameter[] params=new Parameter[xo.getChildCount()];



        for (int i = 0; i < xo.getChildCount(); i++) {
            temp = (Parameter) xo.getChild(i);
            params[i]=temp;}

        BlockUpperTriangularMatrixParameter ltmp=new BlockUpperTriangularMatrixParameter(name, params, diagonalRestriction);
        if(transpose){
            return ltmp.transposeBlock();
        }
        else {
            return ltmp;  //To change body of implemented methods use File | Settings | File Templates.
        }
    }

    @Override
    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }

    private final XMLSyntaxRule[] rules = {
            new ElementRule(Parameter.class, 0, Integer.MAX_VALUE),
            AttributeRule.newBooleanRule(TRANSPOSE, true),
            AttributeRule.newIntegerRule(COLUMN_DIMENSION, true),
            AttributeRule.newBooleanRule(DIAGONAL_RESTRICTION, true),
    };


    @Override
    public String getParserDescription() {
        return "Returns a blockUpperTriangularMatrixParameter which is a compoundParameter which forces the last element to be of full length, the second to last element to be of full length-1, etc.";  //To change body of implemented methods use File | Settings | File Templates.
    }

    @Override
    public Class getReturnType() {
        return BlockUpperTriangularMatrixParameter.class;  //To change body of implemented methods use File | Settings | File Templates.
    }

    @Override
    public String getParserName() {
        return BLOCK_UPPER_TRIANGULAR_MATRIX;  //To change body of implemented methods use File | Settings | File Templates.
    }
}
