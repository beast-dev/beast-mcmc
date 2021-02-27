/*
 * DataFromTreeTipsParser.java
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

package dr.evomodelxml.continuous;

import dr.evomodelxml.treelikelihood.TreeTraitParserUtilities;
import dr.inference.model.MatrixParameter;
import dr.inference.model.MatrixParameterInterface;
import dr.inference.model.Parameter;
import dr.xml.*;

/**
 * @author Max Tolkoff
 * @author Marc Suchard
 */

public class DataFromTreeTipsParser extends AbstractXMLObjectParser {
    public final static String DATA_FROM_TREE_TIPS = "dataFromTreeTips";
    public final static String DATA = "data";
    public static final String CONTINUOUS = "continuous";


    public String getParserName() {
        return DATA_FROM_TREE_TIPS;
    }

    public Object parseXMLObject(XMLObject xo) throws XMLParseException {

        DataAndMissingFromTreeTipsParser dataAndMissingFromTreeTipsParser = new DataAndMissingFromTreeTipsParser();
        TreeTraitParserUtilities.TraitsAndMissingIndices returnValue =
                (TreeTraitParserUtilities.TraitsAndMissingIndices) dataAndMissingFromTreeTipsParser.parseXMLObject(xo);


        MatrixParameterInterface dataParameter;
        if (returnValue.traitParameter instanceof MatrixParameterInterface) {
            dataParameter = (MatrixParameterInterface) returnValue.traitParameter;
        } else {
            dataParameter = MatrixParameter.recast(returnValue.traitParameter.getId(),
                    returnValue.traitParameter);
        }

        if (xo.hasChildNamed(TreeTraitParserUtilities.MISSING)) {
            Parameter missing = (Parameter) xo.getChild(TreeTraitParserUtilities.MISSING).getChild(Parameter.class);
            missing.setDimension(dataParameter.getDimension());

            boolean[] missingIndicators = returnValue.getMissingIndicators();

            for (int i = 0; i < missing.getDimension(); i++) {
                if (missingIndicators[i]) {
                    missing.setParameterValue(i, 1);
                } else {
                    missing.setParameterValue(i, 0);
                }
            }
        }

        return dataParameter;
    }


    public XMLSyntaxRule[] getSyntaxRules() {
        XMLSyntaxRule[] dataAndMissingRules = DataAndMissingFromTreeTipsParser.rules;
        XMLSyntaxRule[] rules = new XMLSyntaxRule[dataAndMissingRules.length + 1];
        System.arraycopy(dataAndMissingRules, 0, rules, 0, dataAndMissingRules.length);
        rules[dataAndMissingRules.length] =
                new ElementRule(TreeTraitParserUtilities.MISSING, new XMLSyntaxRule[]{
                        new ElementRule(Parameter.class)
                }, true);
        return rules;
    }

    @Override
    public String getParserDescription() {
        return "Takes the data from the tips of a tree and puts it into a MatrixParameter";
    }

    @Override
    public Class getReturnType() {
        return MatrixParameter.class;
    }


}
