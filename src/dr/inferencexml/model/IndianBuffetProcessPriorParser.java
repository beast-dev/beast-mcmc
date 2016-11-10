/*
 * IndianBuffetProcessPriorParser.java
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

import dr.inference.model.AdaptableSizeFastMatrixParameter;
import dr.inference.model.IndianBuffetProcessPrior;
import dr.inference.model.MatrixParameter;
import dr.inference.model.Parameter;
import dr.xml.*;

/**
 * @author Max Tolkoff
 */
public class IndianBuffetProcessPriorParser extends AbstractXMLObjectParser {
    public static final String INDIAN_BUFFET_PROCESS="indianBuffetProcess";
    public static final String BETA="beta";
    public static final String ALPHA="alpha";
    public static final String DATA="data";


    @Override
    public Object parseXMLObject(XMLObject xo) throws XMLParseException {
        Parameter alpha=(Parameter) xo.getChild(ALPHA).getChild(0);
        AdaptableSizeFastMatrixParameter data=(AdaptableSizeFastMatrixParameter) xo.getChild(DATA).getChild(0);
        Parameter beta;
        if(xo.hasChildNamed(BETA))
        {
            beta=(Parameter) xo.getChild(BETA).getChild(0);
        }
        else
        {
            beta=new Parameter.Default(1);
        }
        return new IndianBuffetProcessPrior(alpha, beta, data);
    }

    @Override
    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }

    private final XMLSyntaxRule[] rules = {
            new ElementRule(BETA,
                    new XMLSyntaxRule[]{
                            new ElementRule(Parameter.class)
                    }
            , true),
            new ElementRule(ALPHA,
                    new XMLSyntaxRule[]{
                            new ElementRule(Parameter.class)
                    }
            ),
            new ElementRule(DATA,
                    new XMLSyntaxRule[]{
                            new ElementRule(AdaptableSizeFastMatrixParameter.class)
                    }
            )
    };

    @Override
    public String getParserDescription() {
        return "Indian Buffet Process prior on a binary matrix parameter";
    }

    @Override
    public Class getReturnType() {
        return IndianBuffetProcessPrior.class;
    }

    @Override
    public String getParserName() {
        return INDIAN_BUFFET_PROCESS;
    }
}
