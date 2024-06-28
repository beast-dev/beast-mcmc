/*
 * PoissonBranchLengthLikelihoodParser.java
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

package dr.evomodelxml.bigfasttree.thorney;

import dr.evomodel.bigfasttree.thorney.PoissonBranchLengthLikelihoodDelegate;
import dr.xml.*;

public class PoissonBranchLengthLikelihoodParser extends AbstractXMLObjectParser {

    public static final String POISSON_BRANCH_LENGTH_LIKELIHOOD = "poissonBranchLengthLikelihood";

    public String getParserName() {
        return POISSON_BRANCH_LENGTH_LIKELIHOOD;
    }

    public Object parseXMLObject(XMLObject xo) throws XMLParseException {
        double scale = xo.getAttribute("scale",1.0);
        return new PoissonBranchLengthLikelihoodDelegate(POISSON_BRANCH_LENGTH_LIKELIHOOD,scale);
    }

    //************************************************************************
    // AbstractXMLObjectParser implementation
    //************************************************************************

    public String getParserDescription() {
        return "This element provides the likelihood of observing a branchlength given a mutation rate.";
    }

    public Class getReturnType() {
        return PoissonBranchLengthLikelihoodDelegate.class;
    }

    public static final XMLSyntaxRule[] rules = {
            AttributeRule.newDoubleRule("scale",true,"a scale factor to multiply by the rate such as sequence length. default is 1"),

    };

    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }
}
