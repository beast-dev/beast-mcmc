/*
 * CompoundLikelihoodParser.java
 *
 * Copyright (C) 2002-2012 Alexei Drummond, Andrew Rambaut & Marc A. Suchard
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

package dr.app.beagle.evomodel.parsers;

import dr.inference.model.*;
import dr.xml.*;

import java.util.logging.Logger;

/**
 * @author Alexei Drummond
 * @author Andrew Rambaut
 * @author Marc Suchard
 */
public class CompoundLikelihoodParser extends AbstractXMLObjectParser {

    public static final String LIKELIHOOD = "likelihood";


    public String getParserName() {
        return LIKELIHOOD;
    }

    public Object parseXMLObject(XMLObject xo) throws XMLParseException {
   
        ThreadedCompoundLikelihood compoundLikelihood = new ThreadedCompoundLikelihood();

        for (int i = 0; i < xo.getChildCount(); i++) {
            if (xo.getChild(i) instanceof Likelihood) {
                compoundLikelihood.addLikelihood((Likelihood) xo.getChild(i));
            } else {

                Object rogueElement = xo.getChild(i);

                throw new XMLParseException("An element (" + rogueElement + ") which is not a likelihood has been added to the "
                        + LIKELIHOOD + " element");
            }
        }

        Logger.getLogger("dr.evomodel").info("Multithreaded Likelihood, using "
                + compoundLikelihood.getLikelihoodCount() + " threads.");

        return compoundLikelihood;
    }

    //************************************************************************
    // AbstractXMLObjectParser implementation
    //************************************************************************

    public String getParserDescription() {
        return "A likelihood function which is simply the product of its component likelihood functions.";
    }

    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }

    private final XMLSyntaxRule[] rules = {
            new ElementRule(Likelihood.class, 1, Integer.MAX_VALUE),
    };

    public Class getReturnType() {
        return CompoundLikelihood.class;
    }
}