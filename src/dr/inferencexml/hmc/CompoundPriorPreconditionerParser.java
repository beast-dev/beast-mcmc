/*
 * CompoundPriorPreconditionerParser.java
 *
 * Copyright (c) 2002-2017 Alexei Drummond, Andrew Rambaut and Marc Suchard
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

package dr.inferencexml.hmc;

import dr.inference.hmc.CompoundPriorPreconditioner;
import dr.inference.model.PriorPreconditioningProvider;
import dr.xml.*;

import java.util.List;

/**
 * @author Alexander Fisher
 */

public class CompoundPriorPreconditionerParser extends AbstractXMLObjectParser {

    public static final String COMPOUND_PRIOR_PRECONDITIONER = "compoundPriorPreconditioner";

    @Override
    public Object parseXMLObject(XMLObject xo) throws XMLParseException {
        List<PriorPreconditioningProvider> priorList;

        priorList = xo.getAllChildren(PriorPreconditioningProvider.class);

        return new CompoundPriorPreconditioner(priorList);
    }

    @Override
    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }

    private final XMLSyntaxRule[] rules = {
            new ElementRule(PriorPreconditioningProvider.class, 1, Integer.MAX_VALUE),
    };

    @Override
    public String getParserDescription() {
        return null;
    }

    @Override
    public Class getReturnType() {
        return CompoundPriorPreconditioner.class;
    }

    @Override
    public String getParserName() {
        return COMPOUND_PRIOR_PRECONDITIONER;
    }
}
