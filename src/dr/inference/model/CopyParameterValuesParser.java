/*
 * ResetParameterParser.java
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

package dr.inference.model;

import dr.inference.trace.LogFileTraces;
import dr.inference.trace.TraceException;
import dr.math.MathUtils;
import dr.math.distributions.Distribution;
import dr.xml.*;

import java.io.File;
import java.io.IOException;
import java.util.List;

/**
 * Resets a multi-dimensional continuous parameter.
 * 
 * @author Marc A. Suchard
 *
 */
public class CopyParameterValuesParser extends AbstractXMLObjectParser {

    public static final String RESET_PARAMETER = "copyParameterValues";
    public static final String SOURCE = "source";
    public static final String DESTINATION = "destination";

    public String getParserName() {
        return RESET_PARAMETER;
    }

    public Object parseXMLObject(XMLObject xo) throws XMLParseException {

        Parameter source = (Parameter) xo.getChild(SOURCE).getChild(Parameter.class);
        Parameter destination = (Parameter) xo.getChild(DESTINATION).getChild(Parameter.class);

        if (source.getDimension() != destination.getDimension()) {
            throw new XMLParseException("Source (" + source.getDimension() + ") and destination (" +
                    destination.getDimension() + ") dimensions do not match");
        }

        for (int i = 0; i < source.getDimension(); ++i) {
            destination.setParameterValueQuietly(i, source.getParameterValue(i));
        }
        destination.fireParameterChangedEvent();

        return destination;
    }

    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }

    private final XMLSyntaxRule[] rules = {
            new ElementRule(SOURCE, new XMLSyntaxRule[] {
                    new ElementRule(Parameter.class),
            }),
            new ElementRule(DESTINATION, new XMLSyntaxRule[] {
                    new ElementRule(Parameter.class),
            }),
    };


    public String getParserDescription() {
        return "Copy parameter values from source to destination";
    }

    public Class getReturnType() {
        return Parameter.class;
    }

}
