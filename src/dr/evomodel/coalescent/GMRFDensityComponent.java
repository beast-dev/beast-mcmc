/*
 * GMRFDensityComponent.java
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

package dr.evomodel.coalescent;

import dr.inference.model.Likelihood;
import dr.inference.model.Model;
import dr.inference.loggers.LogColumn;
import dr.inference.loggers.NumberColumn;
import dr.xml.*;

import java.util.Arrays;
import java.util.HashSet;
import java.util.Set;

/**
 * @author Marc Suchard
 */
public class GMRFDensityComponent implements Likelihood {

    public static final String PARSER_NAME = "gmrfDensityComponent";
    public static final String COMPONENT = "component";
    public static final String COALESCENT_TAG = "coalescent";
    public static final String GMRF_TAG = "gmrf";

    public GMRFDensityComponent(GMRFSkyrideLikelihood skyride, boolean returnCoalescent) {
        this.skyride = skyride;
        this.returnCoalescent = returnCoalescent;
        if ( returnCoalescent ) {
            tag = "." + COALESCENT_TAG;
        } else {
            tag = "." + GMRF_TAG;
        }
    }

    public Model getModel() {
        return skyride;
    }

    public double getLogLikelihood() {

        double skyrideLogLikelihood = skyride.getLogLikelihood(); // Handles updating flags, etc.

        if ( returnCoalescent ) {
            skyrideLogLikelihood -= skyride.peakLogFieldLikelihood();
        } else {
            skyrideLogLikelihood -= skyride.peakLogCoalescentLikelihood();
        }
        return skyrideLogLikelihood;
    }

    public boolean evaluateEarly() {
        return false;
    }

    public void makeDirty() {
        skyride.makeDirty();
    }

    public String prettyName() {
        return skyride.prettyName() + tag;
    }

    @Override
    public Set<Likelihood> getLikelihoodSet() {
        return new HashSet<Likelihood>(Arrays.asList(this));
    }

    public boolean isUsed() {
        return skyride.isUsed();
    }

    public void setUsed() {
        skyride.setUsed();
    }

     public LogColumn[] getColumns() {
        return new LogColumn[]{
                new LikelihoodColumn(getId())
        };
    }

    protected class LikelihoodColumn extends NumberColumn {
        public LikelihoodColumn(String label) {
            super(label);
        }

        public double getDoubleValue() {
            return getLogLikelihood();
        }
    }

    public String getId() {
        return skyride.getId() + tag;
    }

    public void setId(String id) {
        skyride.setId(id);
    }

    public static XMLObjectParser PARSER = new AbstractXMLObjectParser() {

        public String getParserName() {
            return PARSER_NAME;
        }

        public Object parseXMLObject(XMLObject xo) throws XMLParseException {

            GMRFSkyrideLikelihood skyride = (GMRFSkyrideLikelihood) xo.getChild(GMRFSkyrideLikelihood.class);
            String componentString = (String) xo.getAttribute(COMPONENT);
            boolean returnCoalescent;
            if (componentString.compareToIgnoreCase(COALESCENT_TAG) == 0) {
                returnCoalescent = true;
            } else if (componentString.compareToIgnoreCase(GMRF_TAG) == 0) {
                returnCoalescent = false;
            } else {
                throw new XMLParseException("Unknown component of GMRF Skyride");
            }

            return new GMRFDensityComponent(skyride,returnCoalescent);            
        }

        //************************************************************************
        // AbstractXMLObjectParser implementation
        //************************************************************************

        public String getParserDescription() {
            return "This element represents a single component (coalescent or field) of the GMRF Skyride.";
        }

        public Class getReturnType() {
            return Likelihood.class;
        }

        public XMLSyntaxRule[] getSyntaxRules() {
            return rules;
        }

        private final XMLSyntaxRule[] rules = {
                AttributeRule.newStringRule(COMPONENT),
                new ElementRule(GMRFSkyrideLikelihood.class),
        };
    };


    private GMRFSkyrideLikelihood skyride;
    private boolean returnCoalescent;
    private String tag;
}
