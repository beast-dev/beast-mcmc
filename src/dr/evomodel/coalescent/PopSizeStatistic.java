/*
 * PopSizeStatistic.java
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

import dr.inference.model.Statistic;
import dr.inference.model.Variable;
import dr.xml.*;

/**
 * A statistic that the population size at a particular time...
 *
 * @author Alexei Drummond
 * @version $Id: PopSizeStatistic.java,v 1.1 2005/10/28 02:49:39 alexei Exp $
 */
public class PopSizeStatistic extends Statistic.Abstract {

    public static final String POPSIZE_STATISTIC = "popSizeStatistic";

    public DemographicModel model;
    public double time;

    public PopSizeStatistic(String name, DemographicModel model, double time) {
        super(name);
        this.model = model;
        this.time = time;
    }

    public int getDimension() {
        return 2 + model.getVariableCount();
    }

    /**
     * @return the height of the MRCA node.
     */
    public double getStatisticValue(int dim) {
        if (dim == 0) return model.getDemographicFunction().getDemographic(time);
        if (dim == 1) return model.getDemographicFunction().getIntensity(time);
        return ((Variable<Double>)model.getVariable(dim - 2)).getValue(0);
    }

    public static XMLObjectParser PARSER = new AbstractXMLObjectParser() {

        public String getParserName() {
            return POPSIZE_STATISTIC;
        }

        public Object parseXMLObject(XMLObject xo) throws XMLParseException {

            String name = xo.getAttribute(NAME, xo.getId());
            DemographicModel demo = (DemographicModel) xo.getChild(DemographicModel.class);
            double time = xo.getDoubleAttribute("time");


            return new PopSizeStatistic(name, demo, time);
        }

        //************************************************************************
        // AbstractXMLObjectParser implementation
        //************************************************************************

        public String getParserDescription() {
            return "A statistic that has as its value the height of the most recent common ancestor of a set of taxa in a given tree";
        }

        public Class getReturnType() {
            return PopSizeStatistic.class;
        }

        public XMLSyntaxRule[] getSyntaxRules() {
            return rules;
        }

        private XMLSyntaxRule[] rules = new XMLSyntaxRule[]{
                new ElementRule(DemographicModel.class),
                new StringAttributeRule("name", "A name for this statistic primarily for the purposes of logging", true),
                AttributeRule.newDoubleRule("time", false),
        };
    };
}
