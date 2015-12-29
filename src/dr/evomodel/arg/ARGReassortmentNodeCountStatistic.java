/*
 * ARGReassortmentNodeCountStatistic.java
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

package dr.evomodel.arg;

import dr.inference.model.Statistic;
import dr.xml.*;


/**
 * @author Marc Suchard
 */
public class ARGReassortmentNodeCountStatistic extends Statistic.Abstract {

    public static final String REASSORTMENT_STATISTIC = "argReassortmentNodeCount";


    public ARGReassortmentNodeCountStatistic(String name, ARGModel arg) {
        super(name);
        this.arg = arg;      
    }

    public int getDimension() {
        return 1;
    }

    public double getStatisticValue(int dim) {
        return arg.getReassortmentNodeCount();
    }

     public static XMLObjectParser PARSER = new AbstractXMLObjectParser() {

        public String getParserName() {
            return REASSORTMENT_STATISTIC;
        }

        public Object parseXMLObject(XMLObject xo) throws XMLParseException {

            String name = xo.getAttribute(NAME, xo.getId());
            ARGModel arg = (ARGModel) xo.getChild(ARGModel.class);

            return new ARGReassortmentNodeCountStatistic(name,arg);
        }

        //************************************************************************
        // AbstractXMLObjectParser implementation
        //************************************************************************

        public String getParserDescription() {
            return "A statistic that returns the number of reassortment nodes in an ARG";
        }

        public Class getReturnType() {
            return Statistic.class;
        }

        public XMLSyntaxRule[] getSyntaxRules() {
            return rules;
        }

        private XMLSyntaxRule[] rules = new XMLSyntaxRule[]{
                new StringAttributeRule(NAME, "A name for this statistic for the purpose of logging", true),
                new ElementRule(ARGModel.class),
        };

    };


    private ARGModel arg;

}
