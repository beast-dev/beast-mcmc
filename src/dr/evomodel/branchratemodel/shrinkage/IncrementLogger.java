/*
 * IncrementLogger.java
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

package dr.evomodel.branchratemodel.shrinkage;

import dr.evomodel.branchratemodel.AutoCorrelatedBranchRatesDistribution;
import dr.inference.loggers.LogColumn;
import dr.inference.loggers.Loggable;
import dr.inference.loggers.NumberColumn;
import dr.xml.*;

// Logs increments of an auto-correlated branch-rate model
public class IncrementLogger implements Loggable {

    public IncrementLogger(AutoCorrelatedBranchRatesDistribution distribution) {
        this.distribution = distribution;
    }


    final private AutoCorrelatedBranchRatesDistribution distribution;

    @Override
    public LogColumn[] getColumns() {

        int size = distribution.getDimension();

        LogColumn[] array = new LogColumn[size];

        for (int i = 0; i < size; ++i) {

            final int index = i;
            array[i] = new NumberColumn("increment." + index) {
                @Override
                public double getDoubleValue() {
                    return distribution.getIncrement(index);
                }
            };
        }

        return array;
    }

    private static final String INCREMENT_PARAMETER = "incrementParameter";
    public static XMLObjectParser PARSER = new AbstractXMLObjectParser() {

        @Override
        public Object parseXMLObject(XMLObject xo) throws XMLParseException {
            AutoCorrelatedBranchRatesDistribution distribution = (AutoCorrelatedBranchRatesDistribution) xo.getChild(AutoCorrelatedBranchRatesDistribution.class);
            return new IncrementLogger(distribution);
        }

        @Override
        public String getParserName() {
            return INCREMENT_PARAMETER;
        }

        @Override
        public XMLSyntaxRule[] getSyntaxRules() {
            return rules;
        }

        XMLSyntaxRule[] rules = new XMLSyntaxRule[]{
                new ElementRule(AutoCorrelatedBranchRatesDistribution.class)
        };

        @Override
        public String getParserDescription() {
            return "Logs the increments of auto-correlated branch rate model";
        }

        @Override
        public Class getReturnType() {
            return IncrementLogger.class;
        }
    };
}
