/*
 * SpeciesTreeSimplePriorParser.java
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

package dr.evomodelxml.speciation;

import dr.evomodel.speciation.SpeciesTreeModel;
import dr.evomodel.speciation.SpeciesTreeSimplePrior;
import dr.inference.distribution.ParametricDistributionModel;
import dr.inference.model.Parameter;
import dr.xml.*;

/**
 */
public class SpeciesTreeSimplePriorParser extends AbstractXMLObjectParser {
    private static final String STPRIOR = "speciesTreePopulationPrior";
    public static final String TIPS = "tipsDistribution";

        public String getParserDescription() {
            return "";
        }

        public Class getReturnType() {
            return SpeciesTreeSimplePrior.class;
        }

        public String getParserName() {
            return STPRIOR;
        }

        public Object parseXMLObject(XMLObject xo) throws XMLParseException {
            SpeciesTreeModel st = (SpeciesTreeModel) xo.getChild(SpeciesTreeModel.class);

            //ParametricDistributionModel pr = (ParametricDistributionModel) xo.getChild(ParametricDistributionModel.class);
            Parameter pr = (Parameter)((XMLObject)xo.getChild("sigma")).getChild(Parameter.class);

            final XMLObject cxo = xo.getChild(TIPS);
            final ParametricDistributionModel tipsPrior = (ParametricDistributionModel) cxo.getChild(ParametricDistributionModel.class);

            return new SpeciesTreeSimplePrior(st, pr, tipsPrior);
        }

        public XMLSyntaxRule[] getSyntaxRules() {
            return new XMLSyntaxRule[]{
                    new ElementRule(SpeciesTreeModel.class),
                    new ElementRule(TIPS,
                            new XMLSyntaxRule[] { new ElementRule(ParametricDistributionModel.class) }),
                    //new ElementRule(ParametricDistributionModel.class),
                    new ElementRule("sigma", new XMLSyntaxRule[] { new ElementRule(Parameter.class) })
            };
        }
}
