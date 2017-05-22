/*
 * MarkovJumpsLikelihoodLoggerParser.java
 *
 * Copyright (c) 2002-2016 Alexei Drummond, Andrew Rambaut and Marc Suchard
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

package dr.evomodelxml.treelikelihood;

import dr.evomodel.treelikelihood.MarkovJumpsTraitProvider;
import dr.inference.loggers.LogColumn;
import dr.inference.loggers.Loggable;
import dr.inference.loggers.NumberColumn;
import dr.xml.*;

/**
 * @author Marc A. Suchard
 */
public class MarkovJumpsLikelihoodLoggerParser extends AbstractXMLObjectParser {

    public static final String PARSER_NAME = "dataLikelihood";

    @Override
    public Object parseXMLObject(XMLObject xo) throws XMLParseException {

        final MarkovJumpsTraitProvider dataLikelihood =
                (MarkovJumpsTraitProvider) xo.getChild(MarkovJumpsTraitProvider.class);

        return new Loggable() {
            public LogColumn[] getColumns() {
                return new LogColumn[]{
                        new LikelihoodColumn(dataLikelihood, "dataLike")
                };
            }
        };
    }

    protected class LikelihoodColumn extends NumberColumn {
        final private MarkovJumpsTraitProvider tree;

        public LikelihoodColumn(MarkovJumpsTraitProvider tree, String label) {
            super(label);
            this.tree = tree;
        }

        public double getDoubleValue() {
            return tree.getLogLikelihood();
        }
    }

    @Override
    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }

    private static XMLSyntaxRule[] rules = {
            new ElementRule(MarkovJumpsTraitProvider.class),
    };

    @Override
    public String getParserDescription() {
        return null;
    }

    @Override
    public Class getReturnType() {
        return Loggable.class;
    }

    public String getParserName() {
        return PARSER_NAME;
    }
}
