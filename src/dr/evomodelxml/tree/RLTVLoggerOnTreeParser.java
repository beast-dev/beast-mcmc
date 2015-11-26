/*
 * RLTVLoggerOnTreeParser.java
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

package dr.evomodelxml.tree;

import dr.evomodel.tree.randomlocalmodel.RLTVLoggerOnTree;
import dr.evomodel.tree.randomlocalmodel.RandomLocalTreeVariable;
import dr.xml.*;

/**
 * @author Marc A. Suchard
 */
public class RLTVLoggerOnTreeParser extends AbstractXMLObjectParser {

    private static final String RANDOM_LOCAL_LOGGER = "randomLocalLoggerOnTree";

        public String getParserName() {
            return RANDOM_LOCAL_LOGGER;
        }

        public Object parseXMLObject(XMLObject xo) throws XMLParseException {

            RandomLocalTreeVariable treeVariable = (RandomLocalTreeVariable) xo.getChild(RandomLocalTreeVariable.class);
            return new RLTVLoggerOnTree(treeVariable);
        }

        //************************************************************************
        // AbstractXMLObjectParser implementation
        //************************************************************************

        public String getParserDescription() {
            return "A parser to log changed == 0 : 1 in a tree log";
        }

        public Class getReturnType() {
            return RLTVLoggerOnTree.class;
        }

        public XMLSyntaxRule[] getSyntaxRules() {
            return rules;
        }

        private XMLSyntaxRule[] rules = new XMLSyntaxRule[]{
                new ElementRule(RandomLocalTreeVariable.class),
        };
}
