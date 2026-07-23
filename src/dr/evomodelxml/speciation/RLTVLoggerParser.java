/*
 * RLTVLoggerParser.java
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

import dr.evomodel.tree.TreeModel;
import dr.evomodel.tree.randomlocalmodel.RLTVLogger;
import dr.evomodel.tree.randomlocalmodel.RandomLocalTreeVariable;
import dr.inference.loggers.TabDelimitedFormatter;
import dr.xml.*;

import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;

/**
 * Parses an element from an DOM document into a SpeciationModel. Recognises YuleModel.
 */
public class RLTVLoggerParser extends AbstractXMLObjectParser {

    private static final String RANDOM_LOCAL_LOGGER = "randomLocalLogger";
    private static final String FILENAME = "fileName";
    private static final String LOG_EVERY = "logEvery";

        public String getParserName() {
            return RANDOM_LOCAL_LOGGER;
        }

        public Object parseXMLObject(XMLObject xo) throws XMLParseException {

            TreeModel treeModel = (TreeModel) xo.getChild(TreeModel.class);
            RandomLocalTreeVariable randomLocal =
                    (RandomLocalTreeVariable) xo.getChild(RandomLocalTreeVariable.class);

            String fileName = xo.getStringAttribute(FILENAME);
            int logEvery = xo.getIntegerAttribute(LOG_EVERY);

            TabDelimitedFormatter formatter = null;
            try {
                formatter = new TabDelimitedFormatter(new PrintWriter(new FileWriter(fileName)));
            } catch (IOException e) {
                e.printStackTrace();
            }

            return new RLTVLogger(formatter, logEvery, treeModel, randomLocal);
        }

        //************************************************************************
        // AbstractXMLObjectParser implementation
        //************************************************************************

        public String getParserDescription() {
            return "A speciation model of a Yule process whose rate can evolve down the tree.";
        }

        public Class getReturnType() {
            return RLTVLogger.class;
        }

        public XMLSyntaxRule[] getSyntaxRules() {
            return rules;
        }

        private XMLSyntaxRule[] rules = new XMLSyntaxRule[]{
                new ElementRule(TreeModel.class),
                new ElementRule(RandomLocalTreeVariable.class),
                AttributeRule.newIntegerRule(LOG_EVERY),
                AttributeRule.newStringRule(FILENAME),
        };
}
