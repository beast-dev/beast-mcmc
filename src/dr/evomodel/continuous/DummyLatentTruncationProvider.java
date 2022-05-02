/*
 * DummyLatentTrunctationProvider.java
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


package dr.evomodel.continuous;

/**
 * DummyLatentTruncationProvider - Places un-truncated model into context of latent truncation model.
 *
 * @author Gabriel Hassler
 */


import dr.evomodel.treedatalikelihood.TreeDataLikelihood;
import dr.math.distributions.Distribution;
import dr.xml.ElementRule;
import dr.xml.XMLObject;
import dr.xml.XMLParseException;
import dr.xml.XMLSyntaxRule;

public class DummyLatentTruncationProvider implements LatentTruncation {


    final TreeDataLikelihood likelihood;

    DummyLatentTruncationProvider(TreeDataLikelihood likelihood){

        this.likelihood = likelihood;
    }

    @Override
    public boolean validTraitForTip(int tip) {
        return true;
    }

    @Override
    public double getNormalizationConstant(Distribution working) {
        return 1;
    }

    @Override
    public double getLogLikelihood() {
        return likelihood.getLogLikelihood();
    }

    public static dr.xml.AbstractXMLObjectParser PARSER = new dr.xml.AbstractXMLObjectParser(){

        private final String PARSER_NAME = "dummyLatentTruncationProvider";

        @Override
        public Object parseXMLObject(XMLObject xo) throws XMLParseException {

            TreeDataLikelihood treeDataLikelihood = (TreeDataLikelihood) xo.getChild(TreeDataLikelihood.class);
            return new DummyLatentTruncationProvider(treeDataLikelihood);
        }

        private final XMLSyntaxRule[] rules = {
                new ElementRule(TreeDataLikelihood.class)
        };
        @Override
        public XMLSyntaxRule[] getSyntaxRules() {
            return rules;
        }

        @Override
        public String getParserDescription() {
            return "Provides a dummy truncation model with no truncation.";
        }

        @Override
        public Class getReturnType() {
            return DummyLatentTruncationProvider.class;
        }

        @Override
        public String getParserName() {
            return PARSER_NAME;
        }
    };
}
