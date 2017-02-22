/*
 * AlloppChangeNumHybridizationsParser.java
 *
 * Copyright (c) 2002-2017 Alexei Drummond, Andrew Rambaut and Marc Suchard
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

package dr.evomodel.alloppnet.parsers;

import dr.evomodel.alloppnet.operators.AlloppChangeNumHybridizations;
import dr.evomodel.alloppnet.speciation.AlloppSpeciesBindings;
import dr.evomodel.alloppnet.speciation.AlloppSpeciesNetworkModel;
import dr.inference.operators.MCMCOperator;
import dr.xml.*;

/**
 * Created with IntelliJ IDEA.
 * User: Graham
 * Date: 22/07/12
 * Time: 11:24
 */
public class AlloppChangeNumHybridizationsParser extends AbstractXMLObjectParser {

    public static final String CHANGE_NUM_HYBRIDIZATIONS = "changeNumHybridizations";


    public String getParserName() {
        return CHANGE_NUM_HYBRIDIZATIONS;
    }

    @Override
    public Object parseXMLObject(XMLObject xo) throws XMLParseException {
        AlloppSpeciesBindings apsp = (AlloppSpeciesBindings) xo.getChild(AlloppSpeciesBindings.class);
        AlloppSpeciesNetworkModel apspnet = (AlloppSpeciesNetworkModel) xo.getChild(AlloppSpeciesNetworkModel.class);
        if (apspnet.getOneHybridization()) {
            throw new XMLParseException("oneHybridization is set to true but there is a changeNumHybridizations operator");
        }

        final double weight = xo.getDoubleAttribute(MCMCOperator.WEIGHT);
        return new AlloppChangeNumHybridizations(apspnet, apsp, weight);
    }

    @Override
    public XMLSyntaxRule[] getSyntaxRules() {
        return new XMLSyntaxRule[]{
                AttributeRule.newDoubleRule(MCMCOperator.WEIGHT),
                new ElementRule(AlloppSpeciesBindings.class),
                new ElementRule(AlloppSpeciesNetworkModel.class)
        };
    }

    @Override
    public String getParserDescription() {
        return "Operator which changes the number of tetraploid subtrees (hybridizations) by merging and splitting them.";

    }

    @Override
    public Class getReturnType() {
        return AlloppChangeNumHybridizations.class;
    }



}
