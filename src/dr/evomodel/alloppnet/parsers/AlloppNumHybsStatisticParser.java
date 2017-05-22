/*
 * AlloppNumHybsStatisticParser.java
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

import dr.evomodel.alloppnet.speciation.AlloppNumHybsStatistic;
import dr.evomodel.alloppnet.speciation.AlloppSpeciesNetworkModel;

import dr.xml.*;

/**
 * @author Graham Jones
 * Date: 08/10/2012
 */
public class AlloppNumHybsStatisticParser  extends AbstractXMLObjectParser {
    public static final String NUMHYBS_STATISTIC = "alloppNumHybsStatistic";
    public static final String APSPNETWORK = "apspNetwork";

    public String getParserName() {
        return NUMHYBS_STATISTIC;
    }


    @Override
    public Object parseXMLObject(XMLObject xo) throws XMLParseException {
        final XMLObject asnmxo = xo.getChild(APSPNETWORK);
        AlloppSpeciesNetworkModel aspnet = (AlloppSpeciesNetworkModel) asnmxo.getChild(AlloppSpeciesNetworkModel.class);
        return new AlloppNumHybsStatistic(aspnet);
    }



    private  XMLSyntaxRule[] asnmRules() {
        return new XMLSyntaxRule[]{
                new ElementRule(AlloppSpeciesNetworkModel.class)
        };

    }

    @Override
    public XMLSyntaxRule[] getSyntaxRules() {
        return new XMLSyntaxRule[]{
                new ElementRule(APSPNETWORK, asnmRules()),

        };
    }



    @Override
    public String getParserDescription() {
        return "Statistic for number of hybridizations in allopolyploid network";
    }

    @Override
    public Class getReturnType() {
        return AlloppNumHybsStatistic.class;
    }
}
