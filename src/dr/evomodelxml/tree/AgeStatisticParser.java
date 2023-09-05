/*
 * TMRCAStatisticParser.java
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

import dr.evolution.tree.Tree;
import dr.evolution.tree.TreeUtils;
import dr.evolution.util.Taxa;
import dr.evolution.util.TaxonList;
import dr.evomodel.tree.AgeStatistic;
import dr.evomodel.tree.TMRCAStatistic;
import dr.inference.model.Statistic;
import dr.xml.*;

/**
 *
 * Converts a statistic or parameter into absolute age:
 * <ageStatistic id="age">
 *     <parameter idref="height"/>
 * </ageStatistic>
 *
 * @author Andrew Rambaut
 */
public class AgeStatisticParser extends AbstractXMLObjectParser {

    public static final String AGE_STATISTIC = "ageStatistic";


    public String getParserName() {
        return AGE_STATISTIC;
    }

    public Object parseXMLObject(XMLObject xo) throws XMLParseException {

        String name = xo.getAttribute(Statistic.NAME, xo.getId());
        Statistic heightStatistic = (Statistic) xo.getChild(Statistic.class);
        return new AgeStatistic(name, heightStatistic);
    }

    //************************************************************************
    // AbstractXMLObjectParser implementation
    //************************************************************************

    public String getParserDescription() {
        return "A statistic that converts a height statistic into an absolute age using the date of the most recent tip. ";
    }

    public Class getReturnType() {
        return AgeStatistic.class;
    }

    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }

    private final XMLSyntaxRule[] rules = {
            new ElementRule(Statistic.class)
    };

}
