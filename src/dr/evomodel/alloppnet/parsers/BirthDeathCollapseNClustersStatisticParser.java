/*
 * BirthDeathCollapseNClustersStatisticParser.java
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

import dr.evomodel.alloppnet.speciation.BirthDeathCollapseModel;
import dr.evomodel.alloppnet.speciation.BirthDeathCollapseNClustersStatistic;
import dr.evomodel.speciation.SpeciesTreeModel;
import dr.xml.*;

/**
 * @author Graham Jones
 *         Date: 01/09/2013
 */
public class BirthDeathCollapseNClustersStatisticParser extends AbstractXMLObjectParser {
    public static final String BDC_NCLUSTERS_STATISTIC = "bdcNClustersStatistic";
    public static final String SPECIES_TREE = "speciesTree";
    public static final String COLLAPSE_MODEL = "collapseModel";



    public String getParserName() {
        return BDC_NCLUSTERS_STATISTIC;
    }


    @Override
    public Object parseXMLObject(XMLObject xo) throws XMLParseException {
        xo.getAttribute("name");
        final XMLObject spptreexo = xo.getChild(SPECIES_TREE);
        SpeciesTreeModel spptree = (SpeciesTreeModel) spptreexo.getChild(SpeciesTreeModel.class);
        final XMLObject cmxo = xo.getChild(COLLAPSE_MODEL);
        BirthDeathCollapseModel bdcm = (BirthDeathCollapseModel)cmxo.getChild(BirthDeathCollapseModel.class);
        return new BirthDeathCollapseNClustersStatistic(spptree, bdcm);
    }


    private  XMLSyntaxRule[] spptreeRules() {
        return new XMLSyntaxRule[]{
                new ElementRule(SpeciesTreeModel.class)
        };
    }

    private  XMLSyntaxRule[] bdcmRules() {
        return new XMLSyntaxRule[]{
                new ElementRule(BirthDeathCollapseModel.class)
        };
    }



    @Override
    public XMLSyntaxRule[] getSyntaxRules() {
        return new XMLSyntaxRule[]{
                AttributeRule.newStringRule("name"),
                new ElementRule(SPECIES_TREE, spptreeRules()),
                new ElementRule(COLLAPSE_MODEL, bdcmRules())
        };
    }


    @Override
    public String getParserDescription() {
        return "Statistic for number of collapsed nodes in species tree when using birth-death-collapse model.";
    }

    @Override
    public Class getReturnType() {
        return BirthDeathCollapseNClustersStatistic.class;
    }


}
