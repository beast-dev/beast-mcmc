/*
 * TKF91LikelihoodParser.java
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

package dr.oldevomodelxml.indel;

import dr.evolution.alignment.Alignment;
import dr.oldevomodel.indel.TKF91Likelihood;
import dr.oldevomodel.indel.TKF91Model;
import dr.oldevomodel.sitemodel.GammaSiteModel;
import dr.evomodel.tree.TreeModel;
import dr.xml.*;

/**
 *
 */
public class TKF91LikelihoodParser extends AbstractXMLObjectParser {

    public static final String TKF91_LIKELIHOOD = "tkf91Likelihood";
    public static final String TKF91_DEATH = "deathRate";
    //public static final String MU = "mutationRate";

    public String getParserName() {
        return TKF91_LIKELIHOOD;
    }

    public Object parseXMLObject(XMLObject xo) throws XMLParseException {

        TreeModel tree = (TreeModel) xo.getChild(TreeModel.class);
        Alignment alignment = (Alignment) xo.getChild(Alignment.class);
        GammaSiteModel siteModel = (GammaSiteModel) xo.getChild(GammaSiteModel.class);
        TKF91Model tkfModel = (TKF91Model) xo.getChild(TKF91Model.class);
        return new TKF91Likelihood(tree, alignment, siteModel, tkfModel);
    }

    //************************************************************************
    // AbstractXMLObjectParser implementation
    //************************************************************************

    public String getParserDescription() {
        return "Returns the total likelihood of a single alignment under the TKF91 model, for a given tree. " +
                "In particular all possible ancestral histories of insertions and deletions leading to the " +
                "alignment of sequences at the tips are taken into account.";
    }

    public Class getReturnType() {
        return TKF91Likelihood.class;
    }

    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }

    private final XMLSyntaxRule[] rules = {
            new ElementRule(TreeModel.class),
            new ElementRule(Alignment.class),
            new ElementRule(GammaSiteModel.class),
            new ElementRule(TKF91Model.class)
    };
}
