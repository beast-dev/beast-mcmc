/*
 * MultiSpeciesCoalescentParser.java
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

import dr.evomodel.speciation.MultiSpeciesCoalescent;
import dr.evomodel.speciation.SpeciesBindings;
import dr.evomodel.speciation.SpeciesTreeModel;
import dr.xml.*;

/**
 */
public class MultiSpeciesCoalescentParser extends AbstractXMLObjectParser {
    public static final String SPECIES_COALESCENT = "speciesCoalescent";

    public Object parseXMLObject(XMLObject xo) throws XMLParseException {
        final SpeciesBindings sb = (SpeciesBindings) xo.getChild(SpeciesBindings.class);
        final SpeciesTreeModel tree = (SpeciesTreeModel) xo.getChild(SpeciesTreeModel.class);
        return new MultiSpeciesCoalescent(sb, tree);
    }

    public XMLSyntaxRule[] getSyntaxRules() {
        return new XMLSyntaxRule[]{
                new ElementRule(SpeciesBindings.class),
                new ElementRule(SpeciesTreeModel.class),
        };
    }

    public String getParserDescription() {
        return "Compute coalecent log-liklihood of a set of gene trees embedded inside one species tree.";
    }

    public Class getReturnType() {
        return MultiSpeciesCoalescent.class;
    }

    public String getParserName() {
        return SPECIES_COALESCENT;
    }
}
