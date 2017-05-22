
/*
 * MulMSCoalescentParser.java
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

import dr.evomodel.alloppnet.speciation.MulMSCoalescent;
import dr.evomodel.alloppnet.speciation.MulSpeciesBindings;
import dr.evomodel.alloppnet.speciation.MulSpeciesTreeModel;
import dr.xml.*;

/**
 * 
 * @author Graham Jones
 *         Date: 20/12/2011
 */
public class MulMSCoalescentParser extends AbstractXMLObjectParser {
    public static final String MUL_MS_COALESCENT = "mulMSCoalescent";

    public Object parseXMLObject(XMLObject xo) throws XMLParseException {
        final MulSpeciesBindings sb = (MulSpeciesBindings) xo.getChild(MulSpeciesBindings.class);
        final MulSpeciesTreeModel tree = (MulSpeciesTreeModel) xo.getChild(MulSpeciesTreeModel.class);
        return new MulMSCoalescent(sb, tree);
    }

    public XMLSyntaxRule[] getSyntaxRules() {
        return new XMLSyntaxRule[]{
                new ElementRule(MulSpeciesBindings.class),
                new ElementRule(MulSpeciesTreeModel.class),
        };
    }

    public String getParserDescription() {
        return "Compute coalecent log-liklihood of a set of gene trees embedded inside one species tree.";
    }

    public Class getReturnType() {
        return MulMSCoalescent.class;
    }

    public String getParserName() {
        return MUL_MS_COALESCENT;
    }
}
