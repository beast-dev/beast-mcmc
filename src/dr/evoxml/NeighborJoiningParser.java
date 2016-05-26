/*
 * NeighborJoiningParser.java
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

package dr.evoxml;

import dr.evolution.distance.DistanceMatrix;
import dr.evolution.tree.NeighborJoiningTree;
import dr.evolution.tree.Tree;
import dr.xml.*;

/**
 * @author Alexei Drummond
 * @author Andrew Rambaut
 *
 * @version $Id: NeighborJoiningParser.java,v 1.2 2005/05/24 20:25:59 rambaut Exp $
 */
public class NeighborJoiningParser extends AbstractXMLObjectParser {

    //
    // Public stuff
    //
    public final static String NEIGHBOR_JOINING_TREE = "neighborJoiningTree";

    public String getParserName() { return NEIGHBOR_JOINING_TREE; }

    public Object parseXMLObject(XMLObject xo) throws XMLParseException {

        DistanceMatrix distances = (DistanceMatrix)xo.getChild(DistanceMatrix.class);
        return new NeighborJoiningTree(distances);
    }

    public String getParserDescription() {
        return "This element returns a neighbour-joining tree generated from the given distances.";
    }

    public Class getReturnType() { return Tree.class; }

    public XMLSyntaxRule[] getSyntaxRules() { return rules; }

    private XMLSyntaxRule[] rules = new XMLSyntaxRule[] {
        new ElementRule(DistanceMatrix.class)
    };
}