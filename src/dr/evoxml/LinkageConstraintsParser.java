/*
 * LinkageConstraintsParser.java
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

import dr.evolution.LinkageConstraints;
import dr.evolution.LinkedGroup;
import dr.xml.*;

import java.util.ArrayList;

/**
 * @author Aaron Darling (koadman)
 */
public class LinkageConstraintsParser extends AbstractXMLObjectParser {

	
	public String getParserDescription() {
		return "Data representing metagenome reads that are linked by mate-pair, strobe, or other information";
	}


	public Class getReturnType() {
		return LinkageConstraints.class;
	}


	public XMLSyntaxRule[] getSyntaxRules() {
		return rules;
	}


	public Object parseXMLObject(XMLObject xo) throws XMLParseException {
    	ArrayList<LinkedGroup> groups = new ArrayList<LinkedGroup>();
        for (int i = 0; i < xo.getChildCount(); i++) {
            Object child = xo.getChild(i);
            if (child instanceof LinkedGroup) {
            	groups.add((LinkedGroup)child);
            }
    	}
    	LinkageConstraints lc = new LinkageConstraints(groups);
		return lc;
	}


	public String getParserName() {
		return "LinkageConstraints";
	}

	private XMLSyntaxRule[] rules = new XMLSyntaxRule[]{
			new ElementRule(LinkedGroup.class, 1, Integer.MAX_VALUE),
    };

}
