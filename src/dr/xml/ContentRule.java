/*
 * ContentRule.java
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

package dr.xml;

import java.util.Collections;
import java.util.Set;

/**
 * A syntax rule to ensure that allows one to document arbitrary content.
 */
public class ContentRule implements XMLSyntaxRule {

	/**
	 * Creates a required element rule.
	 */
	public ContentRule(String htmlDescription) {
		this.htmlDescription = htmlDescription;
	}

	/**
	 * @return true
	 */
	public boolean isSatisfied(XMLObject xo) { return true; }

    public boolean containsAttribute(String name) {
        return false;
    }

    /**
	 * @return a string describing the rule.
	 */
	public String ruleString() { return htmlDescription; }

	/**
	 * @return a string describing the rule.
	 */
	public String htmlRuleString(XMLDocumentationHandler handler) {
		return htmlDescription;
	}

	/**
	 * @return a string describing the rule.
	 */
	public String wikiRuleString(XMLDocumentationHandler handler, String prefix) {
		return prefix + ":" + htmlDescription;
	}

	/**
	 * @return a string describing the rule.
	 */
	public String markdownRuleString(XMLDocumentationHandler handler, String prefix) {
		return prefix + ": " + htmlDescription + "\n";
	}

	/**
	 * @return a string describing the rule.
	 */
	public String ruleString(XMLObject xo) { return null; }

	/**
	 * @return a set containing the required types of this rule.
	 */
	public Set<Class> getRequiredTypes() { return Collections.EMPTY_SET; }

    public boolean isLegalElementName(String elementName) {
        return true;
    }

    public boolean isLegalSubelementName(String elementName) {
        return true;
    }

    public boolean isLegalElementClass(Class c) {
        return true;
    }

    public boolean isAttributeRule() { return false; }

	private final String htmlDescription;
}
