/*
 * XMLSyntaxRule.java
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

import java.util.Set;

public interface XMLSyntaxRule {

	/**
	 * Returns true if the rule is satisfied for the given XML object.
	 */
	public boolean isSatisfied(XMLObject object);

    /**
     * Check if rule contains attribute of that name
     *
     * @param name attribute name
     * @return   true if contains attribute
     */
    public boolean containsAttribute(String name);

    /**
	 * Describes the rule in general.
	 */
	public String ruleString();

	/**
	 * Describes the rule in html.
	 */
	public String htmlRuleString(XMLDocumentationHandler handler);

	/**
	 * Describes the rule in wiki.
	 */
	public String wikiRuleString(XMLDocumentationHandler handler, String prefix);

	/**
	 * Describes the rule in markdown.
	 */
	public String markdownRuleString(XMLDocumentationHandler handler, String prefix);

	/**
	 * Describes the rule as pertains to the given object.
	 * In particular if object does not satisfy the rule then how.
	 */
	public String ruleString(XMLObject object);

	/**
	 * @return the classes potentially required by this rule.
	 */
	public Set<Class> getRequiredTypes();

    /**
     *  Check for possible elements: catch typos, old syntax and elements with identical names to global
     *  xml element parsers.
     * @param elementName
     * @return true if rule allows a element with that name
     */
    boolean isLegalElementName(String elementName);

    /**
     *
     * @param c  class type
     * @return true if rule accepts an element which, after parsing, is represented as a class of type 'c'
     */
    boolean isLegalElementClass(Class c);

    /**
     *
     * @param elementName
     * @return true if rule allows a sub-element with that name
     */
    boolean isLegalSubelementName(String elementName);
}
