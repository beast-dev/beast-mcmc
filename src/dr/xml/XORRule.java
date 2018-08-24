/*
 * XORRule.java
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
import java.util.TreeSet;

public class XORRule implements XMLSyntaxRule {

	public XORRule(XMLSyntaxRule a, XMLSyntaxRule b) {
		this(a, b, false);
	}

	public XORRule(XMLSyntaxRule[] rules) {
        this(rules, false);
	}

    public XORRule(XMLSyntaxRule a, XMLSyntaxRule b, boolean optional) {
        this(new XMLSyntaxRule[] {a, b}, optional);
    }

    public XORRule(XMLSyntaxRule[] rules, boolean optional) {
		this.rules = rules;
        this.optional = optional;
	}

	public XMLSyntaxRule[] getRules() { return rules; }

	/**
	 * Returns true if the rule is satisfied for the given XML object.
	 */
	public boolean isSatisfied(XMLObject xo) {

		int satisfiedCount = 0;
        for (XMLSyntaxRule rule : rules) {
            if (rule.isSatisfied(xo)) {
                satisfiedCount += 1;
                if (satisfiedCount > 1) return false;
            }
        }
        return optional || satisfiedCount == 1;
	}

    public boolean containsAttribute(String name) {
        for( XMLSyntaxRule rule : rules ) {
            if( rule.containsAttribute(name) ) {
                return true;
            }
        }
        return false;
    }

    /**
	 * Describes the rule.
	 */
	public String ruleString() {
		StringBuffer buffer = new StringBuffer();
        if (optional) {
            buffer.append("*Optionally, one of \n");
        } else {
            buffer.append("*One of \n");
        }
        for (XMLSyntaxRule rule : rules) {
            buffer.append("*").append(rule.ruleString()).append("\n");
        }
        return buffer.toString();
	}

	/**
	 * Describes the rule.
	 */
	public String htmlRuleString(XMLDocumentationHandler handler) {
		StringBuffer buffer = new StringBuffer("<div class=\"requiredcompoundrule\">One of:\n");
        for (XMLSyntaxRule rule : rules) {
            buffer.append(rule.htmlRuleString(handler));
        }
        buffer.append("</div>\n");
		return buffer.toString();
	}

	/**
	 * Describes the rule.
	 */
	public String wikiRuleString(XMLDocumentationHandler handler, String prefix) {
		StringBuffer buffer = new StringBuffer(prefix + "One of:\n");
        for (XMLSyntaxRule rule : rules) {
            buffer.append(rule.wikiRuleString(handler, prefix + "*"));
        }
        buffer.append("\n");
		return buffer.toString();
	}

    /**
     * Describes the rule.
     */
    public String markdownRuleString(XMLDocumentationHandler handler, String prefix) {
        StringBuffer buffer = new StringBuffer(prefix + "One of:\n");
        for (XMLSyntaxRule rule : rules) {
            buffer.append(rule.markdownRuleString(handler, prefix + "    "));
        }
        buffer.append("\n");
        return buffer.toString();
    }


    /**
	 * Describes the rule.
	 */
	public String ruleString(XMLObject xo) {
		return ruleString();
	}

	/**
	 * @return a set containing the required types of this rule.
	 */
	public Set<Class> getRequiredTypes() {

		Set<Class> set = new TreeSet<Class>(ClassComparator.INSTANCE);

        for (XMLSyntaxRule rule : rules) {
            set.addAll(rule.getRequiredTypes());
        }
        return set;
	}

    public boolean isLegalElementName(String elementName) {
        for (XMLSyntaxRule rule : rules) {
            if( rule.isLegalElementName(elementName) ) {
                return true;
            }
        }
        return false;
    }

    public boolean isLegalElementClass(Class c) {
        for (XMLSyntaxRule rule : rules) {
            if( rule.isLegalElementClass(c) ) {
                return true;
            }
        }
        return false;
    }


    public boolean isLegalSubelementName(String elementName) {
        for (XMLSyntaxRule rule : rules) {
            if( rule.isLegalSubelementName(elementName) ) {
                return true;
            }
        }
        return false;
    }

    private final XMLSyntaxRule[] rules;
    private final boolean optional;
}
