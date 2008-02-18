/*
 * OrRule.java
 *
 * Copyright (C) 2002-2006 Alexei Drummond and Andrew Rambaut
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

public class OrRule implements XMLSyntaxRule {

	public OrRule(XMLSyntaxRule a, XMLSyntaxRule b) {
		rules = new XMLSyntaxRule[] {a,b};
	}
	
	public OrRule(XMLSyntaxRule[] rules) {
		this.rules = rules;
	}

	public XMLSyntaxRule[] getRules() {
		return rules;
	}

	/**
	 * Returns true if the rule is satisfied for the given XML object.
	 */
	public boolean isSatisfied(XMLObject xo) {
		for (int i = 0; i < rules.length; i++) {
			if (rules[i].isSatisfied(xo)) return true;
		}
		return false;
	}

    public boolean containsAttribute(String name) {
        for (int i = 0; i < rules.length; i++) {
            if( rules[i].containsAttribute((name)) ) {
                return true;
            }
        }
        return false;
    }

    /**
	 * Describes the rule.
	 */
	public String ruleString() {
		String ruleString = "(" + rules[0].ruleString();
		for (int i = 1; i < rules.length; i++) {
			ruleString += "| " + rules[i].ruleString();
		}
		return ruleString + ")";
	}
	
	/**
	 * Describes the rule.
	 */
	public String htmlRuleString(XMLDocumentationHandler handler) {
		String html = "<div class=\"requiredcompoundrule\">At least one of:";
		for (int i = 0; i < rules.length; i++) {
			html += rules[i].htmlRuleString(handler);
		}	
		html += "</div>";
		return html;	
	}
	
	/**
	 * Describes the rule.
	 */
	public String ruleString(XMLObject xo) {
	
		for (int i = 0; i < rules.length; i++) {
			if (!rules[i].isSatisfied(xo)) return rules[i].ruleString(xo);
		}
		return ruleString();
	}
	
	/**
	 * @return a set containing the required types of this rule.
	 */
	public Set getRequiredTypes() { 
		Set set = new TreeSet(ClassComparator.INSTANCE);
		for (int i = 0; i < rules.length; i++) {
			set.addAll(rules[i].getRequiredTypes());
		}
		return set;
	}

	XMLSyntaxRule[] rules;
}
