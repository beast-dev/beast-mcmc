/*
 * AttributeRule.java
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

public class AttributeRule implements XMLSyntaxRule {

	AttributeRule() {}

	public void setName(String name) { this.name = name; }
	public void setDescription(String description) { this.description = description; }
	public void setAttributeClass(Class c) { this.c = c; }
	public void setOptional(boolean optional) { this.optional = optional; }

	public boolean getOptional() { return optional; }

	/**
	 * Creates a required attribute rule.
	 */
	private AttributeRule(String name, Class c) {
		this.name = name;
		this.c = c;
		this.optional = false;
		//this.description = null;
	}

	/**
	 * Creates an attribute rule.
	 */
	private AttributeRule(String name, Class c, boolean optional) {
		this.name = name;
		this.c = c;
		this.optional = optional;
	}

	/**
	 * Creates an attribute rule.
	 */
	private AttributeRule(String name, Class c, boolean optional, String description) {
		this.name = name;
		this.c = c;
		this.optional = optional;
		this.description = description;
	}

	public static AttributeRule newIntegerRule(String name) { return new AttributeRule(name, Integer.class); }
    public static AttributeRule newLongIntegerRule(String name) { return new AttributeRule(name, Long.class); }
	public static AttributeRule newDoubleRule(String name) { return new AttributeRule(name, Double.class); }
	public static AttributeRule newDoubleArrayRule(String name) { return new AttributeRule(name, Double[].class); }
	public static AttributeRule newBooleanRule(String name) { return new AttributeRule(name, Boolean.class); }
	public static AttributeRule newStringRule(String name) { return new AttributeRule(name, String.class); }
	public static AttributeRule newStringArrayRule(String name) { return new AttributeRule(name, String[].class); }

	public static AttributeRule newIntegerRule(String name, boolean optional) { return new AttributeRule(name, Integer.class, optional); }
	public static AttributeRule newIntegerArrayRule(String name, boolean optional) { return new AttributeRule(name, Integer[].class, optional); }
    public static AttributeRule newLongIntegerRule(String name, boolean optional) { return new AttributeRule(name, Long.class, optional); }
	public static AttributeRule newDoubleRule(String name, boolean optional) { return new AttributeRule(name, Double.class, optional); }
	public static AttributeRule newDoubleArrayRule(String name, boolean optional) { return new AttributeRule(name, Double[].class, optional); }
	public static AttributeRule newBooleanRule(String name, boolean optional) { return new AttributeRule(name, Boolean.class, optional); }
	public static AttributeRule newStringRule(String name, boolean optional) { return new AttributeRule(name, String.class, optional); }
    public static AttributeRule newStringArrayRule(String name, boolean optional) { return new AttributeRule(name, String[].class, optional); }

	public static AttributeRule newIntegerRule(String name, boolean optional, String description) { return new AttributeRule(name, Integer.class, optional, description); }
    public static AttributeRule newLongIntegerRule(String name, boolean optional, String description) { return new AttributeRule(name, Long.class, optional, description); }
	public static AttributeRule newDoubleRule(String name, boolean optional, String description) { return new AttributeRule(name, Double.class, optional, description); }
	public static AttributeRule newDoubleArrayRule(String name, boolean optional, String description) { return new AttributeRule(name, Double[].class, optional, description); }
	public static AttributeRule newBooleanRule(String name, boolean optional, String description) { return new AttributeRule(name, Boolean.class, optional, description); }
	public static AttributeRule newStringRule(String name, boolean optional, String description) { return new AttributeRule(name, String.class, optional, description); }
    public static AttributeRule newStringArrayRule(String name, boolean optional, String description) { return new AttributeRule(name, String[].class, optional, description); }

	public String getName() { return name; }
	public Class getAttributeClass() { return c; }
	public String getDescription() { return description; }

	public Object getAttribute(XMLObject xo) throws XMLParseException {
		return xo.getAttribute(name);
	}

	public boolean hasDescription() { return description != null; }

	public boolean hasExample() { return false; }

	public String getExample() { return null; }

	/**
	 * @return true if the required attribute of the correct type is present.
	 */
	public boolean isSatisfied(XMLObject xo) {

		if (xo.hasAttribute(name)) {
			try {
				Object obj = xo.getAttribute(name);
				return isCompatible(obj);
			} catch (XMLParseException xpe) { return false; }
		} else if (optional) {
			return true;
		}
		return false;
	}

    public boolean containsAttribute(String name) {
        return name.equals(getName());
    }

    /**
	 * @return a string describing the rule.
	 */
	public String ruleString() {
		String rule = "ATTRIBUTE " + getTypeName() + " " + name;
		if (optional) {
			rule += " OPTIONAL";
		} else {
			rule += " REQUIRED";
		}
		return rule;
	}

	/**
	 * @return a string describing the rule.
	 */
	public String htmlRuleString(XMLDocumentationHandler handler) {
		String rule =
			"<div class=\"" + (optional ? "optional" : "required") + "rule\"> " +
			"Attribute <span class=\"attrname\">" + name + "</span> is " +
			handler.getHTMLForClass(c) + " " +
			(hasDescription() ? "<div class=\"description\">" + description + "</div>" : "") +
			"</div>" ;

		/*if (optional) {
			rule += " <span class=\"optional\">optional</span></div>";
		} else {
			rule += " <span class=\"required\">required</span></div>";
		}*/

		/* <span class=\"symbol\">&#206;</span> */

		return rule;
	}

	public String wikiRuleString(XMLDocumentationHandler handler, String prefix) {
		String rule = prefix + "Attribute " + (optional ? "(optional) " : "") + "<code>" + name + "</code> is " +
			handler.getHTMLForClass(c) + "\n" + prefix + ":" +
			(hasDescription() ? "\"" + description + "\"" : "") + "\n";

		return rule;
	}

	public String markdownRuleString(XMLDocumentationHandler handler, String prefix) {
		String rule = prefix + "* <code>" + name + "</code> " + (optional ? "(optional)" : "") + " is of type " +
				handler.getHTMLForClass(c) + "\n" +
				": " +
				(hasDescription() ? "\"" + description + "\"" : "") + "\n\n";

		return rule;
	}

	/**
	 * @return a string describing the rule.
	 */
	public String ruleString(XMLObject xo) {

		if (xo.hasAttribute(name)) {
			try {
				Object obj = xo.getAttribute(name);
				boolean compatible = isCompatible(obj);
				if (compatible) {
					return ruleString();
				} else {
                    return "ATTRIBUTE " + name + " expected to be of type " + getTypeName();
                }

			} catch (XMLParseException xpe) {
				return xpe.toString();
			}
		}
		return ruleString();
	}

	/**
	 * @return a set containing the required types of this rule.
	 */
	public Set<Class> getRequiredTypes() { return Collections.singleton(c); }

    public boolean isLegalElementName(String elementName) {
        return false;
    }

    public boolean isLegalElementClass(Class c) {
        return false;
    }

    public boolean isLegalSubelementName(String elementName) {
        return false;
    }

    public boolean isAttributeRule() { return true; }

	/**
	 * @return true if the given object is compatible with the required class.
	 */
	private boolean isCompatible(Object o) {

		if (c == null) return true;

		if (c.isInstance(o)) { return true; }

		if (o instanceof String) {

			if (c == Double[].class) {

				return XMLObject.isDoubleArray((String)o, null);
			}
			if (c == Integer[].class) {

				return XMLObject.isIntegerArray((String)o, null);
			}
            if (c == String[].class) {

                return true;
            }
			if (c == Double.class) {
				try {
					Double.parseDouble((String)o);
					return true;
				} catch (NumberFormatException nfe) { return false; }
			}
			if (c == Integer.class) {
				try {
					Integer.parseInt((String)o);
					return true;
				} catch (NumberFormatException nfe) { return false; }
			}
            if (c == Long.class) {
                try {
                    Long.parseLong((String)o);
                    return true;
                } catch (NumberFormatException nfe) { return false; }
            }
			if (c == Float.class) {
				try {
					Float.parseFloat((String)o);
					return true;
				} catch (NumberFormatException nfe) { return false; }
			}
			if (c == Boolean.class) {
				return (o.equals("true") || o.equals("false"));
			}
			if (c == Number.class) {
				try {
					Double.parseDouble((String)o);
					return true;
				} catch (NumberFormatException nfe) { return false; }
			}
		}

		return false;
	}

	/**
	 * @return a pretty name for a class.
	 */
	private String getTypeName() {

		if (c == null) return "Object";
		final String name = c.getName();
        final String cBaseName = name.substring(name.lastIndexOf('.') + 1);
        if( c.isArray() ) {
            return "Array of " + cBaseName;
        }
        return cBaseName;
	}

	private String name;
	private Class c;
	private boolean optional;
	private String description;
}
