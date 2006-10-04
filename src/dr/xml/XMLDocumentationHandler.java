/*
 * XMLDocumentationHandler.java
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

import dr.math.MathUtils;

import java.io.File;
import java.io.FileWriter;
import java.io.PrintWriter;
import java.util.*;

public class XMLDocumentationHandler {

	public XMLDocumentationHandler(XMLParser parser) {
		
		this.parser = parser;
		
		Iterator iterator = parser.getParsers();
		while (iterator.hasNext()) {
			XMLObjectParser xmlparser = (XMLObjectParser)iterator.next();
			
			XMLSyntaxRule[] rules = xmlparser.getSyntaxRules();
			
			if (rules != null) {
				for (int i =0; i < rules.length; i++) {
					Set requiredTypesForRule = rules[i].getRequiredTypes();
					requiredTypes.addAll(requiredTypesForRule);
				}
			}
		}
	}
	
	public void outputElements(File directory) throws java.io.IOException {
	
		PrintWriter writer = new PrintWriter(new FileWriter(new File(directory, "index.html")));
	
		writer.println("<html>");
		writer.println("<head>");  
		writer.println("  <link rel=\"stylesheet\" href=\"../beast.css\">");
		writer.println("  <title>BEAST elements</title>");
		writer.println("</head>");
		writer.println("<body>");
		writer.println("<h1>BEAST elements</h1>");
		writer.println("<p>");
		writer.println("The following is a list of valid elements in a beast file.<br>");
		writer.println("<span class=\"required\">&nbsp;&nbsp;&nbsp;&nbsp;</span> required<br>"); 
		writer.println("<span class=\"optional\">&nbsp;&nbsp;&nbsp;&nbsp;</span> optional<br>"); 
		writer.println("</p>");
		
		Iterator iterator = parser.getParsers();
		while (iterator.hasNext()) {
			XMLObjectParser xmlParser = (XMLObjectParser)iterator.next();
			writer.println(xmlParser.toHTML(this));
			System.out.println("  outputting HTML for element " + xmlParser.getParserName());
		}
		
		writer.println("</body>");
		writer.println("</html>");
		writer.flush();
		writer.close();
	}
	
	/**
	 * Outputs an example of a particular element, using the syntax information.
	 */
	public void outputExampleXML(PrintWriter writer, XMLObjectParser parser) {
		
		writer.println("<pre>");
		if (parser.hasExample()) {
			outputHTMLSafeText(writer, parser.getExample());
		} else {
			outputExampleXML(writer, parser, 0);
		}
		writer.println("</pre>");
	}
	
	public void outputHTMLSafeText(PrintWriter writer, String text) {
		for (int i =0; i < text.length(); i++) {
			char c = text.charAt(i);
			switch (c) {
				case '<': writer.print("&lt;"); break;
				case '>': writer.print("&gt;"); break;
				case '&': writer.print("&amp;"); break;
				default: writer.print(c); break;
			}
		}
	}
	
	/**
	 * Outputs an example of a particular element, using the syntax information.
	 */
	public void outputExampleXML(PrintWriter writer, XMLObjectParser parser, int level) {
		outputElementRules(writer, parser.getParserName(), parser.getSyntaxRules(), level);
	}
	
	public void stochasticCollectRules(XMLSyntaxRule[] allRules, ArrayList attributeList, ArrayList elementList) {
		
		if (allRules != null) {
			for (int i =0; i < allRules.length; i++) {
			
				XMLSyntaxRule rule = allRules[i];
			
				if (rule instanceof AttributeRule) {
					attributeList.add(rule);
				} else if (rule instanceof ElementRule) {
					int min = ((ElementRule)rule).getMin();
					int max = Math.max(min, Math.min(5, ((ElementRule)rule).getMax()));
					
					int numRules = min;
					if (max != min) numRules = MathUtils.nextInt(max-min) + min;
					for (int j =0; j < numRules; j++) {
						elementList.add((ElementRule)rule);
					}
				} else if (rule instanceof XORRule) {
					XORRule xorRule = (XORRule)rule;
					XMLSyntaxRule[] rules = xorRule.getRules();
					int ruleIndex = MathUtils.nextInt(rules.length);
					stochasticCollectRules(new XMLSyntaxRule[] { rules[ruleIndex] }, attributeList, elementList);
				} else if (rule instanceof OrRule) {
					OrRule orRule = (OrRule)rule;
					XMLSyntaxRule[] rules = orRule.getRules();
					int ruleIndex = MathUtils.nextInt(rules.length);
					stochasticCollectRules(new XMLSyntaxRule[] { rules[ruleIndex] }, attributeList, elementList);
				} 
			}
		}
	}
	
	/**
	 * Outputs an example of a rule, using the syntax information.
	 */
	public void outputExampleXML(PrintWriter writer, AttributeRule rule, int level) {
		writer.print(" " + rule.getName()+"=\"");
		if (rule.hasExample()) {
			writer.print(rule.getExample());
		} else {
			outputAttributeValue(writer, rule.getAttributeClass());
		}
		writer.print("\"");
	}
	
	/**
	 * Outputs an example of a rule, using the syntax information.
	 */
	public void outputExampleXML(PrintWriter writer, ElementRule rule, int level) {
		
		if (rule.getElementClass() == null) {
			if (rule.getName() == null) System.err.println(rule + " has a null name");
			outputElementRules(writer, rule.getName(), rule.getRules(), level);
		} else {
			if (rule.hasExample()) {
				writer.println(spaces(level+1) + rule.getExample());
			} else {
				outputExampleXML(writer, rule.getElementClass(), level + 1);
			}
		}
	}
	
	public void outputElementRules(PrintWriter writer, String name, XMLSyntaxRule[] rules, int level) {
		
		ArrayList attributeList = new ArrayList();
		ArrayList elementList = new ArrayList();
		stochasticCollectRules(rules, attributeList, elementList);
		
		writer.print(spaces(level) + "&lt;" + name);
		// write out the attributes
		for (int i =0; i < attributeList.size(); i++) {
			XMLSyntaxRule rule = (XMLSyntaxRule)attributeList.get(i);
			outputExampleXML(writer, (AttributeRule)rule, level + 1); 
		}
		if (elementList.size() > 0) {
			writer.println("&gt;");
			// write out the elements
			for (int i =0; i < elementList.size(); i++) {
				ElementRule rule = (ElementRule)elementList.get(i);
				outputExampleXML(writer, rule, level + 1);
			}
			writer.println(spaces(level) + "&lt;/" + name + "&gt;");
		} else {
			writer.println("/&gt;");
		}
	}
	
	public void outputExampleXML(PrintWriter writer, Class c, int level) {
			
		if (c == String.class) {
			writer.println(spaces(level) + "foo");
		} else if (c == Double.class) {
			writer.println(spaces(level) + "1.0");
		} else if (c == Integer.class) {
			writer.println(spaces(level) + "1");
		} else if (c == Boolean.class) {
			writer.println(spaces(level) + "true");
		} else if (c == Double[].class) {
			writer.println(spaces(level) + "0.5 1.0");
		} else if (c == String[].class) {
			writer.println(spaces(level) + "foo bar");
		} else {
			if (c == null) { throw new RuntimeException("Class is null"); }
			XMLObjectParser randomParser = getRandomParser(c);
			if (randomParser == null) {
				writer.println(spaces(level) + "ERROR!");
			} else {
				if (level > 1) {
					writer.println(spaces(level) + "&lt;" + randomParser.getParserName() + 
						" idref=\"" + randomParser.getParserName() + (MathUtils.nextInt(10)+1) + "\"/&gt;");
				} else {
					outputExampleXML(writer, randomParser, level);
				}
			}
		}
		
	}
	
	public void outputAttributeValue(PrintWriter writer, Class c) {
		if (c == String.class) {
			writer.print("foo");
		} else if (c == Double.class) {
			writer.print("1.0");
		} else if (c == Integer.class) {
			writer.print("1");
		} else if (c == Boolean.class) {
			writer.print("true");
		} else if (c == Double[].class) {
			writer.print("0.5 1.0");
		} else if (c == String[].class) {
			writer.print("foo bar");
		} else throw new RuntimeException("Class " + c + " not allowed as attribute value");
	}
	
	private String spaces(int level) {
		StringBuffer buffer = new StringBuffer("");
		for (int i =0; i < level; i++) {
			buffer.append(' ');
		}
		return buffer.toString();
	} 
	
	public XMLObjectParser getRandomParser(Class c) {
		
		ArrayList matchingParsers = getMatchingParsers(c);	
		
		if (matchingParsers.size() == 0) return null;
		return (XMLObjectParser)matchingParsers.get(MathUtils.nextInt(matchingParsers.size()));
	}
	
	public final ArrayList getMatchingParsers(Class c) {
		
		ArrayList matchingParsers = new ArrayList();	
		// find all parsers that match this required type
		Iterator i = parser.getParsers();
		while (i.hasNext()) {
			XMLObjectParser xmlParser = (XMLObjectParser)i.next();
			Class returnType = xmlParser.getReturnType();
			if (c.isAssignableFrom(returnType)) {
				matchingParsers.add(xmlParser);
			}
		}
		return matchingParsers;
	}
	
	/**
	 * Outputs all types that appear as required attributes or elements in an HTML table to the given writer.
	 */
	public void outputTypes(File directory) throws java.io.IOException {
		
		PrintWriter writer = new PrintWriter(new FileWriter(new File(directory, "types.html"))); 
		writer.println("<html>");
		writer.println("<head>");  
		writer.println("  <link rel=\"stylesheet\" href=\"../beast.css\">");
		writer.println("  <title>BEAST elements</title>");
		writer.println("</head>");
		writer.println("<h1>BEAST types</h1>");
		writer.println("<p>");
		writer.println("The following is a list of generic types that elements represent in a beast file.<br>");
		writer.println("</p>");
	
		
		// iterate through the types
		Iterator iterator = requiredTypes.iterator();
		while (iterator.hasNext()) {
			
			Class requiredType = (Class)iterator.next();
			
				
			if (requiredType != Object.class) {
			
				String name = ClassComparator.getName(requiredType);
			
				System.out.println("  outputting HTML for generic type " + name);
				
				
				TreeSet matchingParserNames = new TreeSet();
				
				// find all parsers that match this required type
				Iterator i = parser.getParsers();
				while (i.hasNext()) {
					XMLObjectParser xmlParser = (XMLObjectParser)i.next();
					Class returnType = xmlParser.getReturnType();
					if (requiredType.isAssignableFrom(returnType)) {
						matchingParserNames.add(xmlParser.getParserName());
					}
				}
				
				// output table row containing the type and the matching parser names
				writer.println("<div id=\"" + name + "\"><h2>" + name + "</h2>");
				writer.println("<p>");
				writer.println("Elements of this type include:");
				writer.println("</p>");
				i = matchingParserNames.iterator();
				while (i.hasNext()) {
					String parserName = (String)i.next();
					writer.println("<div><a href=\"index.html#" + parserName + "\"> &lt;" + parserName + "&gt;</a></div>");
				}
				writer.println("</div>");
			}
			
		}
		
		writer.println("</html>");
		writer.flush();
		writer.close();
	}
	
	
	public Set getParsersForClass(Class returnType) {
		
		TreeSet set = new TreeSet();
		return set;
	}
	
	public String getHTMLForClass(Class c) {
		String name = ClassComparator.getName(c);
		return "<A HREF=\"types.html#" + name + "\">" + name + "</A>";
	}
	
	class SetHash {

	    private HashMap table;
	    
	    public SetHash() { table = new HashMap(); }

	    public final void put(Object key, XMLObjectParser o) {
			Set set = (Set)table.get(key);

			if (set != null) {
			    set.add(o);
			} else {
			    TreeSet newSet = new TreeSet();
			    newSet.add(o);
			    table.put(key, newSet);
			}
	    }

	    public final Set keySet() { return table.keySet(); }

	    public final Object[] getArray(Object key) { return getSortedSet(key).toArray(); }

	    public final SortedSet getSortedSet(Object key) { return (SortedSet)table.get(key); }
	}
	
	class XMLObjectParserComparator implements Comparator {
		
		public int compare(Object c1, Object c2) {
	
			String name1 = ((XMLObjectParser)c1).getParserName().toUpperCase();
			String name2 = ((XMLObjectParser)c2).getParserName().toUpperCase();
		
			return name1.compareTo(name2);
		}
	}
	
	private Set requiredTypes = new TreeSet(ClassComparator.INSTANCE);
	private XMLParser parser = null;
}
