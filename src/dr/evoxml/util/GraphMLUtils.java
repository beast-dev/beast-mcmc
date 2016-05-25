/*
 * GraphMLUtils.java
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

package dr.evoxml.util;

import dr.evomodel.arg.ARGModel;
import org.jdom.Attribute;
import org.jdom.DataConversionException;
import org.jdom.Element;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;

/**
 * @author Marc Suchard
 */
public class GraphMLUtils {

	public static final String GRAPH_NAME = "digraph";
	public static final String SPACE = " ";
	public static final String NEW_LINE = "\n";
	public static final String START_SECTION = "{";
	public static final String END_SECTION = "}";
	public static final String START_ATTRIBUTE = "[";
	public static final String END_ATTRIBUTE = "]";
	public static final String NEXT_ATTRIBUTE = ",";
	public static final String TAB = "\t";
	public static final String END_LINE = ";\n";

	public static final boolean printLengths = false;

	private static void space(StringBuilder sb) {
		sb.append(SPACE);
	}


	private static void endLine(StringBuilder sb) {
		sb.append(END_LINE);
	}

	private static void newLine(StringBuilder sb) {
		sb.append(NEW_LINE);
	}

	private static void tab(StringBuilder sb) {
		sb.append(TAB);
	}

	private static void startSection(StringBuilder sb) {
		sb.append(START_SECTION);
	}

	private static void endSection(StringBuilder sb) {
		sb.append(END_SECTION);
	}

	private static void startAttribute(StringBuilder sb) {
		sb.append(START_ATTRIBUTE);
	}

	private static void endAttribute(StringBuilder sb) {
		sb.append(END_ATTRIBUTE);
	}

	private static void nextAttribute(StringBuilder sb) {
		sb.append(NEXT_ATTRIBUTE);
	}

	public static String dotFormat(Element graphElement) {
		StringBuilder sb = new StringBuilder();
		String graphName = graphElement.getAttributeValue(ARGModel.ID_ATTRIBUTE);
		sb.append(GRAPH_NAME);
		space(sb);
		if (graphName != null) {
			sb.append(graphName);
			space(sb);
		}
		startSection(sb);
		newLine(sb);
		tab(sb);
		sb.append(ARGModel.GRAPH_SIZE);
		endLine(sb);
//		newLine(sb);
		tab(sb);
		sb.append(ARGModel.DOT_EDGE_DEF);
		endLine(sb);
//		newLine(sb);
		tab(sb);
		sb.append(ARGModel.DOT_NODE_DEF);
		endLine(sb);
		// Fill in graph details
		List<Element> nodeList = graphElement.getChildren(ARGModel.NODE_ELEMENT);
		List<String> tipNames = new ArrayList<String>();
		for (Element nodeElement : nodeList) {
			String nodeName = nodeElement.getAttributeValue(ARGModel.ID_ATTRIBUTE);
			tab(sb);
			sb.append(nodeName);
			List<Attribute> attributes = nodeElement.getAttributes();
			int cnt = 1;
			boolean started = false;
			boolean isTip = false;
			int length = attributes.size();
			for (Attribute attribute : attributes) {
				String name = attribute.getName();
				if (name.compareTo(ARGModel.ID_ATTRIBUTE) != 0) {
					if (name.compareTo(ARGModel.IS_TIP) == 0) {
						isTip = true;
						try {
							if (attribute.getBooleanValue()) {
								tipNames.add(nodeName);
							}
						} catch (DataConversionException e) {
							e.printStackTrace();  //To change body of catch statement use File | Settings | File Templates.
						}
						cnt++;
					} else {
						if (!ignore(name)) {
							if (!started) {
								space(sb);
								startAttribute(sb);
								started = true;
							} else
								nextAttribute(sb);
							sb.append(translate(name) + "=" + attribute.getValue());
							cnt++;
//							if (cnt < length)
//								nextAttribute(sb);
//							else {

//							}
						}
					}
				}
			}
			if (!isTip) {
				if (!started) {
					space(sb);
					startAttribute(sb);
					started = true;
				} else {
					nextAttribute(sb);
				}
				sb.append("label=\"\",shape=circle,height=0.02,width=0.2,fontsize=1");


			}
			if (started) {
				endAttribute(sb);


			}
//			if (!isTip) {
//
//			}
			endLine(sb);
		}
		List<Element> edgeList = graphElement.getChildren(ARGModel.EDGE_ELEMENT);
		for (Element edgeElement : edgeList) {
			String fromName = edgeElement.getAttributeValue(ARGModel.EDGE_FROM);
			String toName = edgeElement.getAttributeValue(ARGModel.EDGE_TO);
			tab(sb);
			sb.append(fromName + " -> " + toName);
			String edgeLength = edgeElement.getAttributeValue(ARGModel.EDGE_LENGTH);
			if (edgeLength != null && printLengths) {
				space(sb);
				startAttribute(sb);
				sb.append(ARGModel.EDGE_LENGTH + "=" + edgeLength + ",weight=1000.0");
				endAttribute(sb);
			}
			String partitions = edgeElement.getAttributeValue(ARGModel.EDGE_PARTITIONS);
			if (partitions != null) {
				space(sb);
				startAttribute(sb);
				sb.append("label=\"" + partitions + "\"");
				endAttribute(sb);
			}
			endLine(sb);
		}
		//newLine(sb);
		if (tipNames.size() > 0) {
			tab(sb);
			startSection(sb);
			sb.append("rank=same;");
			for (String name : tipNames) {
				space(sb);
				sb.append(name);
			}
			endSection(sb);
			newLine(sb);
		}

		endSection(sb);
		return sb.toString();
	}

	private static class TranslationMap extends HashMap<String, String> {

		TranslationMap() {
			put("taxonName", "label");
		}

	}

	private static class IgnoreList extends ArrayList<String> {
		IgnoreList() {
			add("nodeHeight");
			add("isRoot");
//			add()
		}
	}

	private static TranslationMap translation = new TranslationMap();

	private static IgnoreList ignoreList = new IgnoreList();

	private static String translate(String text) {
		String newText = null;
		if ((newText = translation.get(text)) == null)
			return text;
		return newText;
//		if ()
//		System.err.println("t size = "+translation.size());
//		return text;
	}

	private static boolean ignore(String text) {
		return ignoreList.contains(text);
	}

}
