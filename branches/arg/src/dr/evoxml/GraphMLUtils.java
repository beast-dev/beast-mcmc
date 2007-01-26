package dr.evoxml;

import dr.evomodel.tree.ARGModel;
import org.jdom.Attribute;
import org.jdom.DataConversionException;
import org.jdom.Element;

import java.util.ArrayList;
import java.util.List;

/**
 * Created by IntelliJ IDEA.
 * User: msuchard
 * Date: Jan 15, 2007
 * Time: 10:11:59 AM
 * To change this template use File | Settings | File Templates.
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
			int length = attributes.size();
			for (Attribute attribute : attributes) {
				String name = attribute.getName();
				if (name.compareTo(ARGModel.ID_ATTRIBUTE) != 0) {
					if (name.compareTo(ARGModel.IS_TIP) == 0) {
						try {
							if (attribute.getBooleanValue()) {
								tipNames.add(nodeName);
							}
						} catch (DataConversionException e) {
							e.printStackTrace();  //To change body of catch statement use File | Settings | File Templates.
						}
						cnt++;
					} else {
						if (!started) {
							space(sb);
							startAttribute(sb);
							started = true;
						}
						sb.append(name + "=" + attribute.getValue());
						cnt++;
						if (cnt < length)
							nextAttribute(sb);
						else
							endAttribute(sb);
					}
				}
			}
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
                sb.append("label=\""+partitions+"\"");
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

}
