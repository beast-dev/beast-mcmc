package dr.app.beagle.tools.parsers;

import java.util.ArrayList;

import dr.app.beagle.tools.BeagleSequenceSimulator;
import dr.app.beagle.tools.Partition;
import dr.evolution.alignment.Alignment;
import dr.xml.AbstractXMLObjectParser;
import dr.xml.AttributeRule;
import dr.xml.ElementRule;
import dr.xml.XMLObject;
import dr.xml.XMLParseException;
import dr.xml.XMLSyntaxRule;

/**
 * @author Filip Bielejec
 * @version $Id$
 */
public class BeagleSequenceSimulatorParser extends AbstractXMLObjectParser {

	public static final String BEAGLE_SEQUENCE_SIMULATOR = "beagleSequenceSimulator";
	public static final String REPLICATIONS = "replications";
	
	public String getParserName() {
		return BEAGLE_SEQUENCE_SIMULATOR;
	}

	@Override
	public String getParserDescription() {
		return "Beagle sequence simulator";
	}

	@Override
	public Class<Alignment> getReturnType() {
		return Alignment.class;
	}

	@Override
	public XMLSyntaxRule[] getSyntaxRules() {
		
		return new XMLSyntaxRule[] {
				AttributeRule.newIntegerRule(REPLICATIONS), 
				new ElementRule(Partition.class, 1, Integer.MAX_VALUE)
				};
		
	}//END: getSyntaxRules

	//TODO: fix parser to work with partitions
	@Override
	public Object parseXMLObject(XMLObject xo) throws XMLParseException {

		int replications = xo.getIntegerAttribute(REPLICATIONS);
		
		ArrayList<Partition> partitionsList = new ArrayList<Partition>();
		for (int i = 0; i < xo.getChildCount(); i++) {

			Partition partition = (Partition) xo.getChild(i);
			
			if (partition.from > replications) {
				throw new XMLParseException(
						"illegal 'from' attribute in " + PartitionParser.PARTITION + " element");
			}

			if (partition.to > replications) {
				throw new XMLParseException(
						"illegal 'to' attribute in " + PartitionParser.PARTITION + " element");
			}

			if (partition.to == -1) {
				partition.to = replications;
			}
			
			partitionsList.add(partition);

		}// END: partitions loop

		BeagleSequenceSimulator s = new BeagleSequenceSimulator(partitionsList, //
				replications //
		);

		return s.simulate();
	}// END: parseXMLObject

}// END: class
