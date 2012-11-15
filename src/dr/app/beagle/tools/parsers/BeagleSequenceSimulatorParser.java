/*
 * BeagleSequenceSimulatorParser.java
 *
 * Copyright (C) 2002-2012 Alexei Drummond, Andrew Rambaut & Marc A. Suchard
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

package dr.app.beagle.tools.parsers;

import java.util.ArrayList;
import java.util.logging.Logger;

import dr.app.beagle.tools.BeagleSequenceSimulator;
import dr.app.beagle.tools.Partition;
import dr.evolution.alignment.Alignment;
import dr.evolution.datatype.Codons;
import dr.evolution.datatype.Nucleotides;
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

	@Override
	public Object parseXMLObject(XMLObject xo) throws XMLParseException {

		String msg = "";
		
		int replications = xo.getIntegerAttribute(REPLICATIONS);
		msg += "\n\t" + replications + ( (replications > 1) ? " replications " : " replication");
		
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
				partition.to = replications - 1;
			}
			
			if (partition.ancestralSequence != null) {

				if (partition.ancestralSequence.getLength() != 3 * replications && partition.freqModel.getDataType() instanceof Codons) {

					throw new RuntimeException("Ancestral codon sequence has "
							+ partition.ancestralSequence.getLength() + " characters "
							+ "expecting " + 3 * replications + " characters");

				} else if (partition.ancestralSequence.getLength() != replications && partition.freqModel.getDataType() instanceof Nucleotides) {

					throw new RuntimeException("Ancestral nuleotide sequence has "
							+ partition.ancestralSequence.getLength() + " characters "
							+ "expecting " + replications + " characters");

				}// END: dataType check
			}// END: ancestralSequence check
			
			// TODO: print partition info
//			msg += "\n\t" + "partition" + (i + 1) + " from " + partition.from + " to " + partition.to + " every " + partition.every;
			
			partitionsList.add(partition);
		}// END: partitions loop

		 if (msg.length() > 0) {
	            Logger.getLogger("dr.app.beagle.tools").info("Using Beagle Sequence Simulator: " + msg);
		 }
		
		BeagleSequenceSimulator s = new BeagleSequenceSimulator(partitionsList, //
				replications //
		);

		return s.simulate();
	}// END: parseXMLObject

}// END: class
