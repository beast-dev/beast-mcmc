/*
 * XMLExporter.java
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

package dr.app.bss;

import java.io.IOException;
import java.io.StringWriter;

import dr.app.beauti.util.XMLWriter;
import dr.evolution.alignment.SimpleAlignment;
import dr.evolution.datatype.DataType;
import dr.evolution.util.Taxon;
import dr.evoxml.AlignmentParser;
import dr.evoxml.SequenceParser;
import dr.evoxml.TaxonParser;
import dr.util.Attribute;
import dr.xml.XMLParser;

public class XMLExporter {

	public XMLExporter() {

	}// END: Constructor
	
	public String exportAlignment(SimpleAlignment alignment) throws IOException {
		
		StringWriter sw = new StringWriter();
		XMLWriter writer = new XMLWriter(sw);

//		TODO: if we keep the taxa element than lets also write dates		
//		writer.writeOpenTag(TaxaParser.TAXA, // tagname
//				new Attribute[] { // attributes[]
//				new Attribute.Default<String>(XMLParser.ID, TaxaParser.TAXA) });
//		
//		for (int i = 0; i < alignment.getTaxonCount(); i++) {
//
//			Taxon taxon = alignment.getTaxon(i);
//
//			writer.writeTag(
//					TaxonParser.TAXON, // tagname
//					new Attribute[] { // attributes[]
//					new Attribute.Default<String>(XMLParser.ID, taxon.getId()) },
//					true // close
//			);
//
////			System.out.println(taxon.getAttribute(Utils.ABSOLUTE_HEIGHT));
//			
////			writer.writeCloseTag(TaxonParser.TAXON);
//			
//		}// END: taxon loop
//		
//		writer.writeCloseTag(TaxaParser.TAXA);
//		
//		writer.writeBlankLine();
		
		writer.writeOpenTag(AlignmentParser.ALIGNMENT, // tagname
				new Attribute[] { // attributes[]
				new Attribute.Default<String>(XMLParser.ID, AlignmentParser.ALIGNMENT),
				new Attribute.Default<String>(DataType.DATA_TYPE, alignment.getDataType().getDescription())
		});
		
		for (int i = 0; i < alignment.getSequenceCount(); i++) {
			
			Taxon taxon = alignment.getTaxon(i);
			
			writer.writeOpenTag(SequenceParser.SEQUENCE);
			writer.writeIDref(TaxonParser.TAXON, taxon.getId());
			writer.writeText(alignment.getSequence(i).getSequenceString());
			writer.writeCloseTag(SequenceParser.SEQUENCE);
			
		}//END: sequences loop
		
		writer.writeCloseTag(AlignmentParser.ALIGNMENT);
		writer.close();
		
		return sw.toString();
	}//END: exportAlignment
	
}//END: class
