package dr.app.bss;

import java.io.BufferedWriter;
import java.io.File;
import java.io.IOException;
import java.io.OutputStreamWriter;
import java.io.PrintStream;
import java.io.PrintWriter;
import java.io.Writer;

import dr.app.beauti.util.XMLWriter;
import dr.evolution.alignment.SimpleAlignment;
import dr.evolution.datatype.DataType;
import dr.evolution.util.Taxon;
import dr.evoxml.AlignmentParser;
import dr.evoxml.SequenceParser;
import dr.evoxml.TaxaParser;
import dr.evoxml.TaxonParser;
import dr.util.Attribute;
import dr.xml.XMLParser;

public class XMLExporter {

	public XMLExporter() {

	}// END: Constructor
	
	public String exportAlignment(SimpleAlignment alignment) throws IOException {
		
		StringBuffer buffer = new StringBuffer();
		File tmp = File.createTempFile("tempfile", ".tmp");
		XMLWriter writer = new XMLWriter(new PrintWriter(tmp));

		writer.writeOpenTag(TaxaParser.TAXA, // tagname
				new Attribute[] { // attributes[]
				new Attribute.Default<String>(XMLParser.ID, TaxaParser.TAXA) });
		
		// TODO: dates
		for (int i = 0; i < alignment.getTaxonCount(); i++) {

			Taxon taxon = alignment.getTaxon(i);

			writer.writeTag(
					TaxonParser.TAXON, // tagname
					new Attribute[] { // attributes[]
					new Attribute.Default<String>(XMLParser.ID, taxon.getId()) },
					false // close
			);

			writer.writeCloseTag(TaxonParser.TAXON);
//			writer.writeBlankLine();
			
		}// END: taxon loop
		
		writer.writeCloseTag(TaxaParser.TAXA);
		
		writer.writeBlankLine();
		
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
		
		buffer.append(writer.toString());
		writer.close();
		
		return buffer.toString();
	}//END: exportAlignment
	
}//END: class
