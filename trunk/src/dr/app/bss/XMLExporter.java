package dr.app.bss;

import java.io.BufferedWriter;
import java.io.OutputStreamWriter;
import java.io.PrintWriter;
import java.io.Writer;

import dr.app.beauti.util.XMLWriter;
import dr.evolution.alignment.SimpleAlignment;
import dr.evolution.util.Taxon;
import dr.evoxml.TaxaParser;
import dr.evoxml.TaxonParser;
import dr.util.Attribute;
import dr.xml.XMLParser;

public class XMLExporter {

	public XMLExporter() {
		
	}//END: Constructor
	
	public String exportAlignment(SimpleAlignment alignment) {
		
		StringBuffer buffer = new StringBuffer();
		XMLWriter writer = new XMLWriter(new BufferedWriter(new OutputStreamWriter(System.out)));

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
			writer.writeBlankLine();
			
		}// END: taxon loop
		
		writer.writeCloseTag(TaxaParser.TAXA);
		
		
		
		
		buffer.append(writer.toString());
		writer.close();
		
		return buffer.toString();
	}//END: exportAlignment
	
}//END: class
