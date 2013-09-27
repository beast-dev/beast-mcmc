package dr.app.bss;

import java.io.BufferedWriter;
import java.io.OutputStreamWriter;
import java.io.PrintWriter;
import java.io.Writer;

import dr.app.beauti.util.XMLWriter;
import dr.evolution.alignment.SimpleAlignment;
import dr.evoxml.TaxaParser;
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
		
		
		buffer.append(writer.toString());
		writer.close();
		
		return buffer.toString();
	}//END: exportAlignment
	
}//END: class
