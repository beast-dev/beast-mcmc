package dr.evolution.util;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.Reader;
import java.io.Writer;
import java.util.ArrayList;
import java.util.logging.Logger;

import dr.math.MathUtils;
import dr.xml.AbstractXMLObjectParser;
import dr.xml.AttributeRule;
import dr.xml.ElementRule;
import dr.xml.XMLObject;
import dr.xml.XMLObjectParser;
import dr.xml.XMLParseException;
import dr.xml.XMLSyntaxRule;

public class RandomTaxaSample extends Taxa {

	public static final String RANDOM_TAXA_SAMPLE = "randomTaxaSample";
	public static final String SAMPLE="sample";
	public static final String PRINT_TAXA = "printTaxa";
	public static final String FILE_NAME= "fileName";

	public RandomTaxaSample(){
	}

	public static XMLObjectParser PARSER = new AbstractXMLObjectParser(){

		public String getParserDescription() {
			return "Randomly samples n taxa from a collection of N taxa";
		}

		public Class getReturnType() {
			return RandomTaxaSample.class;
		}

		public XMLSyntaxRule[] getSyntaxRules() {
			return rules;
		}

		private XMLSyntaxRule[] rules = new XMLSyntaxRule[]{
				AttributeRule.newIntegerRule(SAMPLE),
				new ElementRule(Taxa.class),
		};

		public Object parseXMLObject(XMLObject xo) throws XMLParseException {

			Taxa population = (Taxa)xo.getChild(Taxa.class);
			int n = xo.getIntegerAttribute(SAMPLE);

			int N = population.getTaxonCount();

			if(n <= 0 || n > N){
				throw new XMLParseException("sample must be greater than 0 and less than or equal to the population size");
			}

			RandomTaxaSample sample = new RandomTaxaSample();

			ArrayList<Integer> indexes = new ArrayList<Integer>(N);

			for(int i = 0; i < N; i++)
				indexes.add(i);

			Logger.getLogger("dr.evolution").info("Generating a random taxa sample of size: " + n);

			for(int i = 0; i < n; i++){
				int randomIndex = MathUtils.nextInt(indexes.size());

				Taxon selectedTaxon = population.getTaxon(indexes.get(randomIndex));

				sample.addTaxon(selectedTaxon);

				indexes.remove(randomIndex);
			}

			if(xo.hasAttribute(PRINT_TAXA) && xo.getBooleanAttribute(PRINT_TAXA)){
				String fileName = null;
				if(xo.hasAttribute(FILE_NAME)){
					fileName = xo.getStringAttribute(FILE_NAME);
				}

				if(fileName != null){
					try {
						Writer write;

						File file = new File(fileName);
						String name = file.getName();
						String parent = file.getParent();

						if (!file.isAbsolute()) {
							parent = System.getProperty("user.dir");
						}

						write = new FileWriter(new File(parent, name));

						write.write("<taxa id=\"randomTaxaSample\">\n" );

						for(int i = 0; i < n ;i++){
							write.write("\t<taxon idref=\"" + sample.getTaxonId(i).toString() + "\"/>\n");
						}

						write.write("</taxa id=\"randomTaxaSample\">\n");

						write.flush();
					}catch (IOException fnfe) {
							throw new XMLParseException("File '" + fileName + "' can not be opened for " + getParserName() + " element.");
					}
				}else{
					Logger.getLogger("dr.evomodel").info("<taxa id=\"randomTaxaSample\">");

					for(int i = 0; i < n ;i++){
						Logger.getLogger("dr.evomodel").info("\t<taxon idref=\" " + sample.getTaxonId(i).toString() + " \"> ");
					}
					Logger.getLogger("dr.evomodel").info("</taxa id=\"randomTaxaSample\">");

				}

			}





			return sample;
		}

		public String getParserName() {
			return RANDOM_TAXA_SAMPLE;
		}

	};

}
