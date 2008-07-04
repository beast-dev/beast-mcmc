package dr.evolution.util;

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

public class RandomTaxaSample extends Taxa{
	
	public static final String RANDOM_TAXA_SAMPLE = "randomTaxaSample";
	public static final String SAMPLE="sample";
	
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
			
			Logger.getLogger("dr.evolution").info("Random Taxa Selected");
								
			for(int i = 0; i < n; i++){
				int randomIndex = MathUtils.nextInt(indexes.size());
				
				Taxon selectedTaxon = population.getTaxon(indexes.get(randomIndex));
				
				sample.addTaxon(selectedTaxon);
				
				indexes.remove(randomIndex);
			
				Logger.getLogger("dr.evolution").info(Integer.toString(i + 1) + ": " + selectedTaxon.toString());
			}
			
			return sample;
		}

		public String getParserName() {
			return RANDOM_TAXA_SAMPLE;
		}
		
	};

}
