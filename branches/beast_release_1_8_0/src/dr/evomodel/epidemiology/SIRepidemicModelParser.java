package dr.evomodel.epidemiology;

import dr.evolution.util.Units;
import dr.evoxml.util.XMLUnits;
import dr.inference.model.Parameter;
import dr.xml.*;

/**
 * Parses an element from an XMLObject into SIRepidemicModel. 
 * @author Daniel Wilson
 */
public class SIRepidemicModelParser extends AbstractXMLObjectParser {
	
	public static String POPULATION_SIZE = "populationSize";
	public static String SIREPI_MODEL = "epidemicSIR";
	
	public static String GROWTH_RATE = "growthRate";
	public static String DOUBLING_TIME = "doublingTime";
	public static String TPEAK = "peakTime";
	public static String GAMMA = "clearanceRate";
	
	public static String MIN_PREVALENCE = "minPrevalence";
	
	public String getParserName() { return SIREPI_MODEL; }
	
	public Object parseXMLObject(XMLObject xo) throws XMLParseException {
		Units.Type units = XMLUnits.Utils.getUnitsAttr(xo);
		
		XMLObject cxo = (XMLObject)xo.getChild(POPULATION_SIZE);
		Parameter N0Param = (Parameter)cxo.getChild(Parameter.class);		
		
		boolean usingGrowthRate = true;		
		Parameter rParam = null;		
		if (xo.getChild(GROWTH_RATE) != null) {			
			cxo = (XMLObject)xo.getChild(GROWTH_RATE);
			rParam = (Parameter)cxo.getChild(Parameter.class);
		} else {
			cxo = (XMLObject)xo.getChild(DOUBLING_TIME);
			rParam = (Parameter)cxo.getChild(Parameter.class);
			usingGrowthRate = false;
		}
		
		cxo = (XMLObject)xo.getChild(TPEAK);
		Parameter tpeakParam = (Parameter)cxo.getChild(Parameter.class);
		
		cxo = (XMLObject)xo.getChild(GAMMA);
		Parameter gammaParam = (Parameter)cxo.getChild(Parameter.class);
		
		double minPrevalence = xo.getDoubleAttribute(MIN_PREVALENCE);
		
		return new SIRepidemicModel(N0Param, rParam, tpeakParam, gammaParam, units, usingGrowthRate, minPrevalence);
	}
	
	//************************************************************************
	// AbstractXMLObjectParser implementation
	//************************************************************************
	
	public String getParserDescription() {
		return "SIR epidemic model using RK.";
	}
	
	public Class getReturnType() { return SIRepidemicModel.class; }
	
	public XMLSyntaxRule[] getSyntaxRules() { return rules; }
	
	private XMLSyntaxRule[] rules = new XMLSyntaxRule[] {
	XMLUnits.SYNTAX_RULES[0],
	new ElementRule(POPULATION_SIZE, 
					new XMLSyntaxRule[] { new ElementRule(Parameter.class) },
					"This parameter represents the present day effective population size."),
	new XORRule(
				new ElementRule(GROWTH_RATE, 
								new XMLSyntaxRule[] { new ElementRule(Parameter.class) },
								"This parameter determines the rate of growth during the exponential phase. See exponentialGrowth for details."),
				new ElementRule(DOUBLING_TIME, 
								new XMLSyntaxRule[] { new ElementRule(Parameter.class) },
								"This parameter determines the doubling time at peak growth rate.")
	),
	new ElementRule(TPEAK, 
					new XMLSyntaxRule[] { new ElementRule(Parameter.class) }, 
					"This parameter represents date at which the epidemic reached its peak (positive times indicate the peak was in the past)."),
	new ElementRule(GAMMA, 
					new XMLSyntaxRule[] { new ElementRule(Parameter.class) }, 
					"This parameter represents the rate of loss of infectiousness."),
	AttributeRule.newDoubleRule(MIN_PREVALENCE)
	};
}
