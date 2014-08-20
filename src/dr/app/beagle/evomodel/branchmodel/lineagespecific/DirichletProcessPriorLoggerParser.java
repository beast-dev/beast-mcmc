package dr.app.beagle.evomodel.branchmodel.lineagespecific;

import dr.inference.distribution.ParametricMultivariateDistributionModel;
import dr.inference.model.CompoundParameter;
import dr.inference.model.Parameter;
import dr.xml.AbstractXMLObjectParser;
import dr.xml.XMLObject;
import dr.xml.XMLParseException;
import dr.xml.XMLSyntaxRule;

public class DirichletProcessPriorLoggerParser extends AbstractXMLObjectParser {

	public static final String DPP_LOGGER = "dppLogger";
	public static final String PRECISION = "precision";
	
	@Override
	public Object parseXMLObject(XMLObject xo) throws XMLParseException {

//		ParametricMultivariateDistributionModel baseModel = (ParametricMultivariateDistributionModel) xo.getChild(ParametricMultivariateDistributionModel.class);
		Parameter precisionParameter = (Parameter)xo.getElementFirstChild(PRECISION); 
		CompoundParameter uniquelyRealizedParameters = (CompoundParameter)xo.getChild(CompoundParameter.class);
		Parameter categoriesParameter = (Parameter)xo.getElementFirstChild(DirichletProcessPriorParser.CATEGORIES); 
		
		return new DirichletProcessPriorLogger(precisionParameter, categoriesParameter, uniquelyRealizedParameters);
	}// END: parseXMLObject

	@Override
	public XMLSyntaxRule[] getSyntaxRules() {
		// TODO Auto-generated method stub
		return null;
	}

	@Override
	public String getParserName() {
		return DPP_LOGGER;
	}

	@Override
	public String getParserDescription() {
		return DPP_LOGGER;
	}

	@Override
	public Class getReturnType() {
		return DirichletProcessPriorLogger.class;
	}

}// END: class
