package dr.app.beagle.evomodel.branchmodel.lineagespecific;

import dr.inference.distribution.ParametricMultivariateDistributionModel;
import dr.inference.model.CompoundParameter;
import dr.inference.model.Parameter;
import dr.xml.AbstractXMLObjectParser;
import dr.xml.ElementRule;
import dr.xml.XMLObject;
import dr.xml.XMLParseException;
import dr.xml.XMLSyntaxRule;

public class DirichletProcessPriorParser extends AbstractXMLObjectParser {

	public static final String DIRICHLET_PROCESS_PRIOR = "dirichletProcessPrior";
	public static final String BASE_MODEL = "baseModel";
	public static final String CONCENTRATION = "concentration";
	public static final String CATEGORIES = "categories";
	
	@Override
	public String getParserName() {
		return DIRICHLET_PROCESS_PRIOR;
	}

	@Override
	public Object parseXMLObject(XMLObject xo) throws XMLParseException {

		Parameter categoriesParameter =  (Parameter)xo.getElementFirstChild(CATEGORIES);
		CompoundParameter uniquelyRealizedParameters = (CompoundParameter)xo.getChild(CompoundParameter.class);
		ParametricMultivariateDistributionModel baseModel = (ParametricMultivariateDistributionModel) xo.getElementFirstChild(BASE_MODEL);
		Parameter gamma = (Parameter) xo.getElementFirstChild(CONCENTRATION);

		return new DirichletProcessPrior(categoriesParameter, //
				uniquelyRealizedParameters, //
				baseModel, //
				gamma);
	}//END: parseXMLObject

	@Override
	public XMLSyntaxRule[] getSyntaxRules() {
		return new XMLSyntaxRule[]{
				
				new ElementRule(CATEGORIES,
	                    new XMLSyntaxRule[] { new ElementRule(Parameter.class, false) }), // categories assignments
				
				new ElementRule(CompoundParameter.class, false), // realized parameters
				
		        new ElementRule(BASE_MODEL,
                        new XMLSyntaxRule[] {
                                new ElementRule(ParametricMultivariateDistributionModel.class, 1, Integer.MAX_VALUE),
                        }
		        ), // base models
		        
		        new ElementRule(CONCENTRATION,
	                    new XMLSyntaxRule[] { new ElementRule(Parameter.class, false) }),// gamma
				
		};
	}

	@Override
	public String getParserDescription() {
		return DIRICHLET_PROCESS_PRIOR;
	}

	@Override
	public Class getReturnType() {
		return DirichletProcessPrior.class;
	}
	
}//END: class
