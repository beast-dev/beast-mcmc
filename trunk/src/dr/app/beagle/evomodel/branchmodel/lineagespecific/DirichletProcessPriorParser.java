package dr.app.beagle.evomodel.branchmodel.lineagespecific;

import java.util.ArrayList;
import java.util.List;

import dr.inference.distribution.ParametricMultivariateDistributionModel;
import dr.inference.model.CompoundParameter;
import dr.inference.model.Parameter;
import dr.xml.AbstractXMLObjectParser;
import dr.xml.ElementRule;
import dr.xml.XMLObject;
import dr.xml.XMLParseException;
import dr.xml.XMLSyntaxRule;

public class DirichletProcessPriorParser extends AbstractXMLObjectParser {

	public static final String DIRICHLET_PROCESS_PRIOR = "dirichletProcessModel";
	public static final String BASE_MODELS = "baseModels";
	public static final String CONCENTRATION = "concentration";
	public static final String CATEGORIES = "categories";
	
	public static final String REALIZED_VALUES = "realizedValues";
	
	@Override
	public String getParserName() {
		return DIRICHLET_PROCESS_PRIOR;
	}

	@Override
	public Object parseXMLObject(XMLObject xo) throws XMLParseException {

		Parameter categoriesParameter =  (Parameter)xo.getElementFirstChild(CATEGORIES);
		CompoundParameter realizedParameters = (CompoundParameter)xo.getChild(CompoundParameter.class);

		XMLObject cxo = xo.getChild(BASE_MODELS);
		List<ParametricMultivariateDistributionModel> baseModels = new ArrayList<ParametricMultivariateDistributionModel>();
				for (int i = 0; i < cxo.getChildCount(); i++) {

					ParametricMultivariateDistributionModel baseModel = (ParametricMultivariateDistributionModel) cxo.getChild(i);
					baseModels.add(baseModel);
					
				}//END: base models loop

		Parameter gamma = (Parameter) xo.getElementFirstChild(CONCENTRATION);

		return new DirichletProcessPrior(categoriesParameter, //
				realizedParameters, //
				baseModels, //
				gamma);
	}//END: parseXMLObject

	@Override
	public XMLSyntaxRule[] getSyntaxRules() {
		return new XMLSyntaxRule[]{
				
				new ElementRule(CATEGORIES,
	                    new XMLSyntaxRule[] { new ElementRule(Parameter.class, false) }), // categories assignments
				
				new ElementRule(CompoundParameter.class, false), // realized parameters
				
		        new ElementRule(BASE_MODELS,
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
