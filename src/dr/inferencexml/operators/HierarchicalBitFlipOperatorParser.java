package dr.inferencexml.operators;

import dr.inference.model.Parameter;
import dr.inference.operators.HierarchicalBitFlipOperator;
import dr.inference.operators.MCMCOperator;
import dr.xml.*;

/**
 * @author Gabriela Cybis
 */


	public class HierarchicalBitFlipOperatorParser  extends AbstractXMLObjectParser {

	    public static final String HIERARCHICAL_BIT_FLIP_OPERATOR = "HierarchicalBitFlipOperator";
	    public static final String USES_SUM_PRIOR = "usesPriorOnSum";
	    public static final String H_PARAMETER = "hParameter";
	    public static final String STRATA_PARAMETERS = "strataParameters";
	    
	    
	    
	    public String getParserName() {
	        return HIERARCHICAL_BIT_FLIP_OPERATOR;
	    }

	    public Object parseXMLObject(XMLObject xo) throws XMLParseException {

	        double weight = xo.getDoubleAttribute(MCMCOperator.WEIGHT);
	        boolean usesPriorOnSum = xo.getAttribute(USES_SUM_PRIOR,false);
	        
	        XMLObject cxo = xo.getChild(H_PARAMETER);
	        Parameter hParameter = (Parameter) cxo.getChild(Parameter.class);

	        cxo = xo.getChild(STRATA_PARAMETERS);
	        int NStrata = cxo.getChildCount();
	        
	        Parameter[] strataParameters = new Parameter[NStrata];
	       for (int i = 0; i < NStrata; i++) {
	            Parameter parameter = (Parameter) cxo.getChild(i);
	            strataParameters[i] = parameter;
	       }
	
	       
	        return new HierarchicalBitFlipOperator(hParameter,  strataParameters,NStrata, weight, usesPriorOnSum);
	    }

	    //************************************************************************
	    // AbstractXMLObjectParser implementation
	    //************************************************************************

	    public String getParserDescription() {
	        return "This element returns a bit-flip operator on a set of hierarchical and strata parameters simulatneously.";
	    }

	    public Class getReturnType() {
	        return HierarchicalBitFlipOperator.class;
	    }

	    public XMLSyntaxRule[] getSyntaxRules() {
	        return rules;
	    }

	    private final XMLSyntaxRule[] rules = {
	            AttributeRule.newDoubleRule(MCMCOperator.WEIGHT),
	            AttributeRule.newBooleanRule(USES_SUM_PRIOR,true),
	            new ElementRule(H_PARAMETER,
	                    new XMLSyntaxRule[]{new ElementRule(Parameter.class)}),
	            new ElementRule(STRATA_PARAMETERS,
	                    new XMLSyntaxRule[]{new ElementRule(Parameter.class, 1, Integer.MAX_VALUE)}),
	    };

	}
	
	
	
	

