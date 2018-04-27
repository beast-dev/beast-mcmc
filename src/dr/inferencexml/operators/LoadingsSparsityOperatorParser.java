package dr.inferencexml.operators;

import dr.inference.model.MatrixParameterInterface;
import dr.inference.operators.LoadingsGibbsTruncatedOperator;
import dr.inference.operators.LoadingsSparsityOperator;
import dr.xml.*;

public class LoadingsSparsityOperatorParser extends AbstractXMLObjectParser {
    public static final String LOADINGS_SPARSITY_OPERATOR = "loadingsSparsityOperator";
    public static final String WEIGHT = "weight";

    @Override
    public Object parseXMLObject(XMLObject xo) throws XMLParseException {
        LoadingsGibbsTruncatedOperator operator = (LoadingsGibbsTruncatedOperator) xo.getChild(LoadingsGibbsTruncatedOperator.class);
        MatrixParameterInterface sparseness = (MatrixParameterInterface) xo.getChild(MatrixParameterInterface.class);
        double weight = xo.getDoubleAttribute(WEIGHT);

        return new LoadingsSparsityOperator(weight, operator, sparseness);
    }

    @Override
    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }

    XMLSyntaxRule[] rules = {
            new ElementRule(LoadingsGibbsTruncatedOperator.class),
            new ElementRule(MatrixParameterInterface.class),
    };

    @Override
    public String getParserDescription() {
        return "Sparseness operator for the sparseness on the loadings of a latent factor model";
    }

    @Override
    public Class getReturnType() {
        return LoadingsSparsityOperator.class;
    }

    @Override
    public String getParserName() {
        return LOADINGS_SPARSITY_OPERATOR;
    }
}
