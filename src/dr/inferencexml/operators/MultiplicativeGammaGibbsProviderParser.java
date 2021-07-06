package dr.inferencexml.operators;

import dr.evomodel.continuous.MatrixShrinkageLikelihood;
import dr.inference.model.CompoundParameter;
import dr.inference.model.MatrixParameterInterface;
import dr.inference.model.Parameter;
import dr.inference.operators.repeatedMeasures.GammaGibbsProvider;
import dr.inference.operators.repeatedMeasures.MultiplicativeGammaGibbsHelper;
import dr.xml.*;

import java.util.ArrayList;
import java.util.List;

/**
 * @author Gabriel Hassler
 * @author Marc A. Suchard
 */

public class MultiplicativeGammaGibbsProviderParser extends AbstractXMLObjectParser {

    private static final String MULTIPLICATIVE_PROVIDER = "multiplicativeGammaGibbsProvider";


    @Override
    public Object parseXMLObject(XMLObject xo) throws XMLParseException {
        Parameter multParam = (Parameter) xo.getChild(Parameter.class);

        MultiplicativeGammaGibbsHelper helper =
                (MultiplicativeGammaGibbsHelper) xo.getChild(MultiplicativeGammaGibbsHelper.class);


        int k = multParam.getDimension();

        if (helper.getColumnDimension() != k) {
            throw new XMLParseException("Dimension mismatch: the parameter with id `" + multParam.getId() + "` has dimension " + k +
                    ", while the helper has " + helper.getColumnDimension() + " columns.");
        }


        return new GammaGibbsProvider.MultiplicativeGammaGibbsProvider(multParam, helper);
    }

    @Override
    public XMLSyntaxRule[] getSyntaxRules() {
        return new XMLSyntaxRule[]{
                new ElementRule(Parameter.class),
                new ElementRule(MultiplicativeGammaGibbsHelper.class),
        };
    }

    @Override
    public String getParserDescription() {
        return "Returns sufficient statistics for multiplicative gamma distribution.";
    }

    @Override
    public Class getReturnType() {
        return GammaGibbsProvider.MultiplicativeGammaGibbsProvider.class;
    }

    @Override
    public String getParserName() {
        return MULTIPLICATIVE_PROVIDER;
    }
}
