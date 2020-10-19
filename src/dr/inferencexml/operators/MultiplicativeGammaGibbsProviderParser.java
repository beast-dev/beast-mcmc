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
    private static final String ROW_MULTS = "rowMultipliers";
    private static final String ROW = "row";


    @Override
    public Object parseXMLObject(XMLObject xo) throws XMLParseException {
        XMLObject mxo = xo.getChild(ROW_MULTS);
        List<Parameter> rowMultList = new ArrayList<Parameter>();
        for (Object child : mxo.getChildren()) {
            if (child instanceof Parameter) {
                if (((Parameter) child).getDimension() != 1) {
                    throw new XMLParseException("All parameters in  the " + ROW_MULTS +
                            " XML object must have dimension one.");
                }
                rowMultList.add((Parameter) child);
            }
        }

        MultiplicativeGammaGibbsHelper helper =
                (MultiplicativeGammaGibbsHelper) xo.getChild(MultiplicativeGammaGibbsHelper.class);

//        MatrixParameterInterface matParam = (MatrixParameterInterface) xo.getChild(MatrixParameterInterface.class);
        int row = xo.getIntegerAttribute(ROW);

        int k = rowMultList.size();
//        if (helper.getDimension() != k) {
//            throw new XMLParseException("Dimension mismatch: the `" + ROW_MULTS + "` element has dimension " + k +
//                    ", while the `" + MatrixShrinkageLikelihoodParser.MATRIX_SHRINKAGE + "` element has dimension " +
//                    helper.getDimension() + ".");
//
//        }
        if (helper.getColumnDimension() != k) {
            throw new XMLParseException("Dimension mismatch: the `" + ROW_MULTS + "` element has dimension " + k +
                    ", while the matrix parameter element has " +
                    helper.getColumnDimension() + " columns.");
        }

        Parameter[] rowParams = new Parameter[k];
        for (int i = 0; i < k; i++) {
            rowParams[i] = rowMultList.get(i);
        }
        CompoundParameter multParam = new CompoundParameter(ROW_MULTS, rowParams);

        return new GammaGibbsProvider.MultiplicativeGammaGibbsProvider(multParam, helper, row - 1);
    }

    @Override
    public XMLSyntaxRule[] getSyntaxRules() {
        return new XMLSyntaxRule[]{
                new ElementRule(ROW_MULTS, new XMLSyntaxRule[]{
                        new ElementRule(Parameter.class, 1, Integer.MAX_VALUE)
                }),
                new ElementRule(MultiplicativeGammaGibbsHelper.class),
//                new ElementRule(MatrixParameterInterface.class),
                AttributeRule.newIntegerRule(ROW)
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
