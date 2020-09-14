package dr.inferencexml.distribution;

import dr.inference.distribution.NormalOrthogonalSubspaceDistribution;
import dr.inference.distribution.NormalStatisticsHelpers.IndependentNormalStatisticsProvider;
import dr.inference.model.MatrixParameterInterface;
import dr.xml.*;

public class NormalOrthogonalSubspaceDistributionParser extends AbstractXMLObjectParser {

    private static final String PARSER = "normalOrthogonalSubspaceDistribution";


    @Override
    public Object parseXMLObject(XMLObject xo) throws XMLParseException {
        IndependentNormalStatisticsProvider provider = (IndependentNormalStatisticsProvider)
                xo.getChild(IndependentNormalStatisticsProvider.class);
        MatrixParameterInterface param = (MatrixParameterInterface) xo.getChild(MatrixParameterInterface.class);

        return new NormalOrthogonalSubspaceDistribution(param.getParameterName() + "." + PARSER, provider, param);
    }

    @Override
    public XMLSyntaxRule[] getSyntaxRules() {
        return new XMLSyntaxRule[]{
                new ElementRule(IndependentNormalStatisticsProvider.class),
                new ElementRule(MatrixParameterInterface.class)
        };
    }

    @Override
    public String getParserDescription() {
        return "Applies matrix von Mises-Fisher distribution to the orthogonal subspace of a matrix.";
    }

    @Override
    public Class getReturnType() {
        return NormalOrthogonalSubspaceDistribution.class;
    }

    @Override
    public String getParserName() {
        return PARSER;
    }
}
