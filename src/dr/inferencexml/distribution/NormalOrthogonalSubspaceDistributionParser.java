package dr.inferencexml.distribution;

import dr.inference.distribution.DistributionLikelihood;
import dr.inference.distribution.NormalOrthogonalSubspaceDistribution;
import dr.inference.distribution.NormalStatisticsHelpers.IndependentNormalStatisticsProvider;
import dr.inference.model.MatrixParameterInterface;
import dr.math.distributions.Distribution;
import dr.xml.*;

public class NormalOrthogonalSubspaceDistributionParser extends AbstractXMLObjectParser {

    public static final String PARSER = "normalOrthogonalSubspaceDistribution";


    @Override
    public Object parseXMLObject(XMLObject xo) throws XMLParseException {
        IndependentNormalStatisticsProvider provider =
                (IndependentNormalStatisticsProvider) xo.getChild(IndependentNormalStatisticsProvider.class);
        if (provider == null) {
            DistributionLikelihood like = (DistributionLikelihood) xo.getChild(DistributionLikelihood.class);
            Distribution dist = like.getDistribution();
            if (dist instanceof IndependentNormalStatisticsProvider) {
                provider = (IndependentNormalStatisticsProvider) dist;
            } else {
                throw new XMLParseException("Must supply a normal distribution.");
            }
        }
        MatrixParameterInterface param = (MatrixParameterInterface) xo.getChild(MatrixParameterInterface.class);

        for (int i = 0; i < param.getDimension(); i++) {
            if (provider.getNormalMean(i) != 0) {
                throw new XMLParseException("Not implemented for prior with non-zero mean.");
            }
        }

        return new NormalOrthogonalSubspaceDistribution(param.getParameterName() + "." + PARSER, provider, param);
    }

    @Override
    public XMLSyntaxRule[] getSyntaxRules() {
        return new XMLSyntaxRule[]{
                new XORRule(new XMLSyntaxRule[]{
                        new ElementRule(IndependentNormalStatisticsProvider.class),
                        new ElementRule(DistributionLikelihood.class)
                }),

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
