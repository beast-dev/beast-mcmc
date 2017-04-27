package dr.inferencexml.model;

import dr.inference.distribution.DistributionLikelihood;
import dr.inference.distribution.MultivariateDistributionLikelihood;
import dr.inference.model.*;
import dr.xml.*;

import java.util.ArrayList;
import java.util.List;

/**
 * @author Max Tolkoff
 */
public class AppendedPotentialDerivativeParser extends AbstractXMLObjectParser {

    public final static String SUM_DERIVATIVE = "appendedPotentialDerivative";
    public static final String SUM_DERIVATIVE2 = "compoundGradient";

    @Override
    public String getParserName() {
        return SUM_DERIVATIVE;
    }

    @Override
    public String[] getParserNames() {
        return new String[] { SUM_DERIVATIVE, SUM_DERIVATIVE2 };
    }

    @Override
    public Object parseXMLObject(XMLObject xo) throws XMLParseException {

        List<GradientWrtParameterProvider> gradList = new ArrayList<GradientWrtParameterProvider>();
        List<Likelihood> likelihoodList = new ArrayList<Likelihood>();


        for (int i = 0; i < xo.getChildCount(); ++i) {

            Object obj = xo.getChild(i);

            GradientWrtParameterProvider grad;
            Likelihood likelihood;

            if (obj instanceof DistributionLikelihood) {
                DistributionLikelihood dl = (DistributionLikelihood) obj;
                if (!(dl.getDistribution() instanceof GradientProvider)) {
                    throw new XMLParseException("Not a gradient provider");
                }

                throw new RuntimeException("Not yet implemented");

            } else if (obj instanceof MultivariateDistributionLikelihood) {
                final MultivariateDistributionLikelihood mdl = (MultivariateDistributionLikelihood) obj;
                if (!(mdl.getDistribution() instanceof GradientProvider)) {
                    throw new XMLParseException("Not a gradient provider");
                }

                final GradientProvider provider = (GradientProvider) mdl.getDistribution();
                final Parameter parameter = mdl.getDataParameter();
                likelihood = mdl;

                grad = new GradientWrtParameterProvider() { // Return gradient w.r.t. parameter

                    @Override
                    public Likelihood getLikelihood() {
                        return mdl;
                    }

                    @Override
                    public Parameter getParameter() {
                        return parameter;
                    }

                    @Override
                    public int getDimension() {
                        return parameter.getDimension();
                    }

                    @Override
                    public double[] getGradientLogDensity() {
                        return provider.getGradientLogDensity(parameter.getParameterValues());
                    }
                };

            } else if (obj instanceof GradientWrtParameterProvider) {
                grad = (GradientWrtParameterProvider) obj;
                likelihood = grad.getLikelihood();
            } else {
                throw new XMLParseException("Not a Gaussian process");
            }

            gradList.add(grad);
            likelihoodList.add(likelihood);
        }

        return new CompoundGradient(gradList);
    }

    @Override
    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }

    private final XMLSyntaxRule[] rules = {
            new ElementRule(GradientWrtParameterProvider.class, 1, Integer.MAX_VALUE),
    };

    @Override
    public String getParserDescription() {
        return null;
    }

    @Override
    public Class getReturnType() {
        return CompoundGradient.class;
    }
}
