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
            GradientWrtParameterProvider grad = null;

            Likelihood likelihood = null;

//            int copies = -1;
            if (obj instanceof DistributionLikelihood) {
                DistributionLikelihood dl = (DistributionLikelihood) obj;
                if (!(dl.getDistribution() instanceof GradientProvider)) {
                    throw new XMLParseException("Not a gradient provider");
                }

//                likelihood = dl;
//                grad = (GradientWrtParameterProvider) dl.getDistribution();
                throw new RuntimeException("Not yet implemented");

//                copies = 0;
//
//                for (Attribute<double[]> datum : dl.getDataList()) {
//                    Double draw = (Double) gp.nextRandom();
//                    System.err.println("DL: " + datum.getAttributeName() + " " + datum.getAttributeValue().length + " " + "1");
//                    copies += datum.getAttributeValue().length;
//                }
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

//                grad = (GradientWrtParameterProvider) mdl.getDistribution();

//                copies = 0;
//                double[] draw = (double[]) gp.nextRandom();
//                for (Attribute<double[]> datum : mdl.getDataList()) {
////                    System.err.println("ML: " + datum.getAttributeName() + " " + datum.getAttributeValue().length + " " + draw.length);
//                    copies += datum.getAttributeValue().length / draw.length;
//                }
            } else if (obj instanceof GradientWrtParameterProvider) {
                grad = (GradientWrtParameterProvider) obj;
                likelihood = grad.getLikelihood();

//                copies = 1;
            } else {
                throw new XMLParseException("Not a Gaussian process");
            }
            gradList.add(grad);
            likelihoodList.add(likelihood);
//            copyList.add(copies);
        }

//        for (int i = 0; i < xo.getChildCount(); i++) {
//            derivativeList.add((GradientWrtParameterProvider) xo.getChild(i));
//        }


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
