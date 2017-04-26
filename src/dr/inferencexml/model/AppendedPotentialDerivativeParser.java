package dr.inferencexml.model;

import dr.inference.distribution.DistributionLikelihood;
import dr.inference.distribution.MultivariateDistributionLikelihood;
import dr.inference.model.CompoundGradientProvider;
import dr.inference.model.GradientProvider;
import dr.inference.model.Likelihood;
import dr.math.distributions.GaussianProcessRandomGenerator;
import dr.util.Attribute;
import dr.xml.*;

import java.util.ArrayList;
import java.util.List;

/**
 * @author Max Tolkoff
 */
public class AppendedPotentialDerivativeParser extends AbstractXMLObjectParser {

    public final static String SUM_DERIVATIVE = "appendedPotentialDerivative";
    public static final String SUM_DERIVATIVE2 = "compoundGradientForLikelihood";

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

        List<GradientProvider> gradList = new ArrayList<GradientProvider>();
        List<Likelihood> likelihoodList = new ArrayList<Likelihood>();


        for (int i = 0; i < xo.getChildCount(); ++i) {
            Object obj = xo.getChild(i);
            GradientProvider grad = null;

            Likelihood likelihood = null;

//            int copies = -1;
            if (obj instanceof DistributionLikelihood) {
                DistributionLikelihood dl = (DistributionLikelihood) obj;
                if (!(dl.getDistribution() instanceof GradientProvider)) {
                    throw new XMLParseException("Not a gradient provider");
                }

                likelihood = dl;
                grad = (GradientProvider) dl.getDistribution();

//                copies = 0;
//
//                for (Attribute<double[]> datum : dl.getDataList()) {
//                    Double draw = (Double) gp.nextRandom();
//                    System.err.println("DL: " + datum.getAttributeName() + " " + datum.getAttributeValue().length + " " + "1");
//                    copies += datum.getAttributeValue().length;
//                }
            } else if (obj instanceof MultivariateDistributionLikelihood) {
                MultivariateDistributionLikelihood mdl = (MultivariateDistributionLikelihood) obj;
                if (!(mdl.getDistribution() instanceof GradientProvider)) {
                    throw new XMLParseException("Not a gradient provider");
                }

                likelihood = mdl;
                grad = (GradientProvider) mdl.getDistribution();

//                copies = 0;
//                double[] draw = (double[]) gp.nextRandom();
//                for (Attribute<double[]> datum : mdl.getDataList()) {
////                    System.err.println("ML: " + datum.getAttributeName() + " " + datum.getAttributeValue().length + " " + draw.length);
//                    copies += datum.getAttributeValue().length / draw.length;
//                }
            } else if (obj instanceof GradientProvider) {
                grad = (GradientProvider) obj;
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
//            derivativeList.add((GradientProvider) xo.getChild(i));
//        }


        return new CompoundGradientProvider(gradList, likelihoodList);
    }

    @Override
    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }

    private final XMLSyntaxRule[] rules = {
            new ElementRule(GradientProvider.class, 1, Integer.MAX_VALUE),
    };

    @Override
    public String getParserDescription() {
        return null;
    }

    @Override
    public Class getReturnType() {
        return CompoundGradientProvider.class;
    }
}
