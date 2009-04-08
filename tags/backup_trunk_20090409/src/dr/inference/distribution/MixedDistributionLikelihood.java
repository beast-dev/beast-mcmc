package dr.inference.distribution;

import dr.inference.model.Statistic;
import dr.inference.model.Likelihood;
import dr.inference.model.Model;
import dr.inference.model.CompoundModel;
import dr.xml.*;
import org.w3c.dom.Element;
import org.w3c.dom.Document;

/**
 * A class that returns the log likelihood of a multidimensional statistic
 * being distributed according to a mixture of parametric distributions, where
 * the distribution membership of each element is determined by a separate vector
 * of indices.
 * @author Alexei Drummond
 * @version $Id: DistributionLikelihood.java,v 1.11 2005/05/25 09:35:28 rambaut Exp $
 */

public class MixedDistributionLikelihood extends Likelihood.Abstract {

    public static final String DISTRIBUTION_LIKELIHOOD = "mixedDistributionLikelihood";

    public static final String DISTRIBUTION0 = "distribution0";
    public static final String DISTRIBUTION1 = "distribution1";
    public static final String DATA = "data";
    public static final String INDICATORS = "indicators";

    public MixedDistributionLikelihood(ParametricDistributionModel[] distributions, Statistic data, Statistic indicators) {
        // Mixed Distribution Likelihood contains two distribution models, not necessarily constant.
        // To cater for that, they need to be returned as the "Model" of this likelyhood so that their state is correctly
        // restored when an operation involving their parameters fails.

        super(new CompoundModel("MixedeDistributions"));

        final CompoundModel cm = (CompoundModel)this.getModel();
        for( ParametricDistributionModel m : distributions ) {
          cm.addModel(m);
        }

        this.distributions = distributions;
        this.data = data;
        this.indicators = indicators;

        if (indicators.getDimension() == data.getDimension()-1) {
            impliedOne = true;
        } else if (indicators.getDimension() != data.getDimension()) {
            throw new IllegalArgumentException();
        }
    }

    // **************************************************************
    // Likelihood IMPLEMENTATION
    // **************************************************************

    /**
     * Calculate the log likelihood of the current state.
     * @return the log likelihood.
     */
    public final double calculateLogLikelihood() {

        double logL = 0.0;

        for (int j = 0; j < data.getDimension(); j++) {

            int index;
            if (impliedOne) {
                if (j == 0) {
                    index = 1;
                } else {
                    index = (int)indicators.getStatisticValue(j-1);
                }
            } else {
                index = (int)indicators.getStatisticValue(j);
            }

            logL += distributions[index].logPdf(data.getStatisticValue(j));
        }
        //System.err.println("mixed: " + logL);
        return logL;
    }

    /**
	 * Overridden to always return false.
	 */
	protected boolean getLikelihoodKnown() {
		return false;
	}

    // **************************************************************
    // XMLElement IMPLEMENTATION
    // **************************************************************

    public Element createElement(Document d) {
        throw new RuntimeException("Not implemented yet!");
    }


    /**
     * Reads a distribution likelihood from a DOM Document element.
     */
    public static XMLObjectParser PARSER = new AbstractXMLObjectParser() {

        public String getParserName() { return MixedDistributionLikelihood.DISTRIBUTION_LIKELIHOOD; }

        public Object parseXMLObject(XMLObject xo) {

            XMLObject cxo0 = (XMLObject)xo.getChild(DISTRIBUTION0);
            ParametricDistributionModel model0 = (ParametricDistributionModel)cxo0.getChild(ParametricDistributionModel.class);

            XMLObject cxo1 = (XMLObject)xo.getChild(DISTRIBUTION1);
            ParametricDistributionModel model1 = (ParametricDistributionModel)cxo1.getChild(ParametricDistributionModel.class);

            Statistic data = (Statistic)((XMLObject)xo.getChild(DATA)).getChild(Statistic.class);
            Statistic indicators = (Statistic)((XMLObject)xo.getChild(INDICATORS)).getChild(Statistic.class);

            ParametricDistributionModel[] models = {model0, model1};

            return new MixedDistributionLikelihood(models, data, indicators);
        }

        //************************************************************************
        // AbstractXMLObjectParser implementation
        //************************************************************************

        public XMLSyntaxRule[] getSyntaxRules() { return rules; }

        private final XMLSyntaxRule[] rules = {
            new ElementRule(DISTRIBUTION0,
                new XMLSyntaxRule[] { new ElementRule(ParametricDistributionModel.class) }),
            new ElementRule(DISTRIBUTION1,
                new XMLSyntaxRule[] { new ElementRule(ParametricDistributionModel.class) }),
            new ElementRule(DATA, new XMLSyntaxRule[] {new ElementRule(Statistic.class)}),
            new ElementRule(INDICATORS, new XMLSyntaxRule[] {new ElementRule(Statistic.class)}),
        };

        public String getParserDescription() {
            return "Calculates the likelihood of some data given some mix of parametric distributions.";
        }

        public Class getReturnType() { return Likelihood.class; }
    };

    private final ParametricDistributionModel[] distributions;
    private final Statistic data;
    private final Statistic indicators;
    private boolean impliedOne = false;

    public Model[] getUniqueModels() {
        Model[] m = new Model[distributions[0] == distributions[1] ? 1 : 2];
        m[0] = distributions[0];
        if( m.length > 1 ) m[1] = distributions[1];
        return m;
    }
}
