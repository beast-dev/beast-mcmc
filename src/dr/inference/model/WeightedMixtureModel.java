package dr.inference.model;

import dr.inference.loggers.LogColumn;
import dr.xml.*;

import java.util.ArrayList;
import java.util.List;

/**
 * @author Marc A. Suchard
 */
public class WeightedMixtureModel extends AbstractModelLikelihood {

    public static final String MIXTURE_MODEL = "mixtureModel";
    public static final String MIXTURE_WEIGHTS = "weights";
    public static final String NORMALIZE = "normalize";

    public WeightedMixtureModel(List<Likelihood> likelihoodList, Parameter mixtureWeights) {
        super(MIXTURE_MODEL);
        this.likelihoodList = likelihoodList;
        this.mixtureWeights = mixtureWeights;
        addVariable(mixtureWeights);
    }

    protected void handleModelChangedEvent(Model model, Object object, int index) {
    }

    protected final void handleVariableChangedEvent(Variable variable, int index, Parameter.ChangeType type) {
    }

    protected void storeState() {
    }

    protected void restoreState() {
    }

    protected void acceptState() {
    }

    public Model getModel() {
        return this;
    }

    public double getLogLikelihood() {
        double sum = 0.0;
        for (int i = 0; i < likelihoodList.size(); i++) {
            // todo is there some way to keep everything on the log scale?
            sum += mixtureWeights.getParameterValue(i) * Math.exp(likelihoodList.get(i).getLogLikelihood());
        }
        return Math.log(sum);
    }

    public void makeDirty() {
    }

    public LogColumn[] getColumns() {
        return new LogColumn[0];
    }

    public static XMLObjectParser PARSER = new AbstractXMLObjectParser() {

        public String getParserName() {
            return MIXTURE_MODEL;
        }

        public Object parseXMLObject(XMLObject xo) throws XMLParseException {

            Parameter weights = (Parameter) xo.getChild(Parameter.class);
            List<Likelihood> likelihoodList = new ArrayList<Likelihood>();

            for (int i = 0; i < xo.getChildCount(); i++) {
                if (xo.getChild(i) instanceof Likelihood)
                    likelihoodList.add((Likelihood) xo.getChild(i));
            }

            if (weights.getDimension() != likelihoodList.size()) {
                throw new XMLParseException("Dim of " + weights.getId() + " does not match the number of likelihoods");
            }

            if (xo.hasAttribute(NORMALIZE)) {
                if (xo.getBooleanAttribute(NORMALIZE)) {
                    double sum = 0;
                    for (int i = 0; i < weights.getDimension(); i++)
                        sum += weights.getParameterValue(i);
                    for (int i = 0; i < weights.getDimension(); i++)
                        weights.setParameterValue(i, weights.getParameterValue(i) / sum);
                }
            }

            if (!normalized(weights))
                throw new XMLParseException("Parameter +" + weights.getId() + " must lie on the simplex");

            return new WeightedMixtureModel(likelihoodList, weights);
        }

        private boolean normalized(Parameter p) {
            double sum = 0;
            for (int i = 0; i < p.getDimension(); i++)
                sum += p.getParameterValue(i);
            return (sum == 1.0);
        }

        //************************************************************************
        // AbstractXMLObjectParser implementation
        //************************************************************************

        public String getParserDescription() {
            return "This element represents a discrete mixture of tree likelihood models.";
        }

        public Class getReturnType() {
            return CompoundModel.class;
        }

        public XMLSyntaxRule[] getSyntaxRules() {
            return rules;
        }

        private final XMLSyntaxRule[] rules = {
//			new ElementRule(TreeLikelihood.class,2,2),
//			new ElementRule(Parameter.class)
        };
    };


    private final Parameter mixtureWeights;
    List<Likelihood> likelihoodList;


}
