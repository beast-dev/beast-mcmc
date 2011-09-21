package dr.inference.model;

import dr.inference.loggers.LogColumn;
import dr.math.LogTricks;
import dr.xml.*;

import java.util.ArrayList;
import java.util.List;

/**
 * @author Marc A. Suchard
 * @author Andrew Rambaut
 */
public class WeightedMixtureModel extends AbstractModelLikelihood {

    public static final String MIXTURE_MODEL = "mixtureModel";
//    public static final String MIXTURE_WEIGHTS = "weights";
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
        double logSum = Math.log(mixtureWeights.getParameterValue(0)) + likelihoodList.get(0).getLogLikelihood();
        for (int i = 1; i < likelihoodList.size(); ++i) {
            logSum = LogTricks.logSum(logSum,
                    Math.log(mixtureWeights.getParameterValue(i)) + likelihoodList.get(i).getLogLikelihood());
        }   
        return logSum;                
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
            return "This element represents a finite mixture of likelihood models.";
        }

        public Class getReturnType() {
            return CompoundModel.class;
        }

        public XMLSyntaxRule[] getSyntaxRules() {
            return rules;
        }

        private final XMLSyntaxRule[] rules = {
                AttributeRule.newBooleanRule(NORMALIZE, true),
                new ElementRule(Likelihood.class,2,Integer.MAX_VALUE),
                new ElementRule(Parameter.class)
        };
    };


    private final Parameter mixtureWeights;
    List<Likelihood> likelihoodList;


    public static void main(String[] args) {

        final double l1 = -10;
        final double l2 = -2;

        Likelihood like1 = new Likelihood() {


            public Model getModel() {
                return null;
            }

            public double getLogLikelihood() {
                return l1;
            }

            public void makeDirty() {
            }

            public String prettyName() {
                return null;
            }

            public boolean isUsed() {
                return false;
            }

            public void setUsed() {
            }

            public LogColumn[] getColumns() {
                return new LogColumn[0];
            }

            public String getId() {
                return null;
            }

            public void setId(String id) {
            }
        };

        Likelihood like2 = new Likelihood() {

            public Model getModel() {
                return null;
            }

            public double getLogLikelihood() {
                return l2;
            }

            public void makeDirty() {
            }

            public String prettyName() {
                return null;
            }

            public boolean isUsed() {
                return false;
            }

            public void setUsed() {
            }

            public LogColumn[] getColumns() {
                return new LogColumn[0];
            }

            public String getId() {
                return null;
            }

            public void setId(String id) {                
            }
        };

        List<Likelihood> likelihoodList = new ArrayList<Likelihood>();
        likelihoodList.add(like1);
        likelihoodList.add(like2);

        Parameter weights = new Parameter.Default(2);
        double p1 = 0.05;
        weights.setParameterValue(0, p1);
        weights.setParameterValue(1, 1.0 - p1);

        WeightedMixtureModel mixture = new WeightedMixtureModel(likelihoodList, weights);
        System.err.println("getLogLikelihood() = " + mixture.getLogLikelihood());

        double test = Math.log(p1 * Math.exp(l1) + (1.0 - p1) * Math.exp(l2));
        System.err.println("correct            = " + test);
    }

}
