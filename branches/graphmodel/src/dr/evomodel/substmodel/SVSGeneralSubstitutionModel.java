package dr.evomodel.substmodel;

import dr.evolution.datatype.*;
import dr.inference.model.*;
import dr.inference.loggers.LogColumn;
import dr.inference.loggers.NumberColumn;
import dr.xml.*;

/**
 * <b>A general model of sequence substitution with stochastic variable selection</b>. A general reversible class for any
 * data type.
 *
 * @author Marc Suchard
 * @version $Id: SVSGeneralSubstitutionModel.java,v 1.37 2006/05/05 03:05:10 msuchard Exp $
 */

public class SVSGeneralSubstitutionModel extends GeneralSubstitutionModel implements Likelihood,
        BayesianStochasticSearchVariableSelection {

    public static final String SVS_GENERAL_SUBSTITUTION_MODEL = "svsGeneralSubstitutionModel";
    public static final String INDICATOR = "rateIndicator";
    public static final String ROOT_FREQ = "rootFrequencies";


    public SVSGeneralSubstitutionModel(DataType dataType, FrequencyModel freqModel, Parameter parameter,
                                       Parameter indicator) { //, int relativeTo) {
        super(dataType, freqModel, parameter, 1);

        rateIndicator = indicator;
        addVariable(rateIndicator);
    }

    protected SVSGeneralSubstitutionModel(String name, DataType dataType, FrequencyModel freqModel, int relativeTo) {
        super(name, dataType, freqModel, relativeTo);
    }

    public Parameter getIndicators() {
        return rateIndicator;
    }

    public boolean validState() {
        return !updateMatrix || Utils.connectedAndWellConditioned(probability,this);
    }

    protected void handleVariableChangedEvent(Variable variable, int index, Parameter.ChangeType type) {
        if (variable == ratesParameter && rateIndicator.getParameterValue(index) == 0)
            return; // Does not affect likelihood
        super.handleVariableChangedEvent(variable,index,type);
    }

    /**
     * Get the model.
     *
     * @return the model.
     */
    public Model getModel() {
        return this;
    }

    /**
     * Get the log likelihood.
     *
     * @return the log likelihood.
     */
    public double getLogLikelihood() {
        if (updateMatrix) {
            if (!Utils.connectedAndWellConditioned(probability,this)) {
                return Double.NEGATIVE_INFINITY;
            }
        }
        return 0;
    }

    /**
     * Forces a complete recalculation of the likelihood next time getLikelihood is called
     */
    public void makeDirty() {
       updateMatrix = true;
    }

    /**
     * @return A detailed name of likelihood for debugging.
     */
    public String prettyName() {
        return "SVSGeneralSubstitutionModel-connectedness";
    }

      // **************************************************************
    // Loggable IMPLEMENTATION
    // **************************************************************

    public LogColumn[] getColumns() {
        return new LogColumn[]{
                new LikelihoodColumn(getId())
        };
    }

    protected class LikelihoodColumn extends NumberColumn {
        public LikelihoodColumn(String label) {
            super(label);
        }

        public double getDoubleValue() {
            return getLogLikelihood();
        }
    }

    private double[] probability = null;

    protected void setupRelativeRates() {

        for (int i = 0; i < relativeRates.length; i++) {
            relativeRates[i] = ratesParameter.getParameterValue(i) * rateIndicator.getParameterValue(i);
        }
    }

    void normalize(double[][] matrix, double[] pi) {
        double subst = 0.0;
        int dimension = pi.length;

        //final int dim = rateIndicator.getDimension();
        //int sum = 0;
        //for (int i = 0; i < dim; i++)
        //	sum += rateIndicator.getParameterValue(i);


        for (int i = 0; i < dimension; i++)
            subst += -matrix[i][i] * pi[i];

        for (int i = 0; i < dimension; i++) {
            for (int j = 0; j < dimension; j++) {
                matrix[i][j] = matrix[i][j] / subst; // / sum;
            }
        }
    }

    public static XMLObjectParser PARSER = new AbstractXMLObjectParser() {

        public String[] getParserNames() {
            return new String[] {
                    getParserName(),"beast_"+getParserName()
            };
        }

        public String getParserName() {
            return SVS_GENERAL_SUBSTITUTION_MODEL;
        }

        public Object parseXMLObject(XMLObject xo) throws XMLParseException {

            Parameter ratesParameter;
            Parameter indicatorParameter;

            XMLObject cxo = xo.getChild(FREQUENCIES);
            FrequencyModel freqModel = (FrequencyModel) cxo.getChild(FrequencyModel.class);

            DataType dataType = null;

            if (xo.hasAttribute(DataType.DATA_TYPE)) {
                String dataTypeStr = xo.getStringAttribute(DataType.DATA_TYPE);
                if (dataTypeStr.equals(Nucleotides.DESCRIPTION)) {
                    dataType = Nucleotides.INSTANCE;
                } else if (dataTypeStr.equals(AminoAcids.DESCRIPTION)) {
                    dataType = AminoAcids.INSTANCE;
                } else if (dataTypeStr.equals(Codons.DESCRIPTION)) {
                    dataType = Codons.UNIVERSAL;
                } else if (dataTypeStr.equals(TwoStates.DESCRIPTION)) {
                    dataType = TwoStates.INSTANCE;
                }
            }

            if (dataType == null) dataType = (DataType) xo.getChild(DataType.class);

            cxo = xo.getChild(RATES);

            ratesParameter = (Parameter) cxo.getChild(Parameter.class);

            if (dataType != freqModel.getDataType()) {
                throw new XMLParseException("Data type of " + getParserName() + " element does not match that of its frequencyModel.");
            }

            int rateCount = ((dataType.getStateCount() - 1) * dataType.getStateCount()) / 2;

            if (ratesParameter == null) {

                if (rateCount == 1) {
                    // simplest model for binary traits...
                } else {
                    throw new XMLParseException("No rates parameter found in " + getParserName());
                }
            } else if (ratesParameter.getDimension() != rateCount) {
                throw new XMLParseException("Rates parameter in " + getParserName() + " element should have " + (rateCount) + " dimensions.  However parameter dimension is " + ratesParameter.getDimension());
            }

            cxo = xo.getChild(INDICATOR);

            indicatorParameter = (Parameter) cxo.getChild(Parameter.class);

            if (indicatorParameter == null || ratesParameter == null || indicatorParameter.getDimension() != ratesParameter.getDimension())
                throw new XMLParseException("Rates and indicator parameters in " + getParserName() + " element must be the same dimension.");

            if (xo.hasChildNamed(ROOT_FREQ)) {

                cxo = xo.getChild(ROOT_FREQ);
                FrequencyModel rootFreq = (FrequencyModel) cxo.getChild(FrequencyModel.class);

                if (dataType != rootFreq.getDataType()) {
                    throw new XMLParseException("Data type of " + getParserName() + " element does not match that of its rootFrequencyModel.");
                }

                return new SVSIrreversibleSubstitutionModel(dataType, freqModel, rootFreq, ratesParameter, indicatorParameter);

            }

            return new SVSGeneralSubstitutionModel(dataType, freqModel, ratesParameter, indicatorParameter);
        }

        //************************************************************************
        // AbstractXMLObjectParser implementation
        //************************************************************************

        public String getParserDescription() {
            return "A general reversible model of sequence substitution for any data type with stochastic variable selection.";
        }

        public Class getReturnType() {
            return SubstitutionModel.class;
        }

        public XMLSyntaxRule[] getSyntaxRules() {
            return rules;
        }

        private XMLSyntaxRule[] rules = new XMLSyntaxRule[]{
                new XORRule(
                        new StringAttributeRule(DataType.DATA_TYPE, "The type of sequence data", new String[]{Nucleotides.DESCRIPTION, AminoAcids.DESCRIPTION, Codons.DESCRIPTION, TwoStates.DESCRIPTION}, false),
                        new ElementRule(DataType.class)
                ),
                new ElementRule(FREQUENCIES, FrequencyModel.class),
                new ElementRule(RATES,
                        new XMLSyntaxRule[]{
                                new ElementRule(Parameter.class)}
                ),
                new ElementRule(INDICATOR,
                        new XMLSyntaxRule[]{
                                new ElementRule(Parameter.class)
                        }),
                new ElementRule(ROOT_FREQ,
                        new XMLSyntaxRule[]{
                                new ElementRule(FrequencyModel.class)
                        }, 0, 1)
        };

    };
    
    @Override
    public boolean isUsed() {
        return super.isUsed() && isUsed;
    }

    public void setUsed() {
        isUsed = true;
    }

    private boolean isUsed = false;

    private Parameter rateIndicator;
}
