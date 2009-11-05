package dr.evomodel.substmodel;

import dr.evolution.datatype.DataType;
import dr.evoxml.DataTypeUtils;
import dr.inference.model.Model;
import dr.inference.model.BayesianStochasticSearchVariableSelection;
import dr.inference.distribution.GeneralizedLinearModel;
import dr.inference.distribution.LogLinearModel;
import dr.inference.loggers.LogColumn;
import dr.xml.*;

/**
 * <b>A irreversible class for any data type where
 * rates come from a log-linear model; allows complex eigenstructures.</b>
 *
 * @author Marc A. Suchard
 * @author Alexei J. Drummond
 */
public class GLMSubstitutionModel extends ComplexSubstitutionModel {

    public static final String GLM_SUBSTITUTION_MODEL = "glmSubstitutionModel";

    public GLMSubstitutionModel(String name, DataType dataType, FrequencyModel rootFreqModel,
                                LogLinearModel glm) {

        super(name, dataType, rootFreqModel, null);
        this.glm = glm;
        addModel(glm);
        testProbabilities = new double[stateCount*stateCount];
            
    }

    public double[] getRates() {
        return glm.getXBeta();
    }


    protected void handleModelChangedEvent(Model model, Object object, int index) {
        if (model == glm) {
            updateMatrix = true;
            fireModelChanged();
        }
        else
            super.handleModelChangedEvent(model,object,index);       
    }

    public LogColumn[] getColumns() {
        return glm.getColumns();
    }

    public double getLogLikelihood() {
        double logL = super.getLogLikelihood();
        if (logL == 0 &&
            BayesianStochasticSearchVariableSelection.Utils.connectedAndWellConditioned(testProbabilities,this)) { // Also check that graph is connected
            return 0;
        }
        return Double.NEGATIVE_INFINITY;
    }   

    public static XMLObjectParser PARSER = new AbstractXMLObjectParser() {

        public String getParserName() {
            return GLM_SUBSTITUTION_MODEL;
        }

        public Object parseXMLObject(XMLObject xo) throws XMLParseException {

            DataType dataType = DataTypeUtils.getDataType(xo);

            if (dataType == null) dataType = (DataType) xo.getChild(DataType.class);

            int rateCount = (dataType.getStateCount() - 1) * dataType.getStateCount();

            LogLinearModel glm = (LogLinearModel) xo.getChild(GeneralizedLinearModel.class);

            int length = glm.getXBeta().length;

            if (length != rateCount) {
                throw new XMLParseException("Rates parameter in " + getParserName() + " element should have " + (rateCount) + " dimensions.  However GLM dimension is " + length);
            }

            XMLObject cxo = xo.getChild(ROOT_FREQUENCIES);
            FrequencyModel rootFreq = (FrequencyModel) cxo.getChild(FrequencyModel.class);

            if (dataType != rootFreq.getDataType()) {
                throw new XMLParseException("Data type of " + getParserName() + " element does not match that of its rootFrequencyModel.");
            }

            return new GLMSubstitutionModel(xo.getId(), dataType, rootFreq, glm);
        }

        //************************************************************************
        // AbstractXMLObjectParser implementation
        //************************************************************************

        public String getParserDescription() {
            return "A general model of sequence substitution for any data type where the rates come from the generalized linear model.";
        }

        public Class getReturnType() {
            return SubstitutionModel.class;
        }

        public XMLSyntaxRule[] getSyntaxRules() {
            return rules;
        }

        private XMLSyntaxRule[] rules = new XMLSyntaxRule[]{
                new XORRule(
                        new StringAttributeRule(DataType.DATA_TYPE, "The type of sequence data",
                                DataType.getRegisteredDataTypeNames(), false),
                        new ElementRule(DataType.class)
                ),
                new ElementRule(ROOT_FREQUENCIES, FrequencyModel.class),
                new ElementRule(GeneralizedLinearModel.class),
        };

    };

    private LogLinearModel glm;
    private double[] testProbabilities;    
}
