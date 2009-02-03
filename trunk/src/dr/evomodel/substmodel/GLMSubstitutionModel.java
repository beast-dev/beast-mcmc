package dr.evomodel.substmodel;

import dr.evolution.datatype.DataType;
import dr.evoxml.DataTypeUtils;
import dr.inference.model.Parameter;
import dr.inference.model.Model;
import dr.inference.distribution.GeneralizedLinearModel;
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
                                GeneralizedLinearModel glm) {

        super(name, dataType, rootFreqModel, null);
        this.glm = glm;
        addModel(glm);
            
    }

    public double[] getRates() {
//        double[] xBeta =  glm.getXBeta();
//        for(int i=0; i<xBeta.length; i++)
//            xBeta[i] = Math.exp(xBeta[i]);
//        return xBeta;
        return glm.getXBeta();
    }


    protected void handleModelChangedEvent(Model model, Object object, int index) {
//        System.err.println("a model changed!");
        if (model == glm) {
//            System.err.println("and is glm");
            updateMatrix = true;
            fireModelChanged();
        }
        else
            super.handleModelChangedEvent(model,object,index);       
    }


    public static XMLObjectParser PARSER = new AbstractXMLObjectParser() {

        public String getParserName() {
            return GLM_SUBSTITUTION_MODEL;
        }

        public Object parseXMLObject(XMLObject xo) throws XMLParseException {

            DataType dataType = DataTypeUtils.getDataType(xo);

            if (dataType == null) dataType = (DataType) xo.getChild(DataType.class);

//            XMLObject cxo = (XMLObject) xo.getChild(RATES);
//
//            Parameter ratesParameter = (Parameter) cxo.getChild(Parameter.class);
//
            int rateCount = (dataType.getStateCount() - 1) * dataType.getStateCount();

            GeneralizedLinearModel glm = (GeneralizedLinearModel) xo.getChild(GeneralizedLinearModel.class);

            int length = glm.getXBeta().length;

            if (length != rateCount) {
                throw new XMLParseException("Rates parameter in " + getParserName() + " element should have " + (rateCount) + " dimensions.  However GLM dimension is " + length);
            }

            XMLObject cxo = (XMLObject) xo.getChild(ROOT_FREQUENCIES);
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

    private GeneralizedLinearModel glm;

}
