package dr.evomodel.substmodel;

import dr.evolution.datatype.HiddenNucleotides;
import dr.inference.model.Parameter;
import dr.xml.*;

/**
 * A general time reversible model of nucleotide evolution with
 * covarion hidden rate categories.
 *
 * @author Alexei Drummond
 * @version $Id$
 */
public class CovarionGTR extends AbstractCovarionDNAModel {

    public static final String GTR_COVARION_MODEL = "gtrCovarionModel";

    /**
     * constructor
     *
     * @param dataType        the data type
     * @param freqModel       the equlibrium frequencies
     * @param hiddenRates     the relative rates of the hidden categories
     *                        (first hidden category has rate 1.0 so this parameter
     *                        has dimension one less than number of hidden categories.
     * @param rateACParameter the relative rate of A<->C substitutions
     * @param rateAGParameter the relative rate of A<->G substitutions
     * @param rateATParameter the relative rate of A<->T substitutions
     * @param rateCGParameter the relative rate of C<->G substitutions
     * @param rateCTParameter the relative rate of C<->T substitutions
     * @param rateGTParameter the relative rate of G<->T substitutions
     *                        each hidden category.
     * @param switchingRates  rates of switching between hidden rate classes
     */
    public CovarionGTR(
            HiddenNucleotides dataType,
            Parameter hiddenRates,
            Parameter switchingRates,
            Parameter rateACParameter,
            Parameter rateAGParameter,
            Parameter rateATParameter,
            Parameter rateCGParameter,
            Parameter rateCTParameter,
            Parameter rateGTParameter,
            FrequencyModel freqModel) {

        super(GTR_COVARION_MODEL, dataType, hiddenRates, switchingRates, freqModel);
        if (rateACParameter != null) {
            addParameter(rateACParameter);
            rateACParameter.addBounds(new Parameter.DefaultBounds(Double.POSITIVE_INFINITY, 0.0, 1));
            this.rateACParameter = rateACParameter;
        }

        if (rateAGParameter != null) {
            addParameter(rateAGParameter);
            rateAGParameter.addBounds(new Parameter.DefaultBounds(Double.POSITIVE_INFINITY, 0.0, 1));
            this.rateAGParameter = rateAGParameter;
        }

        if (rateATParameter != null) {
            addParameter(rateATParameter);
            rateATParameter.addBounds(new Parameter.DefaultBounds(Double.POSITIVE_INFINITY, 0.0, 1));
            this.rateATParameter = rateATParameter;
        }

        if (rateCGParameter != null) {
            addParameter(rateCGParameter);
            rateCGParameter.addBounds(new Parameter.DefaultBounds(Double.POSITIVE_INFINITY, 0.0, 1));
            this.rateCGParameter = rateCGParameter;
        }

        if (rateCTParameter != null) {
            addParameter(rateCTParameter);
            rateCTParameter.addBounds(new Parameter.DefaultBounds(Double.POSITIVE_INFINITY, 0.0, 1));
            this.rateCTParameter = rateCTParameter;
        }

        if (rateGTParameter != null) {
            addParameter(rateGTParameter);
            rateGTParameter.addBounds(new Parameter.DefaultBounds(Double.POSITIVE_INFINITY, 0.0, 1));
            this.rateGTParameter = rateGTParameter;
        }
    }

    double[] getRelativeDNARates() {
        double[] rr = new double[6];
        for (int i = 0; i < rr.length; i++) {
            rr[i] = 1.0;
        }

        if (rateACParameter != null) {
            rr[0] = rateACParameter.getParameterValue(0);
        }
        if (rateAGParameter != null) {
            rr[1] = rateAGParameter.getParameterValue(0);
        }
        if (rateATParameter != null) {
            rr[2] = rateATParameter.getParameterValue(0);
        }
        if (rateCGParameter != null) {
            rr[3] = rateCGParameter.getParameterValue(0);
        }
        if (rateCTParameter != null) {
            rr[4] = rateCTParameter.getParameterValue(0);
        }
        if (rateGTParameter != null) {
            rr[5] = rateGTParameter.getParameterValue(0);
        }
        return rr;
    }

    /**
     * Parses an element from an DOM document into a TwoStateCovarionModel
     */
    public static XMLObjectParser PARSER = new AbstractXMLObjectParser() {

        public String getParserName() {
            return GTR_COVARION_MODEL;
        }

        public Object parseXMLObject(XMLObject xo) throws XMLParseException {

            XMLObject cxo = (XMLObject) xo.getChild(FREQUENCIES);
            FrequencyModel freqModel = (FrequencyModel) cxo.getChild(FrequencyModel.class);

            HiddenNucleotides dataType = (HiddenNucleotides) freqModel.getDataType();

            Parameter hiddenRates = (Parameter) xo.getSocketChild(HIDDEN_CLASS_RATES);
            Parameter switchingRates = (Parameter) xo.getSocketChild(SWITCHING_RATES);

            Parameter rateACParameter = null;
            if (xo.hasSocket(GTR.A_TO_C)) {
                rateACParameter = (Parameter) xo.getSocketChild(GTR.A_TO_C);
            }
            Parameter rateAGParameter = null;
            if (xo.hasSocket(GTR.A_TO_G)) {
                rateAGParameter = (Parameter) xo.getSocketChild(GTR.A_TO_G);
            }
            Parameter rateATParameter = null;
            if (xo.hasSocket(GTR.A_TO_T)) {
                rateATParameter = (Parameter) xo.getSocketChild(GTR.A_TO_T);
            }
            Parameter rateCGParameter = null;
            if (xo.hasSocket(GTR.C_TO_G)) {
                rateCGParameter = (Parameter) xo.getSocketChild(GTR.C_TO_G);
            }
            Parameter rateCTParameter = null;
            if (xo.hasSocket(GTR.C_TO_T)) {
                rateCTParameter = (Parameter) xo.getSocketChild(GTR.C_TO_T);
            }
            Parameter rateGTParameter = null;
            if (xo.hasSocket(GTR.G_TO_T)) {
                rateGTParameter = (Parameter) xo.getSocketChild(GTR.G_TO_T);
            }


            if (dataType != freqModel.getDataType()) {
                throw new XMLParseException("Data type of " + getParserName() + " element does not match that of its frequencyModel.");
            }

            return new CovarionGTR(
                    dataType,
                    hiddenRates,
                    switchingRates,
                    rateACParameter,
                    rateAGParameter,
                    rateATParameter,
                    rateCGParameter,
                    rateCTParameter,
                    rateGTParameter,
                    freqModel);
        }

        //************************************************************************
        // AbstractXMLObjectParser implementation
        //************************************************************************

        public String getParserDescription() {
            return "A covarion substitution model of langauge evolution with binary data and a hidden rate state with two rates.";
        }

        public Class getReturnType() {
            return CovarionGTR.class;
        }

        public XMLSyntaxRule[] getSyntaxRules() {
            return rules;
        }

        private XMLSyntaxRule[] rules = new XMLSyntaxRule[]{
                new ElementRule(GTR.A_TO_C, Parameter.class, "relative rate of A<->C substitution", true),
                new ElementRule(GTR.A_TO_G, Parameter.class, "relative rate of A<->G substitution", true),
                new ElementRule(GTR.A_TO_T, Parameter.class, "relative rate of A<->T substitution", true),
                new ElementRule(GTR.C_TO_G, Parameter.class, "relative rate of C<->G substitution", true),
                new ElementRule(GTR.C_TO_T, Parameter.class, "relative rate of C<->T substitution", true),
                new ElementRule(GTR.G_TO_T, Parameter.class, "relative rate of G<->T substitution", true),
                new ElementRule(HIDDEN_CLASS_RATES, Parameter.class),
                new ElementRule(SWITCHING_RATES, Parameter.class),
                new ElementRule(FREQUENCIES, FrequencyModel.class),
        };

    };

    private Parameter rateACParameter = null;
    private Parameter rateAGParameter = null;
    private Parameter rateATParameter = null;
    private Parameter rateCGParameter = null;
    private Parameter rateCTParameter = null;
    private Parameter rateGTParameter = null;
}
