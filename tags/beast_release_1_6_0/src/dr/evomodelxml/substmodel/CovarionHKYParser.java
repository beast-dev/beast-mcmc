package dr.evomodelxml.substmodel;

import dr.evolution.datatype.HiddenNucleotides;
import dr.evomodel.substmodel.AbstractCovarionDNAModel;
import dr.evomodel.substmodel.CovarionHKY;
import dr.evomodel.substmodel.FrequencyModel;
import dr.inference.model.Parameter;
import dr.xml.*;

/**
 * Parses an element from an DOM document into a DemographicModel. Recognises
 * ConstantPopulation and ExponentialGrowth.
 */
public class CovarionHKYParser extends AbstractXMLObjectParser {
    public static final String COVARION_HKY = "CovarionHKYModel";
    public static final String KAPPA = HKYParser.KAPPA;

    public String getParserName() {
        return COVARION_HKY;
    }

    public String getParserDescription() {
        return "A covarion HKY model.";
    }

    public Class getReturnType() {
        return CovarionHKY.class;
    }

    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }

    private XMLSyntaxRule[] rules = new XMLSyntaxRule[]{
            new ElementRule(KAPPA, new XMLSyntaxRule[]{
                    new ElementRule(Parameter.class,
                            "A parameter representing the transition transversion bias")}),
            new ElementRule(AbstractCovarionDNAModel.SWITCHING_RATES, new XMLSyntaxRule[]{
                    new ElementRule(Parameter.class,
                            "A parameter representing the rate of change between the different classes")}),
            new ElementRule(AbstractCovarionDNAModel.HIDDEN_CLASS_RATES, new XMLSyntaxRule[]{
                    new ElementRule(Parameter.class,
                            "A parameter representing the rates of the hidden classes relative to the first hidden class.")})
    };

    public Object parseXMLObject(XMLObject xo) throws XMLParseException {

        Parameter kappaParam;
        Parameter switchingRates;
        Parameter hiddenClassRates;
        FrequencyModel freqModel;

        kappaParam = (Parameter) xo.getElementFirstChild(KAPPA);
        switchingRates = (Parameter) xo.getElementFirstChild(AbstractCovarionDNAModel.SWITCHING_RATES);
        hiddenClassRates = (Parameter) xo.getElementFirstChild(AbstractCovarionDNAModel.HIDDEN_CLASS_RATES);
        freqModel = (FrequencyModel) xo.getElementFirstChild(AbstractCovarionDNAModel.FREQUENCIES);

        if (!(freqModel.getDataType() instanceof HiddenNucleotides)) {
            throw new IllegalArgumentException("Datatype must be hidden nucleotides!!");
        }

        HiddenNucleotides dataType = (HiddenNucleotides) freqModel.getDataType();

        int hiddenStateCount = dataType.getHiddenClassCount();

        int switchingRatesCount = hiddenStateCount * (hiddenStateCount - 1) / 2;

        if (switchingRates.getDimension() != switchingRatesCount) {
            throw new IllegalArgumentException("switching rates parameter must have " +
                    switchingRatesCount + " dimensions, for " + hiddenStateCount +
                    " hidden categories");
        }

        CovarionHKY model = new CovarionHKY(dataType, kappaParam, hiddenClassRates, switchingRates, freqModel);
        System.out.println(model);
        return model;
    }

}
