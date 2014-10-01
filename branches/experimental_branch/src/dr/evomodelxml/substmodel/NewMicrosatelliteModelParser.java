package dr.evomodelxml.substmodel;

import dr.evolution.datatype.Microsatellite;
import dr.evomodel.substmodel.FrequencyModel;
import dr.evomodel.substmodel.NewMicrosatelliteModel;
import dr.xml.*;


/**
 * @author Chieh-Hsi Wu
 *
 */
public class NewMicrosatelliteModelParser extends AbstractXMLObjectParser{
    public static final String NEW_MSAT_MODEL = "newMsatModel";

    public String getParserName() {
        return NEW_MSAT_MODEL;
    }


    //AbstractXMLObjectParser implementation
    public Object parseXMLObject(XMLObject xo) throws XMLParseException {

        //get msat data type
        System.out.println("Using watkins' model");
        Microsatellite msat = (Microsatellite)xo.getChild(Microsatellite.class);
        //get FrequencyModel
        FrequencyModel freqModel = null;
        if(xo.hasChildNamed(FrequencyModelParser.FREQUENCIES)){
            freqModel = (FrequencyModel)xo.getElementFirstChild(FrequencyModelParser.FREQUENCIES);
        }


        return new NewMicrosatelliteModel(msat, freqModel);
    }


    public String getParserDescription() {
        return "This element represents an instance of the stepwise mutation model of microsatellite evolution.";
    }

    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }

    private XMLSyntaxRule[] rules = new XMLSyntaxRule[]{
            new ElementRule(Microsatellite.class),
            new ElementRule(FrequencyModel.class,true)
    };

    public Class getReturnType() {
        return NewMicrosatelliteModel.class;
    }
}
