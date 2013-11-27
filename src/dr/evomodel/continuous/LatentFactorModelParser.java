package dr.evomodel.continuous;

import dr.inference.model.Parameter;
import dr.xml.*;

/**
 * Created with IntelliJ IDEA.
 * User: max
 * Date: 10/28/13
 * Time: 1:02 PM
 * To change this template use File | Settings | File Templates.
 */
public abstract class LatentFactorModelParser extends AbstractXMLObjectParser {
    public final static String LATENT_FACTOR_MODEL="latentFactorModel";
    public final static String NUMBER_OF_FACTORS="factorNumber";
    public final static String FACTORS="factors";
    public final static String DATA="data";
    public final static String LATENT="latent";

    public String getParserName(){
        return LATENT_FACTOR_MODEL;
    }

    public Object parseXMLObject(XMLObject xo) throws XMLParseException {
        Parameter data = null;
        Parameter factor  = null;
        Parameter latent  = null;


        int factors=xo.getAttribute(NUMBER_OF_FACTORS, 4);
        XMLObject dataXML= xo.getChild(DATA);
        if (dataXML.getChild(0) instanceof Parameter) {
            data = (Parameter) dataXML.getChild(Parameter.class);
        }
//        else {
//            data = new Parameter.Default(dataXML.getDoubleChild(0));
//        }
//exit if there's no parameter? What happens?

        XMLObject factorsXML= xo.getChild(FACTORS);
        if (factorsXML.getChild(0) instanceof Parameter) {
            factor = (Parameter) factorsXML.getChild(Parameter.class);
        }
//        else {
//            factor = new Parameter.Default(factorsXML.getDoubleChild(0));
//        }


        XMLObject latentXML= xo.getChild(LATENT);
        if (latentXML.getChild(0) instanceof Parameter) {
            latent = (Parameter) latentXML.getChild(Parameter.class);
        }
//        else {
//            latent = new Parameter.Default(latentXML.getDoubleChild(0));
//        }

        return new LatentFactorModel(data, factor, latent);
    }


}
