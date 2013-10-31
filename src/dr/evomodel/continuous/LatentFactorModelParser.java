package dr.evomodel.continuous;

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
    public final static String ALIGNMENT="Alignment";

    public String getParserName(){
        return LATENT_FACTOR_MODEL;
    }

    public Object parseXMLObject(XMLObject xo) throws XMLParseException {
        int factors=xo.getAttribute(NUMBER_OF_FACTORS, 4);
        XMLObject data= xo.getChild(ALIGNMENT);
        return new LatentFactorModel(factors, data);
    }


}
