package dr.evomodel.continuous;

import dr.geo.cartogram.CartogramMapping;
import dr.inference.model.Parameter;
import dr.xml.*;

import java.awt.geom.Rectangle2D;
import java.util.ArrayList;
import java.util.List;
import java.util.logging.Logger;

/**
 * @author Marc Suchard
 */

public class MixtureCartogramDiffusionModel extends CartogramDiffusionModel {

    public static final String DIFFUSION_PROCESS = "mixtureCartogramDiffusionModel";
    public static final String MIXTURE = "mixture";
    public static final String MAP = "map";

//    public MixtureCartogramDiffusionModel(String name, Parameter precision) {
//        super(name, precision);
//    }

    public MixtureCartogramDiffusionModel(String name, Parameter precision, Parameter mixture) {
        super(name, precision);
        this.mixture = mixture;
        Logger.getLogger("dr.evomodel.continuous").info(
                "\tMixture: " + mixture.getId()
        );
        addVariable(mixture);
    }

    public void addMapping(CartogramMapping mapping) {

        if (mappingList == null)
            mappingList = new ArrayList<CartogramMapping>();

        mappingList.add(mapping);
        Logger.getLogger("dr.evomodel.continuous").info(
            "\tMapping  : " + mapping.toString() + " with density "+mapping.getAverageDensity()
        );
    }

    protected CartogramMapping getMapping() {
        final int index = (int)mixture.getParameterValue(0)-1;
//        System.err.println("index = "+index);
//        if (index > 0)
//        System.exit(-1);

        return mappingList.get(index);
    }

    public static XMLObjectParser PARSER = new AbstractXMLObjectParser() {

           public String getParserName() {
               return DIFFUSION_PROCESS;
           }

           public Object parseXMLObject(XMLObject xo) throws XMLParseException {

               Rectangle2D boundingBox = parseRectangle2D(xo);
         
               Parameter diffusionParam = (Parameter) xo.getChild(Parameter.class);

               Parameter mixtureParam = (Parameter) xo.getChild(MIXTURE).getChild(Parameter.class);

               MixtureCartogramDiffusionModel model = new MixtureCartogramDiffusionModel(xo.getId(),
                       diffusionParam, mixtureParam);

               for(int i=0; i<xo.getChildCount(); i++) {
                   if (xo.getChildName(i).equals(MAP)) {
                       XMLObject cxo = (XMLObject) xo.getChild(i);
                       CartogramMapping mapping = parseCartogramMapping(cxo,boundingBox);
                       model.addMapping(mapping);
                   }
               }

               return model;
           }

           //************************************************************************
           // AbstractXMLObjectParser implementation
           //************************************************************************

           public String getParserDescription() {
               return "Describes a bivariate diffusion process using cartogram distances.";
           }

           public XMLSyntaxRule[] getSyntaxRules() {
               return rules;
           }

           private final XMLSyntaxRule[] rules = {
                   new ElementRule(Parameter.class),
                   new ElementRule(MIXTURE, new XMLSyntaxRule[] {
                           new ElementRule(Parameter.class)
                   }),
//                   new ElementRule(MAP,1,Integer.MAX_VALUE),
//                   AttributeRule.newStringRule(FILENAME,true),
//                   AttributeRule.newIntegerRule(XSIZE),
//                   AttributeRule.newIntegerRule(YSIZE),
                   boundingBoxRules,
           };

           public Class getReturnType() {
               return MultivariateDiffusionModel.class;
           }
       };

      private final Parameter mixture;
      private List<CartogramMapping> mappingList;
}
