package dr.evomodel.continuous;

import dr.inference.model.Statistic;
import dr.xml.*;

/**
 *
 * Gabriela Cybis
 *
 *  */
public class RootTraitStatistic extends Statistic.Abstract{


    private int dimension;
    public static final String SAMPLED_ROOT_TRAITS = "sampledRootTraits";
    private IntegratedMultivariateTraitLikelihood likelihood;

    public RootTraitStatistic( IntegratedMultivariateTraitLikelihood likelihood, String id){
        this.likelihood=likelihood;
        setId(id);

        this.dimension = likelihood.dimTrait;



    }




    public int getDimension() {
        return dimension;
    }
    public String getDimensionName(int dim){
        return "root." + likelihood.getTraitName() + dim;
    }

    public double getStatisticValue(int dim) {

    return    likelihood.getRootNodeTrait()[dim];
    }








public static XMLObjectParser PARSER = new AbstractXMLObjectParser() {

    public Object parseXMLObject(XMLObject xo) throws XMLParseException {

        IntegratedMultivariateTraitLikelihood likelihood = (IntegratedMultivariateTraitLikelihood)
                xo.getChild(IntegratedMultivariateTraitLikelihood.class);

        String id = xo.getId();

        return new RootTraitStatistic(likelihood, id);
    }

    public XMLSyntaxRule[] getSyntaxRules() {
        return new XMLSyntaxRule[] {
                new ElementRule(IntegratedMultivariateTraitLikelihood.class),

        };
    }

    public String getParserDescription() {
        return null;
    }

    public Class getReturnType() {
        return RootTraitStatistic.class;
    }

    public String getParserName() {
        return SAMPLED_ROOT_TRAITS;
    }
};

}