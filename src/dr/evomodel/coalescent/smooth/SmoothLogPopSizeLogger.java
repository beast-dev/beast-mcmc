package dr.evomodel.coalescent.smooth;

import dr.inference.loggers.LogColumn;
import dr.inference.loggers.Loggable;
import dr.inference.loggers.NumberColumn;
import dr.xml.*;


/**
 *
 * @author Xiang Ji
 * @author Marc A. Suchard
 */
public class SmoothLogPopSizeLogger implements Loggable {

    private final NewSmoothSkygridLikelihood likelihood;

    public SmoothLogPopSizeLogger(NewSmoothSkygridLikelihood likelihood) {
        this.likelihood = likelihood;
    }

    @Override
    public LogColumn[] getColumns() {
        int size = likelihood.getGridPointParameter().getDimension();

        LogColumn[] array = new LogColumn[size];

        for (int i = 0; i < size; ++i) {

            final int index = i;
            array[i] = new NumberColumn("log.smooth.pop.size." + (index + 1)) {
                @Override
                public double getDoubleValue() {
                    return -Math.log(likelihood.getSmoothPopSizeInverse(likelihood.getGridPointParameter().getParameterValue(index)));
                }
            };
        }

        return array;
    }

    private static final String SMOOTH_LOG_POP_SIZE_PARAMETER = "smoothLogPopSizeParameter";
    public static XMLObjectParser PARSER = new AbstractXMLObjectParser() {

        @Override
        public Object parseXMLObject(XMLObject xo) throws XMLParseException {
            NewSmoothSkygridLikelihood likelihood = (NewSmoothSkygridLikelihood) xo.getChild(NewSmoothSkygridLikelihood.class);
            return new SmoothLogPopSizeLogger(likelihood);
        }

        @Override
        public String getParserName() {
            return SMOOTH_LOG_POP_SIZE_PARAMETER;
        }

        @Override
        public XMLSyntaxRule[] getSyntaxRules() {
            return rules;
        }

        XMLSyntaxRule[] rules = new XMLSyntaxRule[]{
                new ElementRule(NewSmoothSkygridLikelihood.class)
        };

        @Override
        public String getParserDescription() {
            return "Logs the log of smoothed pop sizes of smooth Skygrid model";
        }

        @Override
        public Class getReturnType() {
            return SmoothLogPopSizeLogger.class;
        }
    };



}
