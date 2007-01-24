package dr.inference.operators;

import dr.xml.*;
import dr.inference.model.Parameter;
import dr.inference.model.Statistic;
import dr.inference.distribution.MixedDistributionLikelihood;
import dr.math.MathUtils;

/**
 * @author Joseph Heled
 * @version $Id$
 */
public class BitSwapOperator extends SimpleMCMCOperator {

    public static final String BIT_SWAP_OPERATOR = "bitSwapOperator";
    private Parameter data;
    private Parameter indicators;
    private final boolean impliedOne;

    public BitSwapOperator(Parameter data, Parameter indicators, int weight) {
        this.data = data;
        this.indicators = indicators;
        setWeight(weight);

         if (indicators.getDimension() == data.getDimension()-1) {
            impliedOne = true;
        } else if (indicators.getDimension() == data.getDimension()) {
             impliedOne = false;
         } else {
            throw new IllegalArgumentException();
        }
    }


    public String getPerformanceSuggestion() {
        return "";
    }

    public String getOperatorName() {
        return BIT_SWAP_OPERATOR;   // todo is that right, seems to conflict with bitSwap
    }

    public double doOperation() throws OperatorFailedException {
        final int dim = indicators.getDimension();
        //double sum = 0.0;
        int[] loc = new int[dim];
        int nLoc = 0;
        for (int i = 0; i < dim; i++) {
            final double value = indicators.getStatisticValue(i);
            if( value > 0 ) {
               //sum += value;
                loc[nLoc] = i;
                ++nLoc;
            }
        }

        if( nLoc > 0 ) {
            double hastingsRatio = 0;
            final int rand = MathUtils.nextInt(nLoc);
            final int pos = loc[rand];
            int direction;
            if( pos == 0 ) {
                direction = +1;
                hastingsRatio = Math.log(2.0);
            } else if( pos == dim -1 ) {
                direction = -1;
                hastingsRatio = Math.log(2.0);
            } else {
                direction = MathUtils.nextInt(2) == 0 ? -1 : 1;
            }
            final int nto = pos + direction;
            if( nto == 0 || nto == dim - 1) {
                hastingsRatio = Math.log(1.0/2.0);
            }

//            System.out.println("swap " + pos + "<->" + nto + "  " +
//                              indicators.getParameterValue(pos) +  "<->" + indicators.getParameterValue(nto) +
//                 "  " +  data.getParameterValue(pos) +  "<->" + data.getParameterValue(nto));

            double vto = indicators.getStatisticValue(nto);
            indicators.setParameterValue(nto, indicators.getParameterValue(pos));
            indicators.setParameterValue(pos, vto);

            final int dataOffset = impliedOne ? 1 : 0;
            final int ntodata = nto + dataOffset;
            final int posdata = pos + dataOffset;
            vto = data.getStatisticValue(ntodata);
            data.setParameterValue(ntodata, data.getParameterValue(posdata));
            data.setParameterValue(posdata, vto);

//            System.out.println("after " + pos + "<->" + nto + "  " +
//                              indicators.getParameterValue(pos) +  "<->" + indicators.getParameterValue(nto) +
//                 "  " +  data.getParameterValue(pos) +  "<->" + data.getParameterValue(nto));

            return hastingsRatio;
        }
        return 0;
    }

    private static final String DATA = MixedDistributionLikelihood.DATA;
    private static final String INDICATORS = MixedDistributionLikelihood.INDICATORS;
    public static XMLObjectParser PARSER = new AbstractXMLObjectParser() {

        public String getParserName() { return BIT_SWAP_OPERATOR; }

        public Object parseXMLObject(XMLObject xo) throws XMLParseException {

            final int weight = xo.getIntegerAttribute(WEIGHT);

            Parameter data = (Parameter )((XMLObject)xo.getChild(DATA)).getChild(Parameter.class);
            Parameter indicators = (Parameter)((XMLObject)xo.getChild(INDICATORS)).getChild(Parameter.class);
            return new BitSwapOperator(data, indicators, weight);
        }

        //************************************************************************
        // AbstractXMLObjectParser implementation
        //************************************************************************

        public String getParserDescription() {
            return "This element returns a bit-swap operator on a given parameter and data.";
        }

        public Class getReturnType() { return MCMCOperator.class; }

        public XMLSyntaxRule[] getSyntaxRules() { return rules; }

        private XMLSyntaxRule[] rules = new XMLSyntaxRule[] {
                AttributeRule.newIntegerRule(WEIGHT),
                new ElementRule(DATA, new XMLSyntaxRule[] {new ElementRule(Statistic.class)}),
                new ElementRule(INDICATORS, new XMLSyntaxRule[] {new ElementRule(Statistic.class)}),
        };

    };
}
