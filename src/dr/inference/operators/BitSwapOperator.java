package dr.inference.operators;

import dr.xml.*;
import dr.inference.model.Parameter;
import dr.inference.model.Statistic;
import dr.inference.distribution.MixedDistributionLikelihood;
import dr.math.MathUtils;

/**
 * Given a values vector (data) and an indicators vector (boolean vector indicating wheather the corrosponding value
 * is used or ignored), this operator explores all possible positions for the used data points while preserving their
 * order.
 * The distribition is uniform on all possible data positions.
 *
 * For example, if data values A and B are used in a vector of dimension 4, each of the following states is visited 1/6
 * of the time.
 *
 * ABcd 1100
 * AcBd 1010
 * AcdB 1001
 * cABd 0110
 * cAdB 0101
 * cdAB 0011
 *
 * The operator works by picking a 1 bit in the indicators and swapping it with a neighbour 0, with the appropriate
 * adjustment to the hastings ratio since a pair of 1,1 and 0,0 are never swapped, and the ends can be swapped in one 
 * direction only.
 *  
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

        final int iDim = indicators.getDimension();
        final int dDim = data.getDimension();
        if (iDim == dDim -1) {
            impliedOne = true;
        } else if (iDim == dDim) {
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
        if( dim < 2 ) {
            throw new OperatorFailedException("no swaps possible");
        }
        int nLoc = 0;
        int[] loc = new int[2*dim];
        double prev = -1;
        int nOnes = 0;
        for (int i = 0; i < dim; i++) {
            final double value = indicators.getStatisticValue(i);
            if( value > 0 ) {
                ++nOnes;
                if( i > 0 && prev == 0 ) {
                    loc[nLoc] = -(i+1);
                    ++nLoc;
                }
                if( i < dim-1 && indicators.getStatisticValue(i+1) == 0 ) {
                   loc[nLoc] = (i+1);
                    ++nLoc;
                }
            }
            prev = value;
        }

        if( nOnes == 0 ) {
            return 0;
        }
        assert nLoc > 0;

        final int rand = MathUtils.nextInt(nLoc);
        int pos = loc[rand];
        int direction = pos < 0 ? -1 : 1;
        pos = (pos < 0 ? -pos : pos) - 1;
        final int maxOut = 2 * nOnes;

        double hastingsRatio = maxOut == nLoc ? 0.0 : Math.log((double)nLoc/maxOut);

//            System.out.println("swap " + pos + "<->" + nto + "  " +
//                              indicators.getParameterValue(pos) +  "<->" + indicators.getParameterValue(nto) +
//                 "  " +  data.getParameterValue(pos) +  "<->" + data.getParameterValue(nto));
        final int nto = pos + direction;
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
