package dr.inference.operators;

import dr.inference.model.ValuesPool;
import dr.inference.model.Variable;
import dr.math.MathUtils;
import dr.xml.*;

/**
 * @author Joseph Heled
 *         Date: 8/09/2009
 */
public class ValuesPoolSwapOperator extends SimpleMCMCOperator {
    public static String VALUESPOOL_OPERATOR = "poolSwapOperator";

    private final ValuesPool pool;

    ValuesPoolSwapOperator(ValuesPool pool) {
        this.pool = pool;
    }

    public String getOperatorName() {
        return VALUESPOOL_OPERATOR + "(" + pool.getModelName() + ")";
    }

    public double doOperation() throws OperatorFailedException {
        final Variable<Double> selector = pool.getSelector();
        final int[] ints = SelectorOperator.intVals(selector);
        int[] c = SelectorOperator.counts_used_m2(ints);

        int n = 0;
        for(int k = 0; k < c.length-1; ++k) {
            if( c[k] == c[k+1] ) {
              ++n;
            }
        }

        if( n == 0 ) {
            throw new OperatorFailedException("No moves");
        }

        int j = MathUtils.nextInt(n);
        int p;
        for(p = 0; p < c.length-1; ++p) {
            if( c[p] == c[p+1] ) {
                if( j == 0 ) {
                    break;
                }
                j -= 1;
            }
        }

        int ip = -1, ip1 = -1;
        // exchange p's and p+1's
        int count = c[p];
        while( count > 0 ) {
            while( ints[++ip] != p );
            while( ints[++ip1] != p+1 );
            selector.setValue(ip, (double)(p+1));
            selector.setValue(ip1, (double)(p));
            --count;
        }
        final Variable<Double> vals = pool.getPool();
        final Double vp = vals.getValue(p);
        final Double vp1 = vals.getValue(p+1);
        vals.setValue(p, vp1);
        vals.setValue(p+1, vp);

        return 0;
    }

    public String getPerformanceSuggestion() {
        return null;
    }


    public static dr.xml.XMLObjectParser PARSER = new dr.xml.AbstractXMLObjectParser() {

        public String getParserName() {
            return VALUESPOOL_OPERATOR;
        }

        public Object parseXMLObject(XMLObject xo) throws XMLParseException {
            final ValuesPool parameter = (ValuesPool) xo.getChild(ValuesPool.class);
            final double weight = xo.getDoubleAttribute(WEIGHT);
            final MCMCOperator op = new ValuesPoolSwapOperator(parameter) ;
            op.setWeight(weight);
            return op;
        }

        //************************************************************************
        // AbstractXMLObjectParser implementation
        //************************************************************************

        public String getParserDescription() {
            return ".";
        }

        public Class getReturnType() {
            return SelectorOperator.class;
        }

        public XMLSyntaxRule[] getSyntaxRules() {
            return new XMLSyntaxRule[] {
                    AttributeRule.newDoubleRule(WEIGHT),
                    new ElementRule(ValuesPool.class),
            };
        }
    };
}
