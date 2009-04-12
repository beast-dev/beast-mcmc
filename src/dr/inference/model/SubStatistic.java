package dr.inference.model;

import dr.xml.AbstractXMLObjectParser;
import dr.xml.AttributeRule;
import dr.xml.ElementRule;
import dr.xml.XMLObject;
import dr.xml.XMLObjectParser;
import dr.xml.XMLParseException;
import dr.xml.XMLSyntaxRule;

public class SubStatistic extends Statistic.Abstract {

	public static final String SUB_STATISTIC = "subStatistic";
	public static final String DIMENSION = "dimension";
	
	private final int[] dimensions;
	private final Statistic statistic;
	
	public SubStatistic(String name, int[] dimensions, Statistic stat){
		super(name);
		this.dimensions = dimensions;
		this.statistic = stat;
	}
	
	public int getDimension() {
		return dimensions.length;
	}

	public double getStatisticValue(int dim) {
		return statistic.getStatisticValue(dimensions[dim]);
	}

    public String getDimensionName(int dim) {
      return statistic.getDimensionName(dimensions[dim]);
    }
    
	public static XMLObjectParser PARSER = new AbstractXMLObjectParser() {

        public String getParserName() {
            return SUB_STATISTIC;
        }

        public Object parseXMLObject(XMLObject xo) throws XMLParseException {

            final String name = xo.getAttribute(NAME, xo.getId());

            final Statistic stat = (Statistic) xo.getChild(Statistic.class);
            
            final int[] values = xo.getIntegerArrayAttribute(DIMENSION);
            
            if( values.length == 0 ){
            	throw new XMLParseException("Must specify at least one dimension");
            }
            
            final int dim = stat.getDimension();
            
            for( int value : values ) {
            	if( value >= dim || value < 0 ) {
            		throw new XMLParseException("Dimension " + value + " is not a valid dimension.");
            	}
            }

            return new SubStatistic(name, values, stat);
        }

        //************************************************************************
        // AbstractXMLObjectParser implementation
        //************************************************************************

        public String getParserDescription() {
            return "Allows you to choose specific dimensions of a given statistic";
        }

        public Class getReturnType() {
            return SubStatistic.class;
        }

        public XMLSyntaxRule[] getSyntaxRules() {
            return new XMLSyntaxRule[]{
            		AttributeRule.newIntegerArrayRule(DIMENSION, false),
            		new ElementRule(Statistic.class),
            };
        }
    };
}
