package dr.evolution.alignment;

import dr.evolution.util.TaxonList;
import dr.util.Attribute;
import dr.xml.*;

import java.util.ArrayList;
import java.util.List;

/**
 * Created by IntelliJ IDEA.
 * User: msuchard
 * Date: Jan 20, 2007
 * Time: 8:30:54 AM
 * To change this template use File | Settings | File Templates.
 */
public class PartitionWithAttributes extends SitePatterns {

	public static final String NAME = "partitionWithAttributes";

	public PartitionWithAttributes(Alignment alignment) {
		super(alignment);
	}

	public PartitionWithAttributes(Alignment alignment, TaxonList taxa) {
		super(alignment, taxa);
	}

	public PartitionWithAttributes(Alignment alignment, int from, int to, int every) {
		super(alignment, from, to, every);
	}

	public PartitionWithAttributes(Alignment alignment, TaxonList taxa, int from, int to, int every) {
		super(alignment, taxa, from, to, every);
	}

	public PartitionWithAttributes(SiteList siteList) {
		super(siteList);
	}

	public PartitionWithAttributes(SiteList siteList, int from, int to, int every) {
		super(siteList, from, to, every);
	}

	public void addAttribute(Attribute a) {
		if (attributes == null)
			attributes = new ArrayList<Attribute>();
		attributes.add(a);
	}



	//************************************************************************
    // AbstractXMLObjectParser implementation
    //************************************************************************

    public static XMLObjectParser PARSER = new AbstractXMLObjectParser() {


	public String getParserName() { return NAME; }

	public Object parseXMLObject(XMLObject xo) throws XMLParseException {

		SitePatterns sites = (SitePatterns)xo.getChild(SitePatterns.class);

//        SitePatterns patterns = new SitePatterns(alignment, taxa, from, to, every
		PartitionWithAttributes partition = new PartitionWithAttributes(
				sites.getSiteList(), sites.getFrom(), sites.getTo(), sites.getEvery()) ;
//
//		) ;
	      for (int i=0; i<xo.getChildCount(); i++) {
		      Object cxo = xo.getChild(i);
		      if (cxo instanceof Attribute)
		        partition.addAttribute((Attribute)cxo);

	      }

		return partition;
	}
    public XMLSyntaxRule[] getSyntaxRules() { return rules; }

    private XMLSyntaxRule[] rules = new XMLSyntaxRule[] {
		new ElementRule(SitePatterns.class, "Must provide a set of site patterns for each partition"),

    };

    public String getParserDescription() {
        return "A weighted list of the unique site patterns (unique columns) in an alignment.";
    }

    public Class getReturnType() { return PatternList.class; }
    };


	private List<Attribute> attributes;

	}
