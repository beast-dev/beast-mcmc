/*
 * MetagenomeData.java
 *
 * Copyright (c) 2002-2015 Alexei Drummond, Andrew Rambaut and Marc Suchard
 *
 * This file is part of BEAST.
 * See the NOTICE file distributed with this work for additional
 * information regarding copyright ownership and licensing.
 *
 * BEAST is free software; you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as
 * published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.
 *
 *  BEAST is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with BEAST; if not, write to the
 * Free Software Foundation, Inc., 51 Franklin St, Fifth Floor,
 * Boston, MA  02110-1301  USA
 */

package dr.evolution;

import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;

import dr.evolution.alignment.Alignment;
import dr.evolution.tree.Tree;
import dr.evolution.util.Taxa;
import dr.evolution.util.Taxon;
import dr.evolution.util.TaxonList;
import dr.util.Attributable;
import dr.util.Identifiable;

/**
 * @author Aaron Darling (koadman)
 */
public class MetagenomeData implements Attributable, Identifiable {

	protected String id = null;
	protected Alignment alignment;
	protected Taxa referenceTaxa;
	protected Taxa reads;
	protected LinkageConstraints constraints;
	protected boolean fixedReferenceTree;

	public MetagenomeData(TaxonList taxa, Alignment alignment)
	{
		this(taxa, alignment, null, false);
	}
	public MetagenomeData(TaxonList taxa, Alignment alignment, LinkageConstraints constraints, boolean fixedReferenceTree)
	{
		this.alignment = alignment;
		this.constraints = constraints;
		this.fixedReferenceTree = fixedReferenceTree;
		referenceTaxa = new Taxa(taxa.asList());
		if(fixedReferenceTree == true && !(taxa instanceof Tree)){
			throw new RuntimeException("Error, a tree must be provided in order to fix the topology of reference taxa");
		}
		reads = new Taxa();
		List<Taxon> alltaxa = alignment.asList();
		for(int i=0; i<alltaxa.size(); i++)
		{
			if(!referenceTaxa.contains(alltaxa.get(i)))
				reads.addTaxon(alltaxa.get(i));
		}
		if(constraints==null){
			this.constraints = new LinkageConstraints(new ArrayList<LinkedGroup>());
		}
	}
	
	public TaxonList getReferenceTaxa(){
		return referenceTaxa;
	}
	public TaxonList getReadsTaxa(){
		return reads;
	}
	
	public Alignment getAlignment() {
		return alignment;
	}

	public void setAlignment(Alignment alignment) {
		this.alignment = alignment;
	}

	public ArrayList<LinkedGroup> getConstraints(){
		return constraints.getLinkedGroups();
	}

	public boolean getFixedReferenceTree(){
		return fixedReferenceTree;
	}
	public void setFixedReferenceTree(boolean fixedReferenceTree){
		this.fixedReferenceTree = fixedReferenceTree;
	}

	public String getId() {
		return id;
	}

	public void setId(String id) {
		this.id = id;
	}
    // **************************************************************
    // Attributable IMPLEMENTATION
    // **************************************************************

    private Attributable.AttributeHelper attributes = null;

    /**
     * Sets an named attribute for this object.
     *
     * @param name  the name of the attribute.
     * @param value the new value of the attribute.
     */
    public void setAttribute(String name, Object value) {
        if (attributes == null)
            attributes = new Attributable.AttributeHelper();
        attributes.setAttribute(name, value);
    }

    /**
     * @param name the name of the attribute of interest.
     * @return an object representing the named attributed for this object.
     */
    public Object getAttribute(String name) {
        if (attributes == null)
            return null;
        else
            return attributes.getAttribute(name);
    }

    /**
     * @return an iterator of the attributes that this object has.
     */
    public Iterator<String> getAttributeNames() {
        if (attributes == null)
            return null;
        else
            return attributes.getAttributeNames();
    }

}
