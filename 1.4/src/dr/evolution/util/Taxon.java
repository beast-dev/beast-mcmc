/*
 * Taxon.java
 *
 * Copyright (C) 2002-2006 Alexei Drummond and Andrew Rambaut
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

package dr.evolution.util;

import dr.util.Attributable;
import dr.util.Identifiable;

import java.util.ArrayList;
import java.util.Iterator;

/**
 * Class for data about a taxon.
 *
 * @version $Id: Taxon.java,v 1.24 2006/09/05 13:29:34 rambaut Exp $
 *
 * @author Andrew Rambaut
 * @author Alexei Drummond
 */
public class Taxon implements Attributable, Identifiable, Comparable {

	public Taxon(String id) {
		setId(id);
	}

	/**
	 * Sets a date for this taxon.
	 */
	public void setDate(Date date) {
		setAttribute("date", date);
	}

	/**
	 * @return a date for this taxon.
	 */
	public Date getDate() {
		Object date = getAttribute("date");
		if (date != null && date instanceof Date) {
			return (Date)date;
		}
		return null;
	}

    // **************************************************************
    // Attributable IMPLEMENTATION
    // **************************************************************

	private Attributable.AttributeHelper attributes = null;

	/**
	 * Sets an named attribute for this object.
	 * @param name the name of the attribute.
	 * @param value the new value of the attribute.
	 */
	public void setAttribute(String name, Object value) {
		if (attributes == null)
			attributes = new Attributable.AttributeHelper();
		attributes.setAttribute(name, value);
	}

	/**
	 * @return an object representing the named attributed for this object.
	 * @param name the name of the attribute of interest.
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
	public Iterator getAttributeNames() {
		if (attributes == null)
			return new ArrayList().iterator();
		else
			return attributes.getAttributeNames();
	}

    // **************************************************************
    // Identifiable IMPLEMENTATION
    // **************************************************************

	protected String id = null;

	/**
	 * @return the id.
	 */
	public String getId() {
		return id;
	}

	/**
	 * Sets the id.
	 */
	public void setId(String id) {
		this.id = id;
	}

	public String toString() { return getId(); }

	// **************************************************************
	// Comparable IMPLEMENTATION
	// **************************************************************

	public int compareTo(Object o) {
		return getId().compareTo(((Taxon)o).getId());
	}

}

