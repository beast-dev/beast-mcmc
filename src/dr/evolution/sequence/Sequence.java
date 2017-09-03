/*
 * Sequence.java
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

package dr.evolution.sequence;

import dr.evolution.datatype.DataType;
import dr.evolution.util.Taxon;
import dr.util.Attributable;
import dr.util.Identifiable;

import java.util.Iterator;

/**
 * Class for storing a molecular sequence.
 *
 * @author Alexei Drummond
 * @author Andrew Rambaut
 * @version $Id: Sequence.java,v 1.35 2005/05/25 09:35:28 rambaut Exp $
 */
public class Sequence implements Identifiable, Attributable {

    /**
     * Empty constructor.
     */
    public Sequence() {
        sequenceString = new StringBuffer();
    }

    /**
     * Constructor with initial sequence string.
     *
     * @param sequence a string representing the sequence
     */
    public Sequence(String sequence) {
        sequenceString = new StringBuffer();
        setSequenceString(sequence);
    }

    /**
     * Clone constructor
     *
     * @param sequence the sequence to clone
     */
    public Sequence(Sequence sequence) {
        // should clone taxon as well!
        this(sequence.getTaxon(), sequence.getSequenceString());
    }

    /**
     * Constructor with taxon and sequence string.
     *
     * @param taxon    the sequence's taxon
     * @param sequence the sequence's symbol string
     */
    public Sequence(Taxon taxon, String sequence) {
        sequenceString = new StringBuffer();
        setTaxon(taxon);
        setSequenceString(sequence);
    }

    /**
     * @return the DataType of the sequences.
     */
    public DataType getDataType() {
        return dataType;
    }

    /**
     * @return the length of the sequences.
     */
    public int getLength() {
        return sequenceString.length();
    }

    /**
     * @return a String containing the sequences.
     */
    public String getSequenceString() {
        return sequenceString.toString();
    }

    /**
     * @return a char containing the state at index.
     */
    public char getChar(int index) {
        return sequenceString.charAt(index);
    }

    /**
     * @return the state at site index.
     */
    public int getState(int index) {
        return dataType.getState(sequenceString.charAt(index));
    }

    /**
     */
    public void setState(int index, int state) {

        sequenceString.setCharAt(index, dataType.getChar(state));
    }

    /**
     * Characters are copied from the sequences into the destination character array dst.
     */
    public void getChars(int srcBegin, int srcEnd, char[] dst, int dstBegin) {
        sequenceString.getChars(srcBegin, srcEnd, dst, dstBegin);
    }

    /**
     * search invalid character in the sequence by given data type, and return its index
     */
    public int getInvalidChar() {
        final char[] validChars = dataType.getValidChars();
        if (validChars != null) {
            String validString = new String(validChars);

            for (int i = 0; i < sequenceString.length(); i++) {
                char c = sequenceString.charAt(i);

                if (validString.indexOf(c) < 0) return i;
            }
        }
        return -1;
    }

    /**
     * Set the DataType of the sequences.
     */
    public void setDataType(DataType dataType) {
        this.dataType = dataType;
    }

    /**
     * Set the DataType of the sequences.
     */
    public DataType guessDataType() {
        return DataType.guessDataType(sequenceString.toString());
    }

    /**
     * Set the sequences using a string.
     */
    public void setSequenceString(String sequence) {
        sequenceString.setLength(0);
        sequenceString.append(sequence.toUpperCase());
    }

    /**
     * Append a string to the sequences.
     */
    public void appendSequenceString(String sequence) {
        sequenceString.append(sequence);
    }

    /**
     * Insert a string into the sequences.
     */
    public void insertSequenceString(int offset, String sequence) {
        sequenceString.insert(offset, sequence);
    }

    /**
     * Sets a taxon for this sequences.
     *
     * @param taxon the taxon.
     */
    public void setTaxon(Taxon taxon) {
        this.taxon = taxon;
    }

    /**
     * @return the taxon for this sequences.
     */
    public Taxon getTaxon() {
        return taxon;
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

    // **************************************************************
    // INSTANCE VARIABLES
    // **************************************************************

    protected Taxon taxon = null;
    protected StringBuffer sequenceString = null;
    protected DataType dataType = null;
}


