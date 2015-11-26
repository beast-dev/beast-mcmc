/*
 * ConvertAlignment.java
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

package dr.evolution.alignment;

import dr.evolution.datatype.*;
import dr.evolution.sequence.Sequence;
import dr.evolution.util.Taxon;
import java.util.*;

/**
 * An alignment class that takes another alignment and converts it on the fly
 * to a different dataType.
 *
 * @author Andrew Rambaut
 * @author Alexei Drummond
 *
 * @version $Id: ConvertAlignment.java,v 1.29 2005/05/24 20:25:55 rambaut Exp $
 */
public class ConvertAlignment extends WrappedAlignment implements dr.util.XHTMLable
{

    /**
     * Constructor.
     */
    public ConvertAlignment(DataType dataType) {
        this(dataType, null, null);
    }

    /**
     * Constructor.
     */
    public ConvertAlignment(DataType dataType, CodonTable codonTable) {
        this(dataType, codonTable, null);
    }

    /**
     * Constructor.
     */
    public ConvertAlignment(DataType dataType, Alignment alignment) {
        this(dataType, null, alignment);
    }

    /**
     * Constructor.
     */
    public ConvertAlignment(DataType dataType, CodonTable codonTable, Alignment alignment) {
        super(alignment);
        setDataType(dataType);
        setCodonTable(codonTable);
        setAlignment(alignment);
    }

    /**
     * Sets the CodonTable of this alignment.
     */
    public void setCodonTable(CodonTable codonTable) {
        this.codonTable = codonTable;
    }

    /**
     * Sets the contained.
     */
    public void setAlignment(Alignment alignment) {
        if (dataType == null) {
            dataType = alignment.getDataType();
        }
            
        this.alignment = alignment;

        int newType = dataType.getType();
        int originalType = alignment.getDataType().getType();

      //TODO: this logic does not work for pibuss
        if (originalType == DataType.NUCLEOTIDES) {
        	
            if (newType != DataType.CODONS && newType != DataType.AMINO_ACIDS) {
                throw new RuntimeException("Incompatible alignment DataType for ConversionAlignment");
            }
            
        } else if (originalType == DataType.CODONS) {
        	
            if (!(newType == DataType.AMINO_ACIDS || newType == DataType.NUCLEOTIDES)) {

                System.err.println("originalType = " + originalType);
                System.err.println("newType = " + newType);
                throw new RuntimeException("Incompatible alignment DataType for ConversionAlignment");
            }
            
        } else {
        	
            throw new RuntimeException("Incompatible alignment DataType for ConversionAlignment");
            
        }//END: original type check
    }

    /**
     * Sets the dataType of this alignment. This can be different from
     * the dataTypes of the contained alignment - they will be translated
     * as required.
     */
    public void setDataType(DataType dataType) {
        this.dataType = dataType;
    }

    /**
     * @return the DataType of this siteList
     */
    public DataType getDataType() {
        return dataType;
    }

    /**
     * @return number of sites
     */
    public int getSiteCount() {
        if (alignment == null) throw new RuntimeException("ConvertAlignment has no alignment");

        int originalType = alignment.getDataType().getType();
        int count = alignment.getSiteCount();

        if (originalType == DataType.NUCLEOTIDES) {
            count /= 3;
        }

        return count;
    }

    /**
     * @return the sequence state at (taxon, site)
     */
    public int getState(int taxonIndex, int siteIndex) {
        if (alignment == null) throw new RuntimeException("ConvertAlignment has no alignment");

        int newType = dataType.getType();
        int originalType = alignment.getDataType().getType();

        int state = 0;

        if (originalType == DataType.NUCLEOTIDES) {
            int siteIndex3 = siteIndex * 3;
            int state1 = alignment.getState(taxonIndex, siteIndex3);
            int state2 = alignment.getState(taxonIndex, siteIndex3 + 1);
            int state3 = alignment.getState(taxonIndex, siteIndex3 + 2);

            if (newType == DataType.CODONS) {
                state = ((Codons)dataType).getState(state1, state2, state3);
            } else { // newType == DataType.AMINO_ACIDS
                state = codonTable.getAminoAcidState(((Codons)dataType).getCanonicalState(((Codons)dataType).getState(state1, state2, state3)));
            }

        } else if (originalType == DataType.CODONS) {
            if (newType == DataType.AMINO_ACIDS) {
                state = codonTable.getAminoAcidState(alignment.getState(taxonIndex, siteIndex));
            } else { // newType == DataType.CODONS
                String string = alignment.getAlignedSequenceString(taxonIndex);
                state = Nucleotides.INSTANCE.getState(string.charAt(siteIndex));
            }
        }

        return state;
    }

    public String toXHTML() {
        String xhtml = "<p><em>Converted Alignment</em> data type = ";
        xhtml += getDataType().getDescription();
        xhtml += ", no. taxa = ";
        xhtml += getTaxonCount();
        xhtml += ", no. sites = ";
        xhtml += getSiteCount();
        xhtml += "</p>";

        xhtml += "<pre>";

        int length, maxLength = 0;
        for (int i =0; i < getTaxonCount(); i++) {
            length = getTaxonId(i).length();
            if (length > maxLength)
                maxLength = length;
        }

        int count, state;
        int type = getDataType().getType();

        for (int i = 0; i < getTaxonCount(); i++) {
            length = getTaxonId(i).length();
            xhtml += getTaxonId(i);
            for (int j = length; j <= maxLength; j++)
                xhtml += " ";

            count = getSiteCount();
            for (int j = 0; j < count; j++) {
                state = getState(i, j);
                if (type == DataType.CODONS)
                    xhtml += Codons.UNIVERSAL.getTriplet(state) + " ";
                else
                    xhtml += AminoAcids.INSTANCE.getTriplet(state) + " ";
            }
            xhtml += "\n";
        }
        xhtml += "</pre>";
        return xhtml;
    }

    private DataType dataType = null;
    private CodonTable codonTable = null;
}
