/*
 * PhylipSequentialImporter.java
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

package dr.evolution.io;

import dr.evolution.alignment.Alignment;
import dr.evolution.alignment.SimpleAlignment;
import dr.evolution.datatype.DataType;
import dr.evolution.sequence.Sequence;
import dr.evolution.sequence.SequenceList;
import dr.evolution.util.Taxon;

import java.io.EOFException;
import java.io.IOException;
import java.io.Reader;
import java.io.Writer;

/**
 * Class for importing PHYLIP sequential file format
 *
 * @version $Id: PhylipSequentialImporter.java,v 1.5 2005/05/24 20:25:56 rambaut Exp $
 *
 * @author Andrew Rambaut
 * @author Alexei Drummond
 */
public class PhylipSequentialImporter extends Importer implements SequenceImporter { 

	/**
	 * Constructor
	 */
	public PhylipSequentialImporter(Reader reader, DataType dataType, int maxNameLength) {
		super(reader);
		setCommentDelimiters('\0', '\0', '\0');
		
		this.dataType = dataType;
		this.maxNameLength = maxNameLength;
	}
	
	public PhylipSequentialImporter(Reader reader, Writer commentWriter, DataType dataType, int maxNameLength) {
		super(reader, commentWriter);
		setCommentDelimiters('\0', '\0', '\0');
		
		this.dataType = dataType;
		this.maxNameLength = maxNameLength;
	}
	
	/**
	 * importAlignment. 
	 */
	public Alignment importAlignment() throws IOException, Importer.ImportException
	{
		SimpleAlignment alignment = null;
		
		try {
		
			int taxonCount = readInteger();
			int siteCount = readInteger();
		
			String firstSeq = null;
		
			for (int i = 0; i < taxonCount; i++) {
				StringBuffer name = new StringBuffer();
				
				char ch = read();
				int n = 0;
				while (!Character.isWhitespace(ch) && (maxNameLength < 1 || n < maxNameLength)) {
					name.append(ch);
					ch = read();
					n++;
				}
				
				StringBuffer seq = new StringBuffer(siteCount);
				readSequence(seq, dataType, "", siteCount, "-", "?", ".", firstSeq);
				
				if (firstSeq == null) { firstSeq = seq.toString(); }
		
		
				if (alignment == null) {
					alignment = new SimpleAlignment();
				}
				
				alignment.addSequence(new Sequence(new Taxon(name.toString()), seq.toString()));
			}
			
			
		} catch (EOFException e) { }	

		return alignment;
	}

	/**
	 * importSequences. 
	 */
	public SequenceList importSequences() throws IOException, ImportException {
		return importAlignment();
	}

	private DataType dataType;
	private int maxNameLength = 10;	
}
