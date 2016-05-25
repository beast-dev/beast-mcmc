/*
 * FastaImporter.java
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
import java.util.StringTokenizer;

/**
 * Class for importing PHYLIP sequential file format
 *
 * @version $Id: PhylipSequentialImporter.java,v 1.5 2005/05/24 20:25:56 rambaut Exp $
 *
 * @author Andrew Rambaut
 * @author Alexei Drummond
 */
public class FastaImporter extends Importer implements SequenceImporter {

    public static final char FASTA_FIRST_CHAR = '>';
    
    /**
     * Constructor
     */
    public FastaImporter(Reader reader, DataType dataType) {
        this(reader, null, dataType);
    }

    public FastaImporter(Reader reader, Writer commentWriter, DataType dataType) {
        super(reader, commentWriter);
        setCommentDelimiters('\0', '\0', '\0');

        this.dataType = dataType;
    }

    /**
     * importAlignment.
     */
    public Alignment importAlignment() throws IOException, ImportException
    {
        SimpleAlignment alignment = null;

        try {
            // find fasta line start
            while (read() != FASTA_FIRST_CHAR) {
            }

            do {
                final String name = readLine().trim();
                StringBuffer seq = new StringBuffer();

                readSequence(seq, dataType, "" + FASTA_FIRST_CHAR, Integer.MAX_VALUE, "-", "?", "", "");

                if (alignment == null) {
                    alignment = new SimpleAlignment();
                }

                alignment.addSequence(new Sequence(new Taxon(name.toString()), seq.toString()));

            } while (getLastDelimiter() == FASTA_FIRST_CHAR);
        } catch (EOFException e) {
            // catch end of file the ugly way.
        }

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