/*
 * Importer.java
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

import dr.evolution.datatype.DataType;

import java.io.*;

/**
 * Base class for phylogenetic file format importers
 *
 * @version $Id: Importer.java,v 1.19 2006/04/25 14:41:08 rambaut Exp $
 *
 * @author Andrew Rambaut
 * @author Alexei Drummond
 */
public abstract class Importer {

	public static class ImportException extends Exception {
		/**
		 *
		 */
		private static final long serialVersionUID = 7858834683324203750L;
		public ImportException() { super(); }
		public ImportException(String message) { super(message); }
	}

    public static class DuplicateFieldException extends ImportException {
		/**
		 *
		 */
		private static final long serialVersionUID = 8047146381348414810L;
		public DuplicateFieldException() { super(); }
		public DuplicateFieldException(String message) { super(message); }
	}

    public static class BadFormatException extends ImportException {
		/**
		 *
		 */
		private static final long serialVersionUID = -8206831989674620748L;
		public BadFormatException() { super(); }
		public BadFormatException(String message) { super(message); }
	}

    public static class UnparsableDataException extends ImportException {
		/**
		 *
		 */
		private static final long serialVersionUID = 5905130039882401006L;
		public UnparsableDataException() { super(); }
		public UnparsableDataException(String message) { super(message); }
	}

    public static class MissingFieldException extends ImportException {
		/**
		 *
		 */
		private static final long serialVersionUID = -7576489210458327552L;
		public MissingFieldException() { super(); }
		public MissingFieldException(String message) { super(message); }
	}

    public static class ShortSequenceException extends ImportException {
		/**
		 *
		 */
		private static final long serialVersionUID = 7460033398106047073L;
		public ShortSequenceException() { super(); }
		public ShortSequenceException(String message) { super(message); }
	}

    public static class TooFewTaxaException extends ImportException {
		/**
		 *
		 */
		private static final long serialVersionUID = -6349041350075169247L;
		public TooFewTaxaException() { super(); }
		public TooFewTaxaException(String message) { super(message); }
	}

    public static class UnknownTaxonException extends ImportException {
		/**
		 *
		 */
		private static final long serialVersionUID = 6611115782536515250L;
		public UnknownTaxonException() { super(); }
		public UnknownTaxonException(String message) { super(message); }
	}

    /**
	 * Constructor
	 */
	public Importer(Reader reader) {
		this.reader = new LineNumberReader(reader);
		this.commentWriter = null;
	}

	public Importer(Reader reader, Writer commentWriter) {
		this.reader = new LineNumberReader(reader);
		this.commentWriter = commentWriter != null ? new BufferedWriter(commentWriter) : null;
	}

	public void setCommentDelimiters(char line) {
		hasComments = true;
		this.lineComment = line;
	}

	public void setCommentDelimiters(char start, char stop) {
		hasComments = true;
		this.startComment = start;
		this.stopComment = stop;
	}

    public void setCommentDelimiters(char start, char stop, char line) {
        hasComments = true;
        this.startComment = start;
        this.stopComment = stop;
        this.lineComment = line;
    }

	public void setCommentDelimiters(char start, char stop, char line, char write, char meta) {
		hasComments = true;
		this.startComment = start;
		this.stopComment = stop;
		this.lineComment = line;
        this.writeComment = write;
        this.metaComment = meta;
	}

	public void setCommentWriter(Writer commentWriter) {
		this.commentWriter = new BufferedWriter(commentWriter);
	}

	public int getLineNumber() {
		return reader.getLineNumber();
	}

	public int getLastDelimiter() {
		return lastDelimiter;
	}

	public char nextCharacter() throws IOException {
		if (lastChar == '\0') {
			lastChar = readCharacter();
		}
		return (char)lastChar;
	}

	public char readCharacter() throws IOException {

		skipSpace();

		char ch = read();

		while (hasComments && (ch == startComment || ch == lineComment)) {
			skipComments(ch);
			skipSpace();
			ch = read();
		}

		return ch;
	}

	public void unreadCharacter(char ch) {
		lastChar = ch;
	}

	public char next() throws IOException {
		if (lastChar == '\0') {
			lastChar = read();
		}
		return (char)lastChar;
	}

	public char read() throws IOException {
		int ch;

		if (lastChar == '\0') {
			ch = reader.read();
			if (ch <= 0) {
				throw new EOFException();
			}
		} else {
			ch = lastChar;
			lastChar = '\0';
		}

		return (char)ch;
	}

	/**
	 * Reads a line, skipping over any comments.
	 */
	public String readLine() throws IOException {

		StringBuffer line = new StringBuffer();

		char ch = read();

		try {

			while (ch != '\n' && ch != '\r') {

				if (hasComments) {
					if (ch == lineComment) {
						skipComments(ch);
						break;
					}
					if (ch == startComment) {
						skipComments(ch);
						ch = read();
					}
				}

				line.append(ch);
				ch = read();
			}

			// accommodate DOS line endings..
			if (ch == '\r') {
				if (next() == '\n') read();
			}

			lastDelimiter = ch;

		} catch (EOFException e) {
			// We catch an EOF and return the line we have so far
		}

		return line.toString();
	}

	/**
	 * Reads sequences, skipping over any comments and filtering using dataType.
	 * @param sequence a StringBuffer into which the sequences is put
	 * @param dataType the dataType of the sequences
	 * @param delimiters list of characters that will stop the reading
	 * @param maxSites maximum number of sites to read
	 */
	public void readSequence(StringBuffer sequence, DataType dataType, 
								String delimiters, int maxSites,
								String gapCharacters, String missingCharacters,
								String matchCharacters, String matchSequence) throws IOException, ImportException {
			
		char ch = read();

		try {
			int n = 0;

			while (n < maxSites && delimiters.indexOf(ch) == -1) {

				if (hasComments && (ch == startComment || ch == lineComment)) {
					skipComments(ch);
					ch = read();
				}

				if (!Character.isWhitespace(ch)) {

					char ch1 = ch;

					if (gapCharacters.indexOf(ch) != -1) {
						ch1 = DataType.GAP_CHARACTER;
					} else if (missingCharacters.indexOf(ch) != -1) {
						ch1 = DataType.UNKNOWN_CHARACTER;
					} else if (matchCharacters.indexOf(ch) != -1) {
						if (matchSequence == null) {
							throw new ImportException("Match character in first sequences");
						}
						if (n >= matchSequence.length()) {
							throw new ImportException("Match sequences too short");
						}

						ch1 = matchSequence.charAt(n);
					}

					sequence.append(ch1);

					n++;
				}

				ch = read();
			}

			lastDelimiter = ch;

			if (Character.isWhitespace((char)lastDelimiter)) {
				ch = nextCharacter();
				if (delimiters.indexOf(ch) != -1) {
					lastDelimiter = readCharacter();
				}
			}

		} catch (EOFException e) {
			// We catch an EOF and return the sequences we have so far
		}
	}

	/**
	 * Reads a line of sequences, skipping over any comments and filtering using dataType.
	 * @param sequence a StringBuffer into which the sequences is put
	 * @param dataType the dataType of the sequences
	 * @param delimiters list of characters that will stop the reading
	 */
	public void readSequenceLine(StringBuffer sequence, DataType dataType,
								String delimiters,
								String gapCharacters, String missingCharacters,
								String matchCharacters, String matchSequence) throws IOException, ImportException {

		char ch = read();

		try {
			int n = 0;

			while (ch != '\r' && ch != '\n' && delimiters.indexOf(ch) == -1) {

				if (hasComments) {
					if (ch == lineComment) {
						skipComments(ch);
						break;
					}
					if (ch == startComment) {
						skipComments(ch);
						ch = read();
					}
				}

				if (ch != ' ' && ch != '\t') {
					char ch1 = ch;

					if (gapCharacters.indexOf(ch) != -1) {
						ch1 = DataType.GAP_CHARACTER;
					} else if (missingCharacters.indexOf(ch) != -1) {
						ch1 = DataType.UNKNOWN_CHARACTER;
					} else if (matchCharacters.indexOf(ch) != -1) {
						if (matchSequence == null) {
							throw new ImportException("Match character in first sequences");
						}
						if (n >= matchSequence.length()) {
							throw new ImportException("Match sequences too short");
						}

						ch1 = matchSequence.charAt(n);
					}

					sequence.append(ch1);

					n++;
				}

				ch = read();
			}

			if (ch == '\r') {
				if (next() == '\n') read();
			}

			lastDelimiter = ch;

			if (Character.isWhitespace((char)lastDelimiter)) {
				ch = nextCharacter();
				if (delimiters.indexOf(ch) != -1) {
					lastDelimiter = readCharacter();
				}
			}

		} catch (EOFException e) {
			// We catch an EOF and return the sequences we have so far
		}
	}

	/**
	 * Attempts to read and parse an integer delimited by whitespace.
	 */
	public int readInteger() throws IOException, ImportException {
		String token = readToken();
		try {
			return Integer.parseInt(token);
		} catch (NumberFormatException nfe) {
			throw new ImportException("Number format error: " + nfe.getMessage());
		}
	}

	/**
	 * Attempts to read and parse an integer delimited by whitespace or by
	 * any character in delimiters.
	 */
	public int readInteger(String delimiters) throws IOException, ImportException {
		String token = readToken(delimiters);
		try {
			return Integer.parseInt(token);
		} catch (NumberFormatException nfe) {
			throw new ImportException("Number format error: " + nfe.getMessage());
		}
	}

	/**
	 * Attempts to read and parse a double delimited by whitespace.
	 */
	public double readDouble() throws IOException, ImportException {
		String token = readToken();
		try {
			return Double.parseDouble(token);
		} catch (NumberFormatException nfe) {
			throw new ImportException("Number format error: " + nfe.getMessage());
		}
	}

	/**
	 * Attempts to read and parse a double delimited by whitespace or by
	 * any character in delimiters.
	 */
	public double readDouble(String delimiters) throws IOException, ImportException {
		String token = readToken(delimiters);
		try {
			return Double.parseDouble(token);
		} catch (NumberFormatException nfe) {
			throw new ImportException("Number format error: " + nfe.getMessage());
		}
	}

	/**
	 * Reads a token stopping when any whitespace or a comment is found.
	 * If the token begins with a quote char then all characters will be
	 * included in token until a matching quote is found (including whitespace or comments).
	 */
	public String readToken() throws IOException {
		return readToken("");
	}

	/**
	 * Reads a token stopping when any whitespace, a comment or when any character
	 * in delimiters is found. If the token begins with a quote char
	 * then all characters will be included in token until a matching
	 * quote is found (including whitespace or comments).
	 */
	public String readToken(String delimiters) throws IOException {
		int space = 0;
		char ch, ch2, quoteChar = '\0';
		boolean done = false, first = true, quoted = false, isSpace;

		nextCharacter();

		StringBuffer token = new StringBuffer();

		while (!done) {
			ch = read();

			try {
				isSpace = Character.isWhitespace(ch);

				if (quoted && ch == quoteChar) { // Found the closing quote
					ch2 = read();

					if (ch == ch2) {
						// A repeated quote character so add this to the token
						token.append(ch);
					} else {
						// otherwise it terminates the token

						lastDelimiter = ' ';
						unreadCharacter(ch2);
						done = true;
						quoted = false;
					}
				} else if (first && (ch == '\'' || ch == '"')) {
					// if the opening character is a quote
					// read everything up to the closing quote
					quoted = true;
					quoteChar = ch;
					first = false;
					space = 0;
				} else if ( ch == startComment || ch == lineComment ) {
					skipComments(ch);
					lastDelimiter = ' ';
					done = true;
				} else {
					if (quoted) {
						// compress multiple spaces into one
						if (isSpace) {
							space++;
							ch = ' ';
						} else {
							space = 0;
						}

						if (space < 2) {
							token.append(ch);
						}
					} else if (isSpace) {
						lastDelimiter = ' ';
						done = true;
					} else if (delimiters.indexOf(ch) != -1) {
						done = true;
						lastDelimiter = ch;
					} else {
						token.append(ch);
						first = false;
					}
				}
			} catch (EOFException e) {
				// We catch an EOF and return the token we have so far
				done = true;
			}
		}

		if (Character.isWhitespace((char)lastDelimiter)) {
			ch = nextCharacter();
			while (Character.isWhitespace(ch)) {
				read();
				ch = nextCharacter();
			}

			if (delimiters.indexOf(ch) != -1) {
				lastDelimiter = readCharacter();
			}
		}

		return token.toString();
	}

	/**
	 * Skips over any comments. The opening comment delimiter is passed.
	 */
	protected void skipComments(char delimiter) throws IOException {

		char ch;
		int n=1;
		boolean write = false;
        StringBuffer meta = null;

		if (nextCharacter() == writeComment) {
			read();
			write = true;
		} else if (nextCharacter() == metaComment) {
			read();
            // combine two consecutive meta comments
            meta = lastMetaComment!= null ? new StringBuffer(lastMetaComment + ";") : new StringBuffer();
		}

        lastMetaComment = null;

		if (delimiter == lineComment) {
			String line = readLine();
			if (write && commentWriter != null) {
				commentWriter.write(line, 0, line.length());
				commentWriter.newLine();
			} else if (meta != null) {
				meta.append(line);
			}
		} else {
			do {
				ch = read();
				if (ch == startComment) {
					n++;
				} else if (ch == stopComment) {
					if (write && commentWriter != null) {
						commentWriter.newLine();
					}
					n--;
				} else if (write && commentWriter != null) {
					commentWriter.write(ch);
				} else if (meta != null) {
					meta.append(ch);
				}
			} while (n > 0);
		}

        if (meta != null) {
            lastMetaComment = meta.toString();
        }
	}

	/**
	 * Skips to the end of the line. If a comment is found then this is read.
	 */
	public void skipToEndOfLine() throws IOException {

		char ch;

		do {
			ch = read();
			if (hasComments) {
				if (ch == lineComment) {
					skipComments(ch);
					break;
				}
				if (ch == startComment) {
					skipComments(ch);
					ch = read();
				}
			}

		} while (ch != '\n' && ch != '\r');

		if (ch == '\r') {
			if (nextCharacter() == '\n') read();
		}
	}

	/**
	 * Skips char any contiguous characters in skip. Will also skip
	 * comments.
	 */
	public void skipWhile(String skip) throws IOException {

		char ch;

		do {
			ch = read();
		} while ( skip.indexOf(ch) > -1 );

		unreadCharacter(ch);
	}

	/**
	 * Skips over any space (plus tabs and returns) in the file. Will also skip
	 * comments.
	 */
	public void skipSpace() throws IOException {
		skipWhile(" \t\r\n");
	}

	/**
	 * Skips over any contiguous characters in skip. Will also skip
	 * comments and space.
	 */
	public void skipCharacters(String skip) throws IOException {
		skipWhile(skip + " \t\r\n");
	}

	/**
	 * Skips over the file until a character from delimiters is found. Returns
	 * the delimiter found. Will skip comments and will ignore delimiters within
	 * comments.
	 */
	public char skipUntil(String skip) throws IOException {
		char ch;

		do {
			ch = readCharacter();
		} while ( skip.indexOf(ch) == -1 );

		return ch;
	}

    public String getLastMetaComment() {
        return lastMetaComment;
    }

    public void clearLastMetaComment() {
        lastMetaComment = null;
    }

	// Private stuff

	private LineNumberReader reader;
	private BufferedWriter commentWriter = null;

	private int lastChar = '\0';
	private int lastDelimiter = '\0';

	private boolean hasComments = false;
	private char startComment = (char)-1;
	private char stopComment = (char)-1;
	private char lineComment = (char)-1;

    private char writeComment = (char)-1;
    private char metaComment = (char)-1;

    private String lastMetaComment = null;
}
