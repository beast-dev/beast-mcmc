/*
 * NexusApplicationImporter.java
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

package dr.app.oldbeauti;

import dr.evolution.io.NexusImporter;

import java.io.IOException;
import java.io.Reader;
import java.io.Writer;

/**
 * Class for importing PAUP, MrBayes and Rhino NEXUS file format
 *
 * @author			Andrew Rambaut
 * @author			Alexei Drummond
 * @version			$Id: NexusApplicationImporter.java,v 1.4 2005/07/11 14:07:25 rambaut Exp $
 */
public class NexusApplicationImporter extends NexusImporter { 

	public static final NexusBlock PAUP_BLOCK = new NexusBlock("PAUP");
	public static final NexusBlock MRBAYES_BLOCK = new NexusBlock("MRBAYES");
	public static final NexusBlock RHINO_BLOCK = new NexusBlock("RHINO");
	public static final NexusBlock BEAST_BLOCK = new NexusBlock("BEAST");

	/**
	 * Constructor
	 */
	public NexusApplicationImporter(Reader reader) {
		super(reader);
		setCommentDelimiters('[', ']', '\0');
	}
	
	public NexusApplicationImporter(Reader reader, Writer commentWriter) {
		super(reader, commentWriter);
		setCommentDelimiters('[', ']', '\0');
	}
	
	/**
	 * This function returns an enum class to specify what the
	 * block given by blockName is. 
	 */
	public NexusBlock findBlockName(String blockName)
	{
		if (blockName.equalsIgnoreCase(PAUP_BLOCK.toString())) {
			return PAUP_BLOCK;
		} else  if (blockName.equalsIgnoreCase(MRBAYES_BLOCK.toString())) {
			return MRBAYES_BLOCK;
		} else  if (blockName.equalsIgnoreCase(BEAST_BLOCK.toString())) {
			return BEAST_BLOCK;
		} else  if (blockName.equalsIgnoreCase(RHINO_BLOCK.toString())) {
			return RHINO_BLOCK;
		} else  {
			return super.findBlockName(blockName);
		}
	}

	/**
	 * Parses a 'PAUP' block. 
	 */
	public void parsePAUPBlock(BeastGenerator options) throws ImportException, IOException
	{
		// PAUP is largely a subset of BEAST block
		readBEASTBlock(options);
	}

	/**
	 * Parses a 'MRBAYES' block. 
	 */
	public void parseMrBayesBlock(BeastGenerator options) throws ImportException, IOException
	{
		// MRBAYES is largely a subset of BEAST block
		readBEASTBlock(options);
	}

	/**
	 * Parses a 'BEAST' block. 
	 */
	public void parseBEASTBlock(BeastGenerator options) throws ImportException, IOException
	{
		readBEASTBlock(options);
	}

	/**
	 * Parses a 'RHINO' block. 
	 */
	public void parseRhinoBlock(BeastGenerator options) throws ImportException, IOException
	{
		readRhinoBlock(options);
	}

	private void readBEASTBlock(BeastGenerator options) throws ImportException, IOException
	{
		boolean done = false;
				
		while (!done) {
			String command = readToken(";");
			if (match("HSEARCH", command, 2)) {
				done = true;
			} else if (match("MCMC", command, 4)) {
				if (getLastDelimiter() != ';') {
					readMCMCCommand(options);
				}
				done = true;
			} else if (match("MCMCP", command, 5)) {
				if (getLastDelimiter() != ';') {
					readMCMCCommand(options);
				}
			} else if (match("LSET", command, 2)) {
				if (getLastDelimiter() != ';') {
					readLSETCommand(options);
				}
			} else if (command.equalsIgnoreCase("ENDBLOCK") || command.equalsIgnoreCase("END")) {
				done = true;
			} else {
						
				System.err.println("The command, '" + command + "', is not used by BEAST and has been ignored");
			}
		}
	}

	private void readLSETCommand(BeautiOptions options) throws ImportException, IOException
	{
		boolean done = false;
				
		while (!done) {
			String subcommand = readToken("=;");
			if (match("NST", subcommand, 2)) {
				int nst = readInteger( ";" );
				if (nst == 1) {
					options.nucSubstitutionModel = BeautiOptions.JC;
				} else if (nst == 2) {
					options.nucSubstitutionModel = BeautiOptions.HKY;
				} else if (nst == 6) {
					options.nucSubstitutionModel = BeautiOptions.GTR;
				} else {
					throw new BadFormatException("Bad value for NST subcommand of LSET command");
				}
			} else if (match("RATES", subcommand, 2)) {
				String token = readToken( ";" );
				
				if (match("EQUAL", token, 1)) {
					options.gammaHetero = false;
					options.invarHetero = false;
				} else if (match("GAMMA", token, 1)) {
					options.gammaHetero = true;
					options.invarHetero = false;
				} else if (match("PROPINV", token, 1)) {
					options.gammaHetero = false;
					options.invarHetero = true;
				} else if (match("INVGAMMA", token, 1)) {
					options.gammaHetero = true;
					options.invarHetero = true;
				} else if (match("ADGAMMA", token, 1)) {
					System.err.println("The option, 'RATES=ADGAMMA', in the LSET command is not used by BEAST and has been ignored");
				} else if (match("SITESPEC", token, 1)) {
					System.err.println("The option, 'RATES=SITESPEC', in the LSET command is not used by BEAST and has been ignored");
				} else {
					throw new BadFormatException("Unknown value, '" + token + "'");
				}
			} else if (match("NGAMMACAT", subcommand, 2)) {
			
				options.gammaCategories = readInteger( ";" );
			} else {
						
				System.err.println("The option, '" + subcommand + "', in the LSET command is not used by BEAST and has been ignored");
			}
			
			if (getLastDelimiter() == ';') {
				done = true;
			}
		}
	}

	private void readMCMCCommand(BeautiOptions options) throws ImportException, IOException
	{
		boolean done = false;
				
		while (!done) {
			String subcommand = readToken("=;");
			if (match("NGEN", subcommand, 2)) {
				options.chainLength = readInteger( ";" );
			} else if (match("SAMPLEFREQ", subcommand, 2)) {
				options.logEvery = readInteger( ";" );
			} else if (match("PRINTFREQ", subcommand, 1)) {
				options.echoEvery = readInteger( ";" );
			} else if (match("FILENAME", subcommand, 1)) {
				options.fileName = readToken( ";" );
			} else if (match("BURNIN", subcommand, 1)) {
				options.burnIn = readInteger( ";" );
			} else if (match("STARTINGTREE", subcommand, 2)) {
				String token = readToken(";");
				if (match("USER", token, 1)) {
					options.userTree = true;
				} else if (match("RANDOM", token, 1)) {
					options.userTree = false;
				} else {
					throw new BadFormatException("Unknown value, '" + token + "'");
				}
			} else {
						
				System.err.println("The option, '" + subcommand + "', in the MCMC command is not used by BEAST and has been ignored");
			}
			
			if (getLastDelimiter() == ';') {
				done = true;
			}
		}
	}

	private void readRhinoBlock(BeastGenerator options) throws ImportException, IOException
	{
		boolean done = false;
				
		while (!done) {
			String command = readToken(";");
			if (match("NUCMODEL", command, 2)) {
				if (getLastDelimiter() != ';') {
					readNUCMODELCommand(options);
				}
			} else if (match("SITEMODEL", command, 2)) {
				if (getLastDelimiter() != ';') {
					readSITEMODELCommand(options);
				}
			} else if (match("TREEMODEL", command, 2)) {
				if (getLastDelimiter() != ';') {
					readTREEMODELCommand(options);
				}
			} else if (match("CPPARTITIONMODEL", command, 2)) {
				if (getLastDelimiter() != ';') {
					readCPPARTITIONMODELCommand(options);
				}
			} else if (match("OPTIMIZE", command, 1)) {
				done = true;
			} else if (command.equalsIgnoreCase("ENDBLOCK") || command.equalsIgnoreCase("END")) {
				done = true;
			} else {
						
				System.err.println("The command, '" + command + "', is not used by BEAST and has been ignored");
			}
		}
	}

	private void readNUCMODELCommand(BeautiOptions options) throws ImportException, IOException
	{
		boolean done = false;
				
		while (!done) {
			String subcommand = readToken("=;");
			if (match("TYPE", subcommand, 1)) {
				String token = readToken(";");
				if (match("HKY", token, 1)) {
					options.nucSubstitutionModel = BeautiOptions.HKY;
				} else if (match("GTR", token, 1)) {
					options.nucSubstitutionModel = BeautiOptions.GTR;
				} else if (match("F84", token, 1)) {
					System.err.println("The option, 'TYPE=F84', in the NUCMODEL command is not used by BEAST and has been ignored");
				} else {
					throw new BadFormatException("Unknown value, '" + token + "'");
				}
			} else {
						
				System.err.println("The option, '" + subcommand + "', in the NUCMODEL command is not used by BEAST and has been ignored");
			}
			
			if (getLastDelimiter() == ';') {
				done = true;
			}
		}
	}

	private void readSITEMODELCommand(BeautiOptions options) throws ImportException, IOException
	{
		boolean done = false;
				
		while (!done) {
			String subcommand = readToken("=;");
			if (match("TYPE", subcommand, 1)) {
				String token = readToken(";");
				if (match("HOMOGENEOUS", token, 1)) {
					options.gammaHetero = false;
					options.invarHetero = false;
				} else if (match("GAMMA", token, 2)) {
					options.gammaHetero = true;
					options.invarHetero = false;
				} else if (match("INVAR", token, 1)) {
					options.gammaHetero = false;
					options.invarHetero = true;
				} else if (match("GI", token, 2)) {
					options.gammaHetero = true;
					options.invarHetero = true;
				} else {
					throw new BadFormatException("Unknown value, '" + token + "'");
				}
			} else if (match("NUMCAT", subcommand, 1)) {
				options.gammaCategories = readInteger( ";" );
			} else {
						
				System.err.println("The option, '" + subcommand + "', in the SITEMODEL command is not used by BEAST and has been ignored");
			}
			
			if (getLastDelimiter() == ';') {
				done = true;
			}
		}
	}

	private void readTREEMODELCommand(BeautiOptions options) throws ImportException, IOException
	{
		boolean done = false;
				
		while (!done) {
			String subcommand = readToken("=;");
			if (match("TYPE", subcommand, 1)) {
				String token = readToken(";");
				if (match("UNCONSTRAINED", token, 1)) {
					System.err.println("The option, 'TYPE=UNCONSTRAINED', in the TREEMODEL command is not used by BEAST and has been ignored");
				} else if (match("CONSTRAINED", token, 1)) {
					// do nothing
				} else if (match("NODEDATES", token, 1)) {
					// do nothing
				} else if (match("TIPDATES", token, 1)) {
					// do nothing
				} else {
					throw new BadFormatException("Unknown value, '" + token + "'");
				}
			} else {
						
				System.err.println("The option, '" + subcommand + "', in the TREEMODEL command is not used by BEAST and has been ignored");
			}
			
			if (getLastDelimiter() == ';') {
				done = true;
			}
		}
	}

	private void readCPPARTITIONMODELCommand(BeautiOptions options) throws ImportException, IOException
	{
		boolean done = false;
		
		options.codonHeteroPattern = null;

		while (!done) {
			String subcommand = readToken("=;");
			if (match("ON", subcommand, 1)) {
				String token = readToken(";");
				if (match("TRUE", token, 1)) {
					options.codonHeteroPattern = "123";
				} else if (match("FALSE", token, 1)) {
					options.codonHeteroPattern = null;
				} else {
					throw new BadFormatException("Unknown value, '" + token + "'");
				}
			} else if (match("SUBSTMODEL", subcommand, 1)) {
				String token = readToken(";");
				if (match("TRUE", token, 1)) {
					options.unlinkedSubstitutionModel = true;
					options.unlinkedHeterogeneityModel = true;
				} else if (match("FALSE", token, 1)) {
					options.unlinkedSubstitutionModel = false;
					options.unlinkedHeterogeneityModel = false;
				} else {
					throw new BadFormatException("Unknown value, '" + token + "'");
				}
			} else if (match("FREQMODEL", subcommand, 1)) {
				String token = readToken(";");
				if (match("TRUE", token, 1)) {
					options.unlinkedFrequencyModel = true;
				} else if (match("FALSE", token, 1)) {
					options.unlinkedFrequencyModel = false;
				} else {
					throw new BadFormatException("Unknown value, '" + token + "'");
				}
			} else {
						
				System.err.println("The option, '" + subcommand + "', in the CPPARTITIONMODEL command is not used by BEAST and has been ignored");
			}
			
			if (getLastDelimiter() == ';') {
				done = true;
			}
		}
	}

	private boolean match(String reference, String target, int min) throws ImportException
	{
		if (target.length() < min) {
			throw new BadFormatException("Ambiguous command or subcommand, '" + target + "'");
		}
		
		return reference.startsWith(target.toUpperCase());
	}
}
