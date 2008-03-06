/*
 * DistanceMatrix.java
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

package dr.evolution.distance;

import dr.evolution.alignment.PatternList;
import dr.evolution.datatype.DataType;
import dr.evolution.util.Taxon;
import dr.evolution.util.TaxonList;
import dr.matrix.Matrix;

/**
 * storage for pairwise distance matrices.<p>
 *
 * @version $Id: DistanceMatrix.java,v 1.23 2005/07/11 14:06:25 rambaut Exp $
 *
 * @author Andrew Rambaut
 * @author Alexei Drummond
 */
public class DistanceMatrix extends Matrix.AbstractMatrix implements TaxonList {

	public static final double MAX_DISTANCE = 1000.0;

	/** constructor */
	public DistanceMatrix()
	{
		super();
	}

	/** constructor taking a dimension */
	public DistanceMatrix(TaxonList taxa)
	{
		super();
		this.taxa = taxa;
		dimension = taxa.getTaxonCount();
		distances = new double[dimension][dimension];
		distancesKnown = true;
	}

	/** constructor taking a pattern source */
	public DistanceMatrix(PatternList patterns)
	{
		super();
		setPatterns(patterns);
	}

	/**
	 * set the pattern source
	 */
	public void setPatterns(PatternList patterns) {
		this.taxa = patterns;
		this.patterns = patterns;
		dimension = patterns.getTaxonCount();
		dataType = patterns.getDataType();
		distancesKnown = false;
	}

	/**
	 * @return the number of rows
	 */
	public int getRowCount() { return dimension; }

	/**
	 * @return the number of columns
	 */
	public int getColumnCount() { return dimension; }

	/**
	 * @return an element
	 */
	public double getElement(int row, int column) {

		if (!distancesKnown) {
			calculateDistances();
		}

		return distances[row][column];
	}

	/**
	 * set an element - this overwrites any existing elements
	 */
	public void setElement(int row, int column, double value) {

		if (!distancesKnown) {
			calculateDistances();
		}

		distances[row][column] = value;
	}

	/**
	 * Calculate the distances
	 */
	public void calculateDistances() {
		distances = new double[dimension][dimension];

		for (int i = 0; i < dimension; i++) {
			for (int j = i + 1; j < dimension; j++) {
				distances[i][j] = calculatePairwiseDistance(i, j);
				distances[j][i] = distances[i][j];
			}
			distances[i][i] = 0.0;
		}

		distancesKnown = true;
	}

	/**
	 * Calculate a pairwise distance
	 */
	protected double calculatePairwiseDistance(int taxon1, int taxon2) {
		int state1, state2;

		int n = patterns.getPatternCount();
		double weight, distance;
		double sumDistance = 0.0;
		double sumWeight = 0.0;

		int[] pattern;

		for (int i = 0; i < n; i++) {
			pattern = patterns.getPattern(i);

			state1 = pattern[taxon1];
			state2 = pattern[taxon2];

			weight = patterns.getPatternWeight(i);
//			sumDistance += dataType.getObservedDistance(state1, state2) * weight;
			if (!dataType.isAmbiguousState(state1) && !dataType.isAmbiguousState(state2) &&
					state1 != state2) {
				sumDistance += weight;
			}
			sumWeight += weight;
		}

		distance = sumDistance / sumWeight;

		return distance;
	}


	/**
	 * Returns the mean pairwise distance of this matrix
	 */
	public double getMeanDistance() {
		if (!distancesKnown) {
		    calculateDistances();
		}

		double dist = 0.0;
		int count = 0;
		for (int i = 0; i < dimension; i++) {
			for (int j = 0; j < dimension; j++) {
				if (i != j) {
					dist += distances[i][j];
					count += 1;
				}
			}
		}
		return dist / (double)count;
	}

	public String toString() {
		try {
			double[] dists = getUpperTriangle();
            StringBuffer buffer = new StringBuffer(String.valueOf(dists[0]));

            for (int i = 1; i < dists.length; i++) {
                buffer.append(", ").append(String.valueOf(dists[i]));
            }
            return buffer.toString();
        } catch(Matrix.NotSquareException e) {
            return e.toString();
        }
    }

    // **************************************************************
    // TaxonList IMPLEMENTATION
    // **************************************************************

	/**
	 * @return a count of the number of taxa in the list.
	 */
	public int getTaxonCount() {
		return taxa.getTaxonCount();
	}

	/**
	 * @return the ith taxon.
	 */
	public Taxon getTaxon(int taxonIndex) {
		return taxa.getTaxon(taxonIndex);
	}

	/**
	 * @return the ID of the ith taxon.
	 */
	public String getTaxonId(int taxonIndex) {
		return taxa.getTaxonId(taxonIndex);
	}

	/**
	 * returns the index of the taxon with the given id.
	 */
	public int getTaxonIndex(String id) {
		return taxa.getTaxonIndex(id);
	}

	/**
	 * returns the index of the given taxon.
	 */
	public int getTaxonIndex(Taxon taxon) {
		return taxa.getTaxonIndex(taxon);
	}

	/**
	 * @return an object representing the named attributed for the given taxon.
	 * @param taxonIndex the index of the taxon whose attribute is being fetched.
	 * @param name the name of the attribute of interest.
	 */
	public Object getTaxonAttribute(int taxonIndex, String name) {
		return taxa.getTaxonAttribute(taxonIndex, name);
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

	/*
	public static XMLObjectParser PARSER = new XMLObjectParser() {

		public String getParserName() { return DISTANCE_MATRIX; }

		public Object parseXMLObject(XMLObject xo, ObjectStore store) throws XMLParseException {

			PatternList patterns = null;
			Taxa taxa = new Taxa();
			ArrayList rows = new ArrayList();

			for (int i = 0; i < xo.getChildCount(); i++) {

				Object child = xo.getChild(i);
				if (child instanceof XMLObject) {

					XMLObject cxo = (XMLObject) child;
					if (cxo.getName().equals(XMLElement.ROW)) {
						if (patterns != null) {
							throw new XMLParseException("distances element must contain either patterns or rows");
						}

						Taxon taxon = null;

						Variate row = new Variate.Double();

						for (int j = 0; j < cxo.getChildCount(); j++) {

							Object child2 = cxo.getChild(j);
							if (child2 instanceof Taxon) {

								taxon = (Taxon)child2;

							} else if (child2 instanceof String) {

								StringTokenizer tokens = new StringTokenizer( (String)child2 );
								while (tokens.hasMoreElements()) {
									String token = tokens.nextToken();
									row.add(Double.parseDouble(token));
								}
							} else throw new XMLParseException("unknown child element found in row element in distances element");
						}
						taxa.addTaxon(taxon);
						rows.add(row);

					} else throw new XMLParseException("unknown child element found in distances element");

				} else if (child instanceof PatternList) {
					if (patterns != null) {
						throw new XMLParseException("multiple pattern list elements in distances element");
					}
					if (rows.size() > 0) {
						throw new XMLParseException("distances element must contain either patterns or rows");
					}

					patterns = (PatternList)child;
				} else throw new XMLParseException("unknown child element found in distances element");
			}

			DistanceMatrix matrix = null;

			if (rows.size() > 0) {

				matrix = new DistanceMatrix(taxa);
				Variate row = null;

				for (int i = 0; i < rows.size(); i++) {

					row = (Variate)rows.get(i);
					if (row.getCount() != rows.size()) {
						throw new XMLParseException("distances element must be square (no. rows = no. columns)");
					}

					for (int j = 0; j < row.getCount(); j++) {
						matrix.setElement(i, j, row.get(j));
					}
				}
			} else if (patterns != null) {

				if (xo.hasAttribute(TYPE)) {
					String type = xo.getStringAttribute(TYPE);
					if (type.equals(Nucleotides.JC)) {
						System.out.print("Creating jukes-cantor distance matrix...");
						System.out.flush();
						matrix = new JukesCantorDistanceMatrix(patterns);
					} else if (type.equals(Nucleotides.F84)) {
						System.out.print("Creating F84 distance matrix...");
						System.out.flush();
						matrix = new F84DistanceMatrix(patterns);
					} else
						throw new XMLParseException("unknown type attribute in distances element");
				} else {
					matrix = new DistanceMatrix(patterns);
				}
			} else {
				throw new XMLParseException("pattern list element (patterns or alignment) or row elements missing from distances element");
			}
			System.out.println("done.");


			return matrix;
		}
	};*/

	//
	// Private stuff
	//

	protected DataType dataType = null;
	int dimension = 0;
	boolean distancesKnown;
	private double[][] distances = null;
	protected PatternList patterns = null;
	private TaxonList taxa = null;
}
