/*
 * TreeImporter.java
 *
 * Copyright © 2002-2024 the BEAST Development Team
 * http://beast.community/about
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
 *
 */

package dr.evolution.io;

import dr.evolution.tree.Tree;
import dr.evolution.util.TaxonList;

import java.io.EOFException;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

/**
 * Interface for importers that do trees
 *
 *
 * @author Andrew Rambaut
 * @author Alexei Drummond
 */
public interface TreeImporter {

	/**
	 * Counts the number of trees in the file without importing them. The file
	 * will need to be reset to the beginning to read the trees.
	 *
	 * @return the number of trees
     */
	int countTrees() throws IOException;

	/**
	 * return whether another tree is available. 
	 */
	boolean hasTree() throws IOException, Importer.ImportException;

	/**
	 * import the next tree. 
	 * return the tree or null if no more trees are available
	 */
	Tree importNextTree() throws IOException, Importer.ImportException;

	/**
	 * import a single tree. 
	 */
	Tree importTree(TaxonList taxonList) throws IOException, Importer.ImportException;

	/**
	 * import an array of all trees. 
	 */
	List<Tree> importTrees(TaxonList taxonList) throws IOException, Importer.ImportException;
}
