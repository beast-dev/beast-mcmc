/*
 * Lattice.java
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

package dr.geo;

import java.awt.*;
import java.util.Random;

/**
 * @author Alexei Drummond
 */
public interface Lattice {

    int latticeWidth();

    int latticeHeight();

    int getState(int i, int j);

    public void paintLattice(Graphics g);

    class Utils {

        public static Location getRandomLocation(Lattice lattice, int state, Random random) {

            int x = random.nextInt(lattice.latticeWidth());
            int y = random.nextInt(lattice.latticeHeight());

            while (lattice.getState(x, y) != state) {
                x = random.nextInt(lattice.latticeWidth());
                y = random.nextInt(lattice.latticeHeight());
            }
            return new Location(x, y);
        }
    }
}
