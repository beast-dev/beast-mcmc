/*
 * IslandLattice.java
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

import javax.swing.*;
import java.awt.*;
import java.util.*;
import java.util.List;

/**
 * @author Alexei Drummond
 */
public class IslandLattice extends JComponent implements Lattice {

    private int[][] lattice;
    private final int width;
    private final int height;

    public IslandLattice(int width, int height) {
        lattice = new int[width][height];
        this.width = width;
        this.height = height;
    }

    public Dimension getPreferredSize() {
        return new Dimension(width, height);
    }

    public void paintLattice(Graphics g) {
        g.setColor(Color.blue);
        g.fillRect(0, 0, width, height);
        g.setColor(Color.green);

        for (int i = 0; i < lattice.length; i++) {
            for (int j = 0; j < lattice.length; j++) {
                if (lattice[i][j] > 0) {
                    g.drawRect(i, j, 1, 1);
                }
            }
        }

    }

    public void paintComponent(Graphics g) {
        paintLattice(g);
    }

    public final int latticeWidth() {
        return width;
    }

    public final int latticeHeight() {
        return height;
    }

    public void setState(int i, int j, int state) {
        lattice[i][j] = state;
    }

    public int getState(int i, int j) {
        return lattice[i][j];
    }

    public int getState(Location loc) {
        return lattice[loc.i][loc.j];
    }


    public void smooth(int boundary) {

        int[][] newLattice = new int[width][height];

        for (int i = 0; i < lattice.length; i++) {
            for (int j = 0; j < lattice.length; j++) {

                int ncount = getNeighbourSum(i, j);

                if (ncount > boundary) {
                    newLattice[i][j] = 1;
                } else {
                    newLattice[i][j] = 0;
                }
            }
        }
        lattice = newLattice;
    }


    /**
     * @param n the number of pixels to be flipped
     * @return the locations of unique islands formed, since some may merge
     */
    public List<Location> addIslands(int n, Random random) {

        List<Location> locations = new ArrayList<Location>();

        for (int k = 0; k < n; k++) {

            int i = random.nextInt(width);
            int j = random.nextInt(height);

            if (getState(i, j) == 0) {

                setState(i, j, 1);
                if (getNeighbourSum(i, j) == 0) {
                    locations.add(new Location(i, j));
                }
            }
        }
        return locations;
    }

    /**
     * @param n      the number of pixels to move from 0 -> 1
     * @param random
     * @return the number of attempts taken to grow the current islands the
     *         perscribed amount of land
     */
    public int growIslands(int n, Random random) {
        int count = 0;
        int attempts = 0;
        while (count < n) {
            int i = random.nextInt(width);
            int j = random.nextInt(height);

            if (getState(i, j) == 0) {

                if (getNeighbourSum(i, j) > 0) {
                    setState(i, j, 1);
                    count += 1;
                }
            }
            attempts += 1;

        }
        return attempts;
    }


    public int sum() {
        int sum = 0;
        for (int i = 0; i < lattice.length; i++) {
            for (int j = 0; j < lattice.length; j++) {

                sum += lattice[i][j];
            }
        }
        return sum;
    }


    int getNeighbourSum(int i, int j) {

        int sum = 0;
        if (i > 0) {
            sum += getState(i - 1, j);
        }
        if (i < width - 1) {
            sum += getState(i + 1, j);
        }

        if (j > 0) {
            sum += getState(i, j - 1);
        }
        if (j < height - 1) {
            sum += getState(i, j + 1);
        }
        return sum;
    }

    private List<Location> getNeighbours(Location loc) {

        List<Location> neighbours = new ArrayList<Location>();

        if (loc.i > 0) {
            neighbours.add(new Location(loc.i - 1, loc.j));
        }
        if (loc.i < width - 1) {
            neighbours.add(new Location(loc.i + 1, loc.j));
        }
        if (loc.j > 0) {
            neighbours.add(new Location(loc.i, loc.j - 1));
        }
        if (loc.j < height - 1) {
            neighbours.add(new Location(loc.i, loc.j + 1));
        }
        return neighbours;
    }


    private Set<Location> getEmptyNeighboursOfIsland(Location loc) {
        List<Location> toProcess = new ArrayList<Location>();
        Set<Location> processed = new HashSet<Location>();
        Set<Location> empty = new HashSet<Location>();
        toProcess.add(loc);

        while (toProcess.size() > 0) {
            Location l = toProcess.remove(toProcess.size() - 1);
            List<Location> neighbours = getNeighbours(l);
            for (Location neighbour : neighbours) {
                if (getState(neighbour) == 0) {
                    empty.add(neighbour);
                } else if (!processed.contains(neighbour)) {
                    toProcess.add(neighbour);
                }
            }
            processed.add(l);
        }

        return empty;
    }

    public static void main(String[] args) {

        Random random = new Random();

        IslandLattice lattice = new IslandLattice(500, 500);
        List<Location> islands = lattice.addIslands(50, random);
        System.out.println("Created " + islands.size() + " islands");

        int growth = 4950;

        int attempts = lattice.growIslands(growth, random);
        System.out.println("Took " + attempts + " attempts to grow islands by " + growth);

        lattice.smooth(2);

        JFrame frame = new JFrame();
        frame.getContentPane().add(lattice, BorderLayout.CENTER);
        frame.pack();
        frame.setVisible(true);

    }
}
