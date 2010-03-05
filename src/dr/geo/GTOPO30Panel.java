/*
 * GTOPO30Panel.java
 *
 * Copyright (C) 2002-2010 Alexei Drummond and Andrew Rambaut
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
 * BEAST is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with BEAST; if not, write to the
 * Free Software Foundation, Inc., 51 Franklin St, Fifth Floor,
 * Boston, MA  02110-1301  USA
 */

package dr.geo;

import dr.gui.ColorFunction;

import javax.swing.*;
import java.awt.*;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

/**
 * @author Alexei Drummond
 */
public class GTOPO30Panel extends JPanel {

    GTOPO30Tile[][] tiles;

    int minLat, maxLat;
    int minLong, maxLong;

    double scale = 0.5;

    public GTOPO30Panel(String[] filenames, ColorFunction colorFunction) throws IOException {

        List<GTOPO30Tile> tileList = new ArrayList<GTOPO30Tile>();

        for (String filename : filenames) {
            GTOPO30Tile tile = new GTOPO30Tile(filename, colorFunction);

            if (tile.getMaxLatitude() > maxLat) maxLat = tile.getMaxLatitude();
            if (tile.getMaxLongitude() > maxLong) maxLong = tile.getMaxLongitude();
            if (tile.getMinLatitude() < minLat) minLat = tile.getMinLatitude();
            if (tile.getMinLongitude() < minLong) minLong = tile.getMinLongitude();

            tileList.add(tile);
        }

        tiles = new GTOPO30Tile[(maxLat - minLat) / 50][(maxLong - minLong) / 40];

        for (GTOPO30Tile tile : tileList) {
            tiles[(maxLat - tile.getMaxLatitude()) / 50][(tile.getMinLongitude() - minLong) / 40] = tile;
        }

        setLayout(new GridLayout(tiles.length, tiles[0].length));

        for (GTOPO30Tile[] tileRow : tiles) {
            for (GTOPO30Tile tile : tileRow) {
                if (tile != null) {
                    add(tile);
                } else {
                    add(new JPanel());
                }
            }
        }
    }

    public Dimension getPreferredSize() {

        return new Dimension(
                (int) Math.round(tiles[0].length * GTOPO30Tile.NCOLS * scale),
                (int) Math.round(tiles.length * GTOPO30Tile.NROWS * scale));

    }

    /**
     * @param y latitudinal pixel, increasing south from north edge
     * @param x longitudinal pixel index, increasing east from west edge
     * @return the height of the pixel at y,x
     */
    public short getHeight(int y, int x) {
        int tileY = y / GTOPO30Tile.NCOLS;
        int tileX = x / GTOPO30Tile.NROWS;

        int ty = y % GTOPO30Tile.NCOLS;
        int tx = x % GTOPO30Tile.NROWS;

        return tiles[tileY][tileX].getHeight(ty, tx);
    }

    public void setScale(double scale) {
        this.scale = scale;
        setSize(getPreferredSize());
        repaint();
    }
}
