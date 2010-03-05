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

    public static void main(String[] args) throws IOException {

        ColorFunction function = new ColorFunction(
                new Color[]{Color.blue, Color.yellow, Color.green.darker(), Color.orange.darker(), Color.white, Color.pink},
                new float[]{-410, 0, 100, 1500, 4000, 8800});

        GTOPO30Panel panel = new GTOPO30Panel(args, function);

        JFrame frame = new JFrame("GTOPO30");
        frame.getContentPane().add(BorderLayout.CENTER, panel);
        frame.setSize(1200, 1000);
        frame.setVisible(true);
    }
}
