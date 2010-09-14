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
