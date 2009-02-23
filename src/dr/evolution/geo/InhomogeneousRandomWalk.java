package dr.evolution.geo;

import javax.swing.*;
import java.awt.*;
import java.util.Random;

/**
 * @author Alexei Drummond
 */
public class InhomogeneousRandomWalk extends JComponent {

    Lattice lattice;
    double[][] rates;
    int i, j;
    Random random;
    double time;
    int[][] sample;

    public InhomogeneousRandomWalk(
            Lattice lattice, Location loc, Random random, double[][] rates) {

        this.lattice = lattice;
        this.rates = rates;
        this.i = loc.i;
        this.j = loc.j;
        this.random = random;

        sample = new int[lattice.latticeWidth()][lattice.latticeHeight()];

    }

    public Dimension getPreferredSize() {
        return new Dimension(lattice.latticeWidth(), lattice.latticeHeight());
    }

    public void paintComponent(Graphics g) {

        lattice.paintComponent(g);

        g.setColor(Color.black);

        for (int i = 0; i < lattice.latticeWidth(); i++) {
            for (int j = 0; j < lattice.latticeHeight(); j++) {
                if (sample[i][j] > 0) {
                    if (lattice.getState(i, j) == 0) {
                        g.setColor(Color.black);
                    } else {
                        g.setColor(Color.red);
                    }
                    g.drawRect(i, j, 1, 1);
                }
            }
        }
    }


    public void simulate(int steps) {
        for (int i = 0; i < steps; i++) {
            step();
        }
    }

    public void step() {

        int fromState = lattice.getState(i, j);

        double moveUp = 0.0;
        double moveLeft = 0.0;
        double moveRight = 0.0;
        double moveDown = 0.0;

        if (i > 0) {
            moveLeft = rates[fromState][lattice.getState(i - 1, j)];
        }

        if (i < lattice.latticeWidth() - 1) {
            moveRight = rates[fromState][lattice.getState(i + 1, j)];
        }

        if (j > 0) {
            moveUp = rates[fromState][lattice.getState(i, j - 1)];
        }

        if (j < lattice.latticeHeight() - 1) {
            moveDown = rates[fromState][lattice.getState(i, j + 1)];
        }

        double totalRate =
                moveUp + moveDown + moveLeft + moveRight;

        double U = random.nextDouble() * totalRate;
        double V = random.nextDouble();
        double dt = -Math.log(V) / totalRate;

        if (U < moveUp) {
            // move up
            j -= 1;
        } else {
            U -= moveUp;
            if (U < moveRight) {
                // move right
                i += 1;
            } else {
                U -= moveRight;
                if (U < moveDown) {
                    // move down
                    j += 1;
                } else {
                    // move left
                    i -= 1;
                }
            }
        }
        time += dt;
        sample[i][j] += 1;
    }

    public static void main(String[] args) {

        Random random = new Random();

        Lattice lattice = new Lattice(500, 500);
        java.util.List<Location> islands = lattice.addIslands(50, random);
        System.out.println("Created " + islands.size() + " islands");

        int growth = 4950;

        int attempts = lattice.growIslands(growth, random);
        System.out.println("Took " + attempts + " attempts to grow islands by " + growth);

        lattice.smooth(2);

        double[][] rates = new double[2][2];
        rates[0][0] = 1.0;
        rates[0][1] = 20.0;
        rates[1][0] = 0.05;
        rates[1][1] = 0.5;

        Location loc = new Location(250, 250);

        InhomogeneousRandomWalk walk = new InhomogeneousRandomWalk(lattice, loc, random, rates);
        walk.simulate(100000000);

        JFrame frame = new JFrame();
        frame.getContentPane().add(walk, BorderLayout.CENTER);
        frame.pack();
        frame.setVisible(true);

    }


}
