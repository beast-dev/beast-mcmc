package dr.geo;

import dr.app.gui.ColorFunction;

import javax.swing.*;
import java.awt.*;
import java.awt.geom.Rectangle2D;
import java.util.Random;

/**
 * @author Alexei Drummond
 */
public class InhomogeneousRandomWalk extends JComponent {

    Lattice lattice;
    RateMatrix rates;
    int i, j;
    Random random;
    double time;
    int[][] sample;
    int maxSample = 0;

    ColorFunction cf = new ColorFunction(
            new Color[]{Color.white, Color.blue, Color.magenta, Color.red},
            new float[]{0.0f, 0.10f, 0.20f, 1.0f}
    );

    public InhomogeneousRandomWalk(
            Lattice lattice, Location loc, Random random, RateMatrix rates) {

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

        lattice.paintLattice(g);

        for (int i = 0; i < lattice.latticeWidth(); i++) {
            for (int j = 0; j < lattice.latticeHeight(); j++) {
                if (sample[i][j] > 0) {
                    if (lattice.getState(i, j) >= 0) {
                        float intensity = (float) sample[i][j] / (float) maxSample;
                        g.setColor(cf.getColor(intensity));
                        g.drawRect(i, j, 1, 1);
                    }
                }
            }
        }
    }


    public double computeP(Location from, Location to, double t, int n) {


        int count = 0;
        for (int rep = 0; rep < n; rep++) {
            simulate(from, t);
            if (i == to.i && j == to.j) {
                count += 1;
            }
        }
        return (double) count / (double) n;
    }


    public void simulate(Location from, double t) {
        this.i = from.i;
        this.j = from.j;
        time = 0;
        while (time < t) {
            step(t);
        }
    }

    public void simulate(int steps) {
        for (int i = 0; i < steps; i++) {
            step(Double.MAX_VALUE);
        }
    }

    public void step(double toTime) {

        int fromState = lattice.getState(i, j);

        double moveUp = 0.0;
        double moveLeft = 0.0;
        double moveRight = 0.0;
        double moveDown = 0.0;

        if (i > 0) {
            moveLeft = rates.getRate(fromState, lattice.getState(i - 1, j));
        }

        if (i < lattice.latticeWidth() - 1) {
            moveRight = rates.getRate(fromState, lattice.getState(i + 1, j));
        }

        if (j > 0) {
            moveUp = rates.getRate(fromState, lattice.getState(i, j - 1));
        }

        if (j < lattice.latticeHeight() - 1) {
            moveDown = rates.getRate(fromState, lattice.getState(i, j + 1));
        }

        double totalRate =
                moveUp + moveDown + moveLeft + moveRight;


        double V = random.nextDouble();
        double dt = -Math.log(V) / totalRate;
        if (time + dt > toTime) {
            time = toTime;
            return;
        }

        double U = random.nextDouble() * totalRate;

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
        if (sample[i][j] > maxSample) maxSample = sample[i][j];
    }

    public static void main(String[] args) {

        String kmlFile = args[0];
        double left = Double.parseDouble(args[1]);
        double top = Double.parseDouble(args[2]);
        double width = Double.parseDouble(args[3]);
        double height = Double.parseDouble(args[4]);

        Rectangle2D bounds = new Rectangle2D.Double(left, top, width, height);

        Random random = new Random();


        final double[][] rates = new double[2][2];
        rates[0][0] = 0.0;
        rates[0][1] = 0.0;
        rates[1][0] = 0.0;
        rates[1][1] = 2.0;

        KMLRenderer renderer = new KMLRenderer(kmlFile, Color.white, Color.black);
        renderer.setBounds(bounds);
        renderer.render(1000);

        RateMatrix matrix = new RateMatrix() {

            public double getRate(int i, int j) {
                return rates[i][j];
            }
        };

        InhomogeneousRandomWalk walk = new InhomogeneousRandomWalk(renderer, new Location(0, 0), random, matrix);

        for (int i = 0; i < 5; i++) {
            Location start = Lattice.Utils.getRandomLocation(renderer, 1, random);
            for (int j = 0; j < 1000; j++) {
                walk.simulate(start, 4000);
            }
        }
        JFrame frame = new JFrame();

        frame.getContentPane().add(BorderLayout.CENTER, walk);
        frame.pack();
        frame.setVisible(true);
    }


}
