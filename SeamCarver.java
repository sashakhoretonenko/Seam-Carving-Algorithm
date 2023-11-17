import edu.princeton.cs.algs4.IndexMinPQ;
import edu.princeton.cs.algs4.Picture;

// finds optimal seams in images
public class SeamCarver {
    // stores picture
    private Picture pic;

    // create a seam carver object based on the given picture
    public SeamCarver(Picture picture) {
        if (picture == null)
            throw new IllegalArgumentException("null argument in constructor");
        // deep copy of the picture
        pic = new Picture(picture);
    }

    // current picture
    public Picture picture() {
        Picture returnPic = new Picture(pic);
        return returnPic;
    }

    // width of current picture
    public int width() {
        return pic.width();
    }

    // height of current picture
    public int height() {
        return pic.height();
    }

    // energy of pixel at column x and row y
    public double energy(int x, int y) {
        if (!inBounds(x, y))
            throw new IllegalArgumentException("Outside of prescribed range");

        int width = width();
        int height = height();

        int xAbove = pic.getRGB(x, (y + height - 1) % height);
        int xBelow = pic.getRGB(x, (y + 1) % height());
        int yRight = pic.getRGB((x + width - 1) % width, y);
        int yLeft = pic.getRGB((x + 1) % width, y);

        return Math.sqrt(energyHelper(xBelow, xAbove) +
                                 energyHelper(yLeft, yRight));
    }

    // helper method for finding the energy of a pixel
    private int energyHelper(int less, int more) {
        int rLess = (less >> 16) & 0xFF;
        int rMore = (more >> 16) & 0xFF;
        int redDiff = (rMore - rLess) * (rMore - rLess);

        int gLess = (less >> 8) & 0xFF;
        int gMore = (more >> 8) & 0xFF;
        int greenDiff = (gMore - gLess) * (gMore - gLess);

        int bLess = less & 0xFF;
        int bMore = more & 0xFF;
        int blueDiff = (bMore - bLess) * (bMore - bLess);

        return redDiff + greenDiff + blueDiff;
    }

    // sequence of indices for horizontal seam
    public int[] findHorizontalSeam() {
        int width = width();
        int height = height();

        // IndexMinPQ with keys that are the shortest path length
        // and values of the pixel in the picture. To assign values to the 2d
        // array, we use a formula to convert it to a 1d array
        IndexMinPQ<Double> pq = new IndexMinPQ<Double>(width * height);

        // keeps track of distances of shortest paths to each pixel
        double[][] distances = new double[width][height];
        // initializes all distances to except left column to -1
        for (int i = 1; i < width; i++) {
            for (int j = 0; j < height; j++) {
                distances[i][j] = -1.0;
            }
        }

        // enqueues all of the vertices in the first column into the pq
        // also initializes the first column of distances
        for (int i = 0; i < height; i++) {
            double energy = energy(0, i);
            pq.insert(oneDIndex(0, i), energy);
            distances[0][i] = energy;
        }

        // keeps track whether pixel has been dequeued from the array
        boolean[][] marked = new boolean[width][height];

        // stores previous y coordinate of pixel in sp to i, j
        int[][] pixelTo = new int[width][height];

        // stores the current item we have dequeued
        int curr = pq.delMin();
        marked[x(curr)][y(curr)] = true;

        // runs Djikstra's algorithm until we reach the right side of the pic
        while (x(curr) < width - 1) {
            int x = x(curr);
            int y = y(curr);

            // relax upper right
            relax(curr, oneDIndex(x + 1, y - 1), pq, marked, distances,
                  pixelTo, true);
            // relax direct right
            relax(curr, oneDIndex(x + 1, y), pq, marked, distances,
                  pixelTo, true);
            // relax lower right
            relax(curr, oneDIndex(x + 1, y + 1), pq, marked, distances,
                  pixelTo, true);

            curr = pq.delMin();
            marked[x(curr)][y(curr)] = true;
        }

        // calculates the seam
        int[] seam = new int[width];

        seam[width - 1] = y(curr);
        for (int i = width - 1; i > 0; i--) {
            seam[i - 1] = pixelTo[i][seam[i]];
        }

        return seam;
    }

    // relax helper method
    private void relax(int fromIndex, int toIndex, IndexMinPQ<Double> pq,
                       boolean[][] marked, double[][] distances,
                       int[][] pixelTo, boolean hor) {

        // corner case
        if (toIndex < 0) return;

        int xTo = x(toIndex);
        int yTo = y(toIndex);
        // checks that vertex we're relaxing is in bounds
        if (!inBounds(xTo, yTo)) return;
        // if vertex has already been marked, no need to relax it
        if (marked[xTo][yTo]) return;

        int xFrom = x(fromIndex);
        int yFrom = y(fromIndex);

        double origDist = distances[xTo][yTo];
        double relaxDist = distances[xFrom][yFrom] + energy(xTo, yTo);

        // case if sp to the toIndex hasn't even been calculated yet
        if (origDist < 0) {
            distances[xTo][yTo] = relaxDist;
            pq.insert(toIndex, relaxDist);
            if (hor) pixelTo[xTo][yTo] = yFrom;
            else pixelTo[xTo][yTo] = xFrom;
        }
        // case if sp to the toIndex already has a value
        else {
            // case if relaxed distance is less than distance stored
            if (relaxDist < origDist) {
                pq.changeKey(toIndex, relaxDist);
                distances[xTo][yTo] = relaxDist;
                if (hor) pixelTo[xTo][yTo] = yFrom;
                else pixelTo[xTo][yTo] = xFrom;
            }
        }
    }


    // helper method that checks if index is within our bounds
    private boolean inBounds(int x, int y) {
        if (x < 0 || x >= width()) return false;
        if (y < 0 || y >= height()) return false;
        return true;
    }

    // helper method
    // formula for converting a pixel into its 1-dimensional array index
    private int oneDIndex(int x, int y) {
        // corner case if x or y is out of bounds
        if (!inBounds(x, y)) return -1;
        return width() * y + x;
    }

    // returns the row of a one-dimensional index
    private int x(int index) {
        return index % width();
    }

    // returns the column of a one-dimensional index
    private int y(int index) {
        return index / width();
    }

    // sequence of indices for vertical seam
    public int[] findVerticalSeam() {
        int width = width();
        int height = height();

        // IndexMinPQ with keys that are the shortest path length
        // and values of the pixel in the picture. To assign values to the 2d
        // array, we use a formula to convert it to a 1d array
        IndexMinPQ<Double> pq = new IndexMinPQ<Double>(width * height);

        // keeps track of distances of shortest paths to each pixel
        double[][] distances = new double[width][height];
        // initializes all distances to -1 except for the top row
        for (int i = 0; i < width; i++) {
            for (int j = 1; j < height; j++) {
                distances[i][j] = -1.0;
            }
        }

        // enqueues all of the vertices in the first row into the pq
        // also initializes all the distances of the top row
        for (int i = 0; i < width; i++) {
            double energy = energy(i, 0);
            pq.insert(oneDIndex(i, 0), energy);
            distances[i][0] = energy;
        }

        // keeps track whether pixel has been dequeued from the array
        boolean[][] marked = new boolean[width][height];

        // stores previous x coordinate of pixel in sp to i, j
        int[][] pixelTo = new int[width][height];

        // stores the current item we have dequeued
        int curr = pq.delMin();
        marked[x(curr)][y(curr)] = true;

        // runs Djikstra's algorithm until we reach the bottom of the pic

        while (y(curr) < height - 1) {
            int x = x(curr);
            int y = y(curr);

            // relax lower left
            relax(curr, oneDIndex(x - 1, y + 1), pq, marked, distances,
                  pixelTo, false);
            // relax directly below
            relax(curr, oneDIndex(x, y + 1), pq, marked, distances,
                  pixelTo, false);
            // relax lower right
            relax(curr, oneDIndex(x + 1, y + 1), pq, marked, distances,
                  pixelTo, false);

            curr = pq.delMin();
            marked[x(curr)][y(curr)] = true;
        }

        // calculates the seam
        int[] seam = new int[height];

        seam[height - 1] = x(curr);
        for (int i = height - 1; i > 0; i--) {
            seam[i - 1] = pixelTo[seam[i]][i];
        }

        return seam;
    }


    // remove horizontal seam from current picture
    public void removeHorizontalSeam(int[] seam) {
        if (seam == null)
            throw new IllegalArgumentException("null argument");

        if (height() == 1)
            throw new IllegalArgumentException("Bro you're deleting the image");

        if (seam.length != width())
            throw new IllegalArgumentException("incorrect seam length");

        // corner case
        if (seam.length == 0) return;


        int height = height();
        Picture newPic = new Picture(width(), height() - 1);

        int previousSplittingPoint = seam[0];

        for (int i = 0; i < seam.length; i++) {
            int splittingPoint = seam[i];

            // checks that splittingPoint is in bounds
            if (splittingPoint < 0 || splittingPoint >= height ||
                    Math.abs(splittingPoint - previousSplittingPoint) > 1)
                throw new IllegalArgumentException("not a valid seam");

            // above to the seam
            for (int j = 0; j < splittingPoint; j++)
                newPic.setRGB(i, j, pic.getRGB(i, j));

            // below the seam
            for (int j = splittingPoint; j < height - 1; j++)
                newPic.setRGB(i, j, pic.getRGB(i, j + 1));

            previousSplittingPoint = splittingPoint;
        }

        // updates picture and energy arrays
        pic = newPic;
    }

    // remove vertical seam from current picture
    public void removeVerticalSeam(int[] seam) {
        if (seam == null)
            throw new IllegalArgumentException("null argument");

        if (width() == 1)
            throw new IllegalArgumentException("Bro you're deleting the image");

        if (seam.length != height())
            throw new IllegalArgumentException("incorrect seam length");

        if (seam.length == 0) return;

        int width = width();
        Picture newPic = new Picture(width() - 1, height());

        int previousSplittingPoint = seam[0];

        for (int i = 0; i < seam.length; i++) {
            int splittingPoint = seam[i];

            // checks that splittingPoint is in bounds
            if (splittingPoint < 0 || splittingPoint >= width ||
                    Math.abs(splittingPoint - previousSplittingPoint) > 1)
                throw new IllegalArgumentException("not a valid seam");

            // left of the seam
            for (int j = 0; j < splittingPoint; j++)
                newPic.setRGB(j, i, pic.getRGB(j, i));

            // right of the seam
            for (int j = splittingPoint; j < width - 1; j++)
                newPic.setRGB(j, i, pic.getRGB(j + 1, i));

            previousSplittingPoint = splittingPoint;
        }

        // updates picture and energy arrays
        pic = newPic;
    }

    //  unit testing (required)
    public static void main(String[] args) {
        Picture pic = new Picture(args[0]);
        int hRemove = Integer.parseInt(args[1]);
        int vRemove = Integer.parseInt(args[2]);

        SeamCarver s = new SeamCarver(pic);


        for (int i = 0; i < hRemove; i++) {
            int[] hSeam = s.findHorizontalSeam();
            s.removeHorizontalSeam(hSeam);
        }

        for (int i = 0; i < vRemove; i++) {
            int[] vSeam = s.findVerticalSeam();
            s.removeVerticalSeam(vSeam);
        }

        Picture p = s.picture();
        p.show();
    }
}
