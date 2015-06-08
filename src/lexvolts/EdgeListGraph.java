
/**
 * Created by spielman on 8/31/14.
 *
 * This might be silly.  But, I wanted a graph class that converted very quickly from matlab.
 * So, it just decorates the ijv that matlab adds to it
 * To construct it, one should type EdgeListGraph(i,j,v)
 * where [i,j,v] = find(a)
 *
 * Yes, it does want every edge twice
 */


package lexvolts;
import java.io.*;
import java.util.Scanner;

public class EdgeListGraph {

    public int[] ai;
    // Array of 'destinations' of edges
    public int[] aj;
    // Array of 'sources' of edges. This list is in sorted order
    public double[] av;
    // Lengths of edges?
    public int ne;
    public int nv;

    // nbrsStart[i] is where vertex i starts in ai
    // we don't really need aj to generate edges
    // we could find the nbrs of a vertex only using nbrsStart
    // but, we will sometimes use deg, because we will remove vertices
    public int[] nbrsStart;
    public int[] deg;

    /* aj should appear in sorted order,
       we will use ai to find the nbrs of each node
     */
    public EdgeListGraph(int[] ai, int[] aj, double[] av) {
        ne = aj.length;

        nv = aj[ne-1]+1;

        nbrsStart = new int[nv+1];
        nbrsStart[0] = 0;
        nbrsStart[nv] = ne;

        this.ai = new int[ne];
        this.aj = new int[ne];
        this.av = new double[ne];

        int curVert = 0;
        for (int i = 0; i < ne; i++) {
            this.ai[i] = ai[i];
            this.aj[i] = aj[i];
            this.av[i] = av[i];

            if (curVert != aj[i]) {
                curVert = aj[i];
                nbrsStart[curVert] = i;
            }
        }

        deg = new int[nv];

        for (int i = 0; i < nv; i++)
            deg[i] = nbrsStart[i+1] - nbrsStart[i];

    }

    public EdgeListGraph(String filename) throws java.io.IOException {

        this(new File(filename));
    }

    public EdgeListGraph(File file) throws java.io.IOException {
        Scanner sc = new Scanner(file);
        nv = sc.nextInt();
        int edgesIn = sc.nextInt();
        ne = 2*edgesIn;

        // these are what we read from the file.
        // we will sort them into appropriate order later
        int[] ai0 = new int[ne];
        int[] aj0 = new int[ne];
        double[] av0 = new double[ne];

        deg = new int[nv];

        int j = 0;
        for (int i = 0; i < edgesIn; i++) {
            int from = sc.nextInt();
            int to = sc.nextInt();
            double weight = sc.nextDouble();
            ai0[j] = from;
            aj0[j] = to;
            av0[j] = weight;
            j++;
            ai0[j] = to;
            aj0[j] = from;
            av0[j] = weight;
            j++;

            deg[from]++;
            deg[to]++;
        }

//        for (int i = 0; i < ne; i++)
//            System.out.print(ai0[i] + " " + aj0[i] + " " + av0[i] + "\n");

        nbrsStart = new int[nv+1];
        nbrsStart[0] = 0;
        int degCount = 0;
        for (int i = 0; i < nv; i++) {
            degCount = degCount + deg[i];
            nbrsStart[i + 1] = degCount;
        }

//        for (int i = 0; i < nv; i++)
//            System.out.print(deg[i] + " " + nbrsStart[i] + "\n");

        ai = new int[ne];
        aj = new int[ne];
        av = new double[ne];

        // now, let's sort edges by aj
        int[] nodePtrs = new int[nv];
        for (int i = 0; i < nv; i++)
            nodePtrs[i] = nbrsStart[i];

        for (int i = 0; i < ne; i++) {

            int to = aj0[i];
            int ptr = nodePtrs[to]++;
            ai[ptr] = ai0[i];
            aj[ptr] = to;
            av[ptr] = av0[i];
        }
    }


    // to make a copy
    public EdgeListGraph(EdgeListGraph g) {
        ne = g.ne;
        nv = g.nv;
        ai = g.ai.clone();
        aj = g.aj.clone();
        av = g.av.clone();
        nbrsStart = g.nbrsStart.clone();
        deg = g.deg.clone();
    }




    // This runs a seeded version of dijkstra: we give it the initial values
    // stored in what at where.  Note that lower values can be returned
    public double[] dijkstra(int[] where, double[] what) {
        double[] minVal = new double[nv];
        boolean[] seen = new boolean[nv];

        //Logger logger = new Logger();
        //logger.start("dijkstra.log");
        for (int i = 0; i < nv; i++) {
            minVal[i] = Double.POSITIVE_INFINITY;
            seen[i] = false;
        }

        NodeHeap pq = new NodeHeap(nv);
        for (int i = 0; i < where.length; i++) {
            pq.add(where[i], what[i]);
            //logger.write("(" + where[i] + ", " + what[i] + "), ");
        }

        while (pq.hasMore()) {
            int node = pq.popMin();
            seen[node] = true;
            double val = pq.keys[node];
            minVal[node] = val;
            //logger.write(node + " : " + val);

            for (int j = nbrsStart[node]; j < nbrsStart[node] + deg[node]; j++) {
                int nbr = ai[j];
                if (!seen[nbr]) {
                    double newVal = val + av[j];
                    pq.add(nbr, newVal);
                }
            }
        }

        return minVal;
    }


    public int[] back;

    // This runs a seeded version of dijkstra: we give it the initial values
    // stored in what at where.  Note that lower values can be returned
    // and, this one computes the tree as back-pointers, in back
    public double[] dijkstraTree(int[] where, double[] what) {
        double[] minVal = new double[nv];
        boolean[] seen = new boolean[nv];

        back = new int[nv];
        for (int i = 0; i < nv; i++)
            back[i] = -1;

        //Logger logger = new Logger();
        //logger.start("dijkstra.log");

        NodeHeap pq = new NodeHeap(nv);
        for (int i = 0; i < where.length; i++) {
            pq.add(where[i], what[i]);
            //logger.write("(" + where[i] + ", " + what[i] + "), ");
        }

        while (pq.hasMore()) {
            int node = pq.popMin();
            seen[node] = true;
            double val = pq.keys[node];
            minVal[node] = val;
            //logger.write(node + " : " + val);

            for (int j = nbrsStart[node]; j < nbrsStart[node] + deg[node]; j++) {
                int nbr = ai[j];
                if (!seen[nbr]) {
                    double newVal = val + av[j];
                    if (pq.inQueue(nbr)) {
                        if (newVal < pq.keys[nbr]) {
                            // Changed this changeKey to decreaseKey
                            pq.decreaseKey(nbr, newVal);
                            back[nbr] = node;
                        }
                    } else {
                        pq.add(nbr, newVal);
                        back[nbr] = node;
                    }

                }
            }
        }

        return minVal;
    }

}
