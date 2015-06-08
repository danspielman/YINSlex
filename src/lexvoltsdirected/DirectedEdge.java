package lexvoltsdirected;

/**
 * Created by sushant on 1/31/15.
 */
public class DirectedEdge {
    // This is an edge
    // Note that the edge could be either from src to dest, **OR from dest to src**
    // The direction is determined by the boolean outgoing
    public int src; // The source
    public int dest;
    public double len;
    public boolean outgoing;
    // If outgoing = true, then the edge goes from src to dest, o/w it goes from dest to src

    public int nextEdge;
    public int prevEdge;

    public int symEdge;

    public int nextActEdge;
    public int prevActEdge;

    //    public Edge(int src, int dest, double len, double pressure) {
    public DirectedEdge(int src, int dest, double len, boolean outgoing) {
        this.src = src;
        this.dest = dest;
        this.len = len;
        this.outgoing = outgoing;
        nextEdge = -1;
        prevEdge = -1;
        symEdge = -1;
        nextActEdge = -1;
        prevActEdge = -1;
    }
}
