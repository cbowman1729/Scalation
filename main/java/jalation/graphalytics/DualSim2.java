
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
/** @author  Matthew Saltz, John Miller
 *  @version 1.1
 *  @date    Sun May 18 11:23:45 EDT 2014
 *  @see     LICENSE (MIT style license file).
 */

package jalation.graphalytics;

import java.util.Arrays;
import java.util.HashSet;
import java.util.Set;

import static java.lang.System.out;

//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
/** The `DualSim2` class provides a second implementation for Dual Graph
 *  Simulation. It differs from DualSim by not using inverse adjacency sets
 *  ('par') in order to save space.
 */
public class DualSim2
{
    /** the data graph G(V, E, l)
     */
    private final Graph g;

    /** the query graph Q(U, D, k)
     */
    private final Graph q;

    /** empty match
     */
    @SuppressWarnings("unchecked")
    private final Set <Integer> [] EMPTY = (Set <Integer> []) new HashSet <?> [0];

    /** indication of self-loops
     */
    private static final boolean SELF_LOOPS = false;

    //::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    /** Construct a Dual Graph Simulation pattern matcher.
     *
     *  @param g  the data graph G(V, E, l)
     *  @param q  the query graph Q(U, D, k)
     */
    public DualSim2 (Graph g, Graph q)
    {
        this.g = g;
        this.q = q;
    } // constructor

    //::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    /** Return the non-destructive intersection of sets 'a' and 'b'.
     *
     *  @param a  the first set
     *  @param b  the second set
     */
    public static Set <Integer> intersect (Set <Integer> a, Set <Integer> b)
    {
        Set <Integer> c = new HashSet <> (a);
        c.retainAll (b);
        return c;
    } // intersect

    //::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    /** Create an initial array of feasible mappings 'phi' from each query vertex
     *  'u' to the corresponding set of data graph vertices '{v}' whose label
     *  matches 'u's.
     */
    @SuppressWarnings("unchecked")
    public Set <Integer> [] feasibleMates ()
    {
        Set <Integer> [] phi = (Set <Integer> []) new HashSet <?> [q.size ()];
        for (int i = 0; i < q.size (); i++) {
            phi [i] = new HashSet <> (g.getVerticesWithLabel (q.label [i]));
        } // for
        return phi;
    } // feasibleMates

    //::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    /** Apply the Dual Graph Simulation pattern matching algorithm to find the
     *  mappings from the query graph 'q' to the data graph 'g'. These are
     *  represented by a multi-valued function 'phi' that maps each query graph
     *  vertex 'u' to a set of data graph vertices '{v}'.
     */
    public Set <Integer> [] mappings ()
    {
        return saltzDualSim (feasibleMates ());
    } // mappings

    //::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    /** Given the mappings 'phi' produced by the 'feasibleMates' method,
     *  eliminate mappings 'u -> v' when (1) v's children fail to match u's or
     *  (2) v's parents fail to match u's.
     *
     *  @param phi  array of mappings from a query vertex u to { graph vertices v }
     */
    public Set <Integer> [] saltzDualSim (Set <Integer> [] phi)
    {
        Set <Integer> phiTemp = null;
        Set <Integer> newPhi  = null;
        boolean alter = true;

        while (alter) {
            alter = false;
            // loop over query vertices u and u's children u_c
            for (int u = 0; u < q.size (); u++) {

                for (int u_c: q.adj [u]) {
                    newPhi = new HashSet <> ();                             // subset of phi(u_c) having a parent in phi(u)
                    for (int v: phi [u]) {                                  // data vertex v matching u's label
                        phiTemp = intersect (g.adj [v], phi [u_c]);         // children of v contained in phi(u_c)
                        if (phiTemp.isEmpty ()) {
                            phi [u].remove (v);                             // remove vertex v from phi(u)
                            if (phi [u].isEmpty ()) return EMPTY;           // no match for vertex u => no overall match
                            alter = true;
                        } // if
                        // build newPhi to contain only those vertices in phi(u_c) which also have a parent in phi(u)
                        newPhi.addAll (phiTemp);
                    } // for

                    if (newPhi.isEmpty ()) return EMPTY;                    // empty newPhi => no match 
                    if (newPhi.size () < phi [u_c].size ()) alter = true;   // since newPhi is smaller than phi(u_c)
                    if (SELF_LOOPS && u_c == u) {
                        phi [u_c] = intersect (phi [u_c], newPhi);
                    } else {
                        phi [u_c] = newPhi;
                    } // if
                } // for

            } // for
        } // while
        return phi;
    } // saltzDualSim

    //::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    /** Test the `DualSim2` class.
     *
     *  @param args  unused command-line arguments
     */
    @SuppressWarnings("unchecked")
    public static void main (String [] args)
    {
        int gSize    = 1000;      // size of the data graph
        int qSize    = 10;        // size of the query graph
        int nLabels  = 100;       // number of distinct labels
        int gAvDegree = 5;        // average vertex out degree for data graph
        int qAvDegree = 2;        // average vertex out degree for query graph

//      Graph g = GraphGenerator.genRandomGraph (gSize, nLabels, gAvDegree)
//      Graph q = GraphGenerator.genBFSQuery (qSize, qAvDegree, g)

        Set <Integer> [] adj_g = (Set <Integer> []) new Set <?> [4];
        adj_g [0] = new HashSet <> (Arrays.asList (1, 2));
        adj_g [1] = new HashSet <> (Arrays.asList (2, 3));
        adj_g [2] = new HashSet <> (Arrays.asList (3));
        adj_g [3] = new HashSet <> (Arrays.asList (1));
        Integer [] label_g = {2, 1, 2, 1};
        Graph g = new Graph (adj_g, label_g, false);

        Set <Integer> [] adj_q = (Set <Integer> []) new Set <?> [3];
        adj_q [0] = new HashSet <> (Arrays.asList (1, 2));
        adj_q [1] = new HashSet <> (Arrays.asList (2));
        adj_q [2] = new HashSet <> ();
        Integer [] label_q = {2, 1, 2};
        Graph q = new Graph (adj_q, label_q, false);

        DualSim2 matcher = new DualSim2 (g, q);                      // Dual Graph Simulation Pattern Matcher
        Set <Integer> [] phi = matcher.mappings ();                  // run the matcher
        for (int i = 0; i < phi.length; i++) out.println ("u_" + i + ": " + phi [i]);
    } // main

} // DualSim2 class

