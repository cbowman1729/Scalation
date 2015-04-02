
//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
/** @author  John Miller, Michael Cotterell
 *  @version 1.1
 *  @date    Sun May  5 13:13:42 EDT 2013
 *  @see     LICENSE (MIT style license file).
 *
 *  @see     http://www.scholarpedia.org/article/Nelder-Mead_algorithm
 */

package scalation.minima

import math.{abs, max, pow}
import util.control.Breaks.{breakable, break}
import util.Sorting.stableSort

import scalation.calculus.Calculus.FunctionV2S
import scalation.linalgebra.VectorD
import scalation.util.Error

//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
/** This class solves Non-Linear Programming (NLP) problems using the 
 *  Nelder-Mead Simplex algorithm.  Given a function 'f' and its dimension
 *  'n',  the algorithm moves a simplex defined by n + 1 points in order to find
 *  an optimal solution.  The algorithm is derivative-free. 
 *
 *  minimize    f(x)
 *
 *  @param f  the vector-to-scalar objective function
 *  @param n  the dimension of the search space
 */
class NelderMeadSimplex (f: FunctionV2S, n: Int)
      extends Minimizer with Error
{
    type Vertex = Tuple2 [VectorD, Double]            // point and its functional value

    private val DEBUG   = true                        // debug flag
    private val WEIGHT  = 1000.0                      // weight on penalty for constraint violation
    private val np1     = n + 1                       // number of vertices/points in simplex
    private val simplex = Array.ofDim [Vertex] (np1)  // simplex used for search

    private val alpha = 1.0                           // alpha (> 0)  parameter for reflection
    private val beta  = 0.5                           // beta  (0, 1) parameter for constraction
    private val gamma = 2.0                           // gamma (> 1)  parameter for expansion
    private val delta = 0.5                           // delta (0, 1) parameter for shrinkage

    private def x_h = simplex(0)._1                   // worst point
    private def x_s = simplex(1)._1                   // second worst point
    private def x_l = simplex(n)._1                   // best point

    private def f_h = simplex(0)._2                   // worst point functional value
    private def f_s = simplex(1)._2                   // second worst point functional value
    private def f_l = simplex(n)._2                   // best point functional value

    private var g: FunctionV2S = null                 // constraint function
    var evals                  = 0                    // number of objective function evaluations

    if (n < 2) flaw ("constructor", "requires at least a 2-dimensional problem")

    //::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    /** Constrined Nelder-Mead Simplex Algorith.
     * 
     *  minimize    f(x)
     *  subject to  g(x) <= 0    [ optionally g(x) == 0 ]
     * 
     *  @param f  the vector-to-scalar objective function
     *  @param g  the constraint function to be satisfied, if any
     *  @param n  the dimension of the search space
     */
    def this (f: FunctionV2S, g: FunctionV2S, n: Int)
    {
        this (f, n)
        this.g = g
    } // this

    //:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    /** The objective function f plus a weighted penalty based on the constraint
     *  function g.
     *  @param x  the coordinate values of the current point
     */
    override def fg (x: VectorD): Double =
    {
        val f_x = f(x); evals += 1
        if (g == null) {                  // unconstrained
            f_x
        } else {                          // constrained, g(x) <= 0
            val penalty = max (g(x), 0.0)
            f_x + abs (f_x) * WEIGHT * penalty * penalty
        } // if
    } // fg

    //::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    /** Initialize the search simplex by setting n + 1 vertices and computing
     *  their functional values.
     *  @param x0  the given starting point
     *  @param step  the step size
     */
    def initSimplex (x0: VectorD, step: Double)
    {
        simplex(0) = (x0, fg(x0))                 // given starting point and its functional value
        for (i <- 1 to n) {
            val x = x0 + x0.oneAt (i-1) * step
            simplex(i) = (x, fg(x))
        } // for
        sort ()                                   // order vertices high to low
    } // initSimplex

    //::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    /** Sort the vertices in non-increasing order (high to low).  Then the key
     *  indices are worst/highest (h=0), second worst (s=1), and best/lowest (l=n).
     */
    def sort () = stableSort (simplex, (v1: Vertex, v2: Vertex) => v1._2 > v2._2)

    //::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    /** Calculate the centroid of the best-side of the simplex (excluding h=0),
     *  returning it and its functional value.
     */
    def centroid (): Vertex =
    {
        val c = new VectorD (n)                   // the centroid of the simplex
        for (i <- 1 to n) c += simplex(i)._1      // add vertex points, except h=0
        val x_c = c / n.toDouble                  // divide by # vertices - 1
        (x_c, fg(x_c))
    } // centroid

    //::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    /** Reflect: compute the reflection point of the worst point (h=0) across
     *  the centroid.
     *  @param x_c  the best-side centroid of the simplex
     */
    def reflect (x_c: VectorD): Vertex =
    {
        val x_r = x_c + (x_c - x_h) * alpha
        (x_r, fg(x_r))
    } // reflect

    //::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    /** Expand:  compute the expansion point beyond the reflection point.
     *  @param x_c  the best-side centroid of the simplex
     *  @param x_r  the reflection point
     */
    def expand (x_c: VectorD, x_r: VectorD): Vertex = 
    {
        val x_e = x_c + (x_r - x_c) * gamma
        (x_e, fg(x_e))
    } // expand

    //::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    /** Contract: compute the outer contraction point between x_r and x_c.
     *  @param x_c  the best-side centroid of the simplex
     *  @param x_r  the reflection point
     */
    def contractOut (x_c: VectorD, x_r: VectorD): Vertex =
    {
        val x_co = x_c + (x_r - x_c) * beta
        val out  = (x_co, fg(x_co))
        if (DEBUG) println("contractOut: out = " + out)
        out
    } // contractOut

    //::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    /** Contract: compute the inner contraction point between x_h and x_c.
     *  @param x_c  the best-side centroid of the simplex
     *  @param x_r  the reflection point
     */
    def contractIn (x_c: VectorD, x_r: VectorD): Vertex =
    {
        val x_ci = x_c + (x_h - x_c) * beta
        val in   = (x_ci, fg(x_ci))
        if (DEBUG) println("contractOut: in = " + in)
        in
    } // contractIn

    //::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    /** Shrink: fixing the best/lowest point (l=n), move the rest of the points
     *  toward it.
     */
    def shrink ()
    {
        if (DEBUG) println ("shrink: before = " + simplex.deep)
        for (i <- 0 until n) {
            val x = x_l + (simplex(i)._1 - x_l) * delta
            simplex(i) = (x, fg(x))                          // updated vertex
        } // for
        sort ()
        if (DEBUG) println ("shrink: after  = " + simplex.deep)
    } // shrink

    //::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    /** Perform an exact (e.g., GoldenSectionLS) or inexact (e.g., WolfeLS) line search.
     *  Search in direction 'dir', returning the distance 'z' to move in that direction.
     *  Currently NOT USED, but may be used to find a better point to add to simplex.
     *  @param x     the current point
     *  @param dir   the direction to move in
     *  @param step  the initial step size
     */
    def lineSearch (x: VectorD, dir: VectorD, step: Double = STEP): Double = 0.0

    //::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    /** Replace the worst vertex (h=0) in the simplex with the new point.
     *  @param  x_n  the new replacement point
     */
    def replace (x_n: VectorD)
    {
        simplex(0) = (x_n, fg(x_n))                                   // put new vertex at index h=0
        if (DEBUG) println ("replace: " + 0 + " with " + simplex(0))
        sort ()                                                       // re-establish the vertex order
    } // replace

    //::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    /** Improve the simplex by replacing the worst/highest vertex (x_h) with a
     *  a better one found on the line containing x_h and the centroid (x_c).
     *  Try the reflection, expansion, outer contraction and inner contraction
     *  points, in that order.  If none succeeds, shrink the simplex and iterate.
     *  Return both distance and difference between x_h (worst) and x_l (best).
     *  @param toler  the tolerence used for termination
     */
    def improveSimplex (toler: Double = EPSILON): Tuple2 [Double, Double] =
    {
        var dist = (simplex(0)._1 - simplex(n)._1).norm          // distance between x_h and x_l
        var diff =  simplex(0)._2 - simplex(n)._2                // difference between f_h and f_l

        breakable { for (k <- 1 to MAX_ITER) {

            val (x_c, f_c) = centroid ()                         // compute best-side centroid of simplex
            val (x_r, f_r) = reflect (x_c)                       // compute reflection point

            val smaller = f_r <  f_l                             // f_r smaller than best
            val larger  = f_r >= f_s                             // f_r at least as large as second worst

            if (smaller) {
                val (x_e, f_e) = expand (x_c, x_r)               // expand beyond reflection point
                if (f_e < f_r) { replace (x_e); break }          // replace worst x_h with x_e
                else           { replace (x_r); break }          // replace worst x_h with x_r
            } // if

            if (larger) {                                        // contract back from reflection point
                if (f_r < f_h) {                                 // f_r between second worst and worst
                    val (x_co, f_co) = contractOut (x_c, x_r)
                    if (f_co <= f_r) { replace (x_co); break }   // replace worst x_h with x_co
                } else {                                         // f_r at least as large as worst
                    val (x_ci, f_ci) = contractIn (x_c, x_r)
                    if (f_ci <= f_h) { replace (x_ci); break }   // replace worst x_h with x_ci
                } // if
            } // if

            shrink ()

            dist = (x_h - x_l).norm          // recompute the distance
            diff = f_h - f_l                 // recompute the difference

            if (dist < toler || diff < toler) break              // check termination condition

        }} // for

        (dist, diff)                                             // return distance and difference
    } // improveSimplex

    //::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    /** Solve the Non-Linear Programming (NLP) problem using the Nelder-Mead
     *  Simplex algorithm.
     *  @param x0     the given starting point
     *  @param step   the initial step size
     *  @param toler  the tolerence used for termination
     */
    def solve (x0: VectorD, step: Double = STEP, toler: Double = EPSILON): VectorD =
    {
        evals = 0
        var dist = Double.PositiveInfinity                     // distance between worst and best vertices
        var diff = Double.PositiveInfinity                     // difference between their functional values
        initSimplex (x0, step)
        if (DEBUG) println ("solve: " + 0 + ":\tdist = " + dist + " diff = " + diff + "\n\tsimplex = " + simplex.deep)

        breakable { for (k <- 1 to MAX_ITER) {
            val (dist, diff) = improveSimplex ()
            if (DEBUG) println ("solve: " + k + ":\tdist = " + dist + " diff = " + diff + "\n\tsimplex = " + simplex.deep)
            if (dist < toler || diff < toler) break            // check termination condition
        }} // for

        println ("solve: optimal vertex = " + x_l)
        x_l                                                    // return the best vertex point
    } // solve

} // NelderMeadSimplex class


//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
/** This object is used to test the NelderMeadSimplex class.
 */
object NelderMeadSimplexTest extends App
{
    var x0 = VectorD (1.0, 1.0)                           // starting point

    println ("\nProblem 1: (x_0 - 2)^2 + (x_1 - 3)^2 + 1") 
    def f (x: VectorD): Double = (x(0) - 2.0) * (x(0) - 2.0) + (x(1) - 3.0) * (x(1) - 3.0) + 1.0
    val solver = new NelderMeadSimplex (f, 2)
    val x = solver.solve (x0)                             // optimal point
    println ("optimal solution = " + x + " at " + f(x))

} // NelderMeadSimplexTest

//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
/** This object is used to test the NelderMeadSimplex class.
 */
object NelderMeadSimplexTest2 extends App
{
    val x0 = new VectorD (2)

    println (":::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::")
    println
    println ("Minimize: (x_0 - 3)^2 + (x_1 - 4)^2 + 1")

    def f (x: VectorD): Double = (x(0) - 3.0) * (x(0) - 3.0) + (x(1) - 4.0) * (x(1) - 4.0) + 1.0
    var solver = new NelderMeadSimplex (f, x0.dim)
    var x = solver.solve (x0)

    println
    println ("Minimize: (x_0 - 3)^2 + (x_1 - 4)^2 + 1")
    println ("optimal solution x = " + x + " with an objective value f(x) = " + f(x))
    println ("evals = " + solver.evals)
    println

    println (":::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::")
    println
    println ("\nMinimize: x_0^4 + (x_0 - 3)^2 + (x_1 - 4)^2 + 1")

    def g (x: VectorD): Double = pow (x(0), 4.0) + (x(0) - 3.0) * (x(0) - 3.0) + (x(1) - 4.0) * (x(1) - 4.0) + 1.0
    solver = new NelderMeadSimplex (g, x0.dim)
    x = solver.solve (x0)

    println
    println ("Minimize: x_0^4 + (x_0 - 3)^2 + (x_1 - 4)^2 + 1")
    println ("optimal solution x = " + x + " with an objective value g(x) = " + g(x))
    println ("evals = " + solver.evals)

} // NelderMeadSimplexTest2

