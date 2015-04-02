
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
/** @author  Michael Cotterell
 *  @version 1.1
 *  @date    Fri Sep 30 13:37:32 EDT 2015
 *  @see     LICENSE (MIT style license file).
 *------------------------------------------------------------------------------
 *  @see https://www.cs.cmu.edu/~ggordon/10725-F12/slides/25-coord-desc.pdf
 */

package scalation.minima

import math.{abs, max, pow}

import scalation.calculus.Calculus.FunctionV2S
import scalation.linalgebra.VectorD
import scalation.util.Error

//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
/** Coordinate Descent (CD) algorithm using a simple cyclical scheme.  CD splits
 *  the coordinates up into roughly equal sized blocks (default is one block). 
 *  It finds a minima of the objective function by performing a line search 
 *  along one coordinate direction per block per iteration.   
 * 
 *  minimize    f(x)
 *  subject to  g(x) <= 0    [ optionally g(x) == 0 ]
 *
 *  @param f     the objective function to be minimized
 *  @param g     the constraint function to be satisfied, if any
 *  @param ineq  whether the constraint function must satisfy inequality or equality
 *  @param exactLS  whether to use exact (e.g., GoldenLS)
 *                            or inexact (e.g., WolfeLS) Line Search
 *  @param nBlocks the number of blocks (default = 1) 
 */
class CoordinateDescent (f: FunctionV2S, g: FunctionV2S = null,
                         ineq: Boolean = true, exactLS: Boolean = true, nBlocks: Int = 1)
      extends Minimizer with Error
{
    private val DEBUG    = true                  // the debug flag
    private var WEIGHT   = 10.0                // weight on penalty for constraint violation
    private val BETA     = 1.5
    private var df: Array [FunctionV2S] = null   // array of partials

    //:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    /** The objective function f plus a weighted penalty based on the constraint
     *  function g.
     *  @param x  the coordinate values of the current point
     */
    override def fg (x: VectorD): Double =
    {
        val f_x = f(x)
        if (g == null) {                  // unconstrained
            f_x
        } else {                          // constrained, g(x) <= 0
            val penalty = if (ineq) max (g(x), 0.0) else abs (g(x))
//            println ("penalty = " + penalty)
            WEIGHT *= BETA
            val ret = f_x + WEIGHT * penalty                           // abs (f_x) * WEIGHT * penalty * penalty
//            println ("fg(x) = " + ret)
            ret
        } // if
    } // fg

    //::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    /** Perform an exact (GoldenSectionLS) or inexact (WolfeLS) line search.
     *  Search along the coordinate, returning the distance 'z' to move that 
     *  coordinate.
     *  @param x     the current point
     *  @param i     the coordinate in which to do the line search
     *  @param step  the initial step size
     */
    def lineSearchC (x: VectorD, dir: VectorD, step: Double, i: Int): Double =
    {
        val y = new VectorD(x)
        def f_1D (z: Double): Double = {                 // create a 1D function
            y(i) = y(i) + dir(i) * z; fg(y)              // only change the current coordinate
        } // f_1D
        val ls = if (exactLS) new GoldenSectionLS (f_1D) // Golden Section line search
                 else new WolfeLS (f_1D, .0001, .1)      // Wolfe line search (c1 = .0001, c2 = .1)
        ls.search (STEP)                                  // perform a line search
    } // lineSearchC

    //::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    @deprecated ("For `CoordinateDescent` use `lineSearchC`", "1.1.1")
    def lineSearch (x: VectorD, dir: VectorD, step: Double = STEP): Double = 0.0

    //::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    /** Returns an iterator with the coordinate blocks based on `nBlocks`.
     *  @param x  an input vector
     */
    def getBlocks (x: VectorD) = (0 until x.dim).grouped (x.dim / nBlocks).toList

    //::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    /** Solve the Non-Linear Programming (NLP) problem using the CD algorithm.
     *  Currently, a cyclical approach is taken for choosing the next 
     *   coordinate.
     *  @param x0     the starting point 
     *  @param step   the initial step size
     *  @param toler  the tolerence
     */
    def solve (x0: VectorD, step: Double = STEP, toler: Double = TOL): VectorD =
    {

        var x      = x0                                     // current point
        var f_x    = fg(x)                                  // objective function at current point
        var y:   VectorD = null                             // next point
        var y_P: VectorD = null                             // potential next point (positive line search)
        var y_N: VectorD = null                             // potential next point (negative line search)
        var f_y_P  = 0.0                                    // objective function after positive line search
        var f_y_N  = 0.0                                    // objective function after negative line search
        var f_y    = 0.0                                    // objective function at next point
        val dir_P  = new VectorD (x.dim); dir_P.set ( 1)    // positive direction
        val dir_N  = new VectorD (x.dim); dir_N.set (-1)    // negative direction
        var dist   = 1.0                                    // distance between current and next point
        var down   = true                                   // moving down flag
        var check  = true                                   // did we change a coordinate
        val blocks = getBlocks (x0)                         // coordinate blocks

        for (k <- 1 to MAX_ITER if dist > toler && check) {

            check  = false                                  // reset check flag

            for (block <- blocks; i <- block) {

                if (DEBUG) println ("solve: k = " + k + ", block = " + block + ", coordinate i = " + i)

                y_P    = new VectorD (x)                    // positive line search
                y_P(i) = y_P(i) + lineSearchC (x, dir_P, step, i)
                f_y_P  = fg(y_P)

                if (DEBUG) println ("solve: k = " + k + ", i = " + i + ", y_P = " + y_P + ", f_y_P = " + f_y_P)

                y_N    = new VectorD (x)                    // negative line search
                y_N(i) = y_N(i) + lineSearchC (x, dir_N, step, i)
                f_y_N  = fg(y_N)

                if (DEBUG) println ("solve: k = " + k + ", i = " + i + ", y_N = " + y_N + ", f_y_N = " + f_y_N)

                if (f_y_P < f_y_N) {                       // determine obj function value for next point
                    y   = y_P
                    f_y = f_y_P
                } else {
                    y   = y_N
                    f_y = f_y_N
                } // if

                if (DEBUG) println ("solve: k = " + k + ", y = " + y + ", f_y = " + f_y + ", i = " + i)

                dist = (x - y).normSq                      // calc the distance between current and next point
                down = f_y < f_x                           // still moving down?

                if (down) {                                // make the next point, the current point
                    x     = y
                    f_x   = f_y
                    check = true
                } // if

                if (DEBUG) println ("solve: k = " + k + ", x = " + x + ", f_x = " + f_x + ", i = " + i)

            } // for 

            if (DEBUG) println ("solve: k = " + k + ", x = " + x + ", f_x = " + f_x)

        } // for
        x                                                  // return the current point
    } // solve

} // CoordinateDescent class


//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
/** This object is used to test the CoordinateDescent class.
 */
object CoordinateDescentTest extends App
{
    val x0 = new VectorD (2)

    println ("\nMinimize: (x_0 - 3)^2 + (x_1 - 4)^2 + 1")
    def f (x: VectorD): Double = (x(0) - 3.0) * (x(0) - 3.0) + (x(1) - 4.0) * (x(1) - 4.0) + 1.0
    var solver = new CoordinateDescent (f)
    var x = solver.solve (x0)
    println ("optimal solution x = " + x + " with an objective value f(x) = " + f(x))

    println ("\nMinimize: x_0^4 + (x_0 - 3)^2 + (x_1 - 4)^2 + 1")
    def g (x: VectorD): Double = pow (x(0), 4.0) + (x(0) - 3.0) * (x(0) - 3.0) + (x(1) - 4.0) * (x(1) - 4.0) + 1.0
    solver = new CoordinateDescent (g)
    x = solver.solve (x0)
    println ("optimal solution x = " + x + " with an objective value g(x) = " + g(x))

} // CoordinateDescentTest object

