
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
/** @author  Michael Cotterell
 *  @version 1.1
 *  @date    Sun Mar 28 16:32:04 EDT 2015
 *  @see     LICENSE (MIT style license file).
 */

package scalation.util

import scalation.calculus.Calculus.FunctionV2S
import scalation.linalgebra.VectorD

import collection.mutable.HashMap

//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
/** Cached version of FunctionV2S. Supports averaging new evaluations of a point
 *  with all previous evaluations of the same point.
 *  @param f    the function
 *  @param avg  whether to enable averaging (default = false)
 */
class CachedFunctionV2S (f: FunctionV2S, avg: Boolean = false) extends FunctionV2S
{

    private val cache = HashMap.empty [VectorD, (Double, Int)]

    override def apply (x: VectorD): Double =
    {
        var (value, n) = cache.getOrElseUpdate (x, (f(x), 1))
        if (avg && n > 0) {                          
            value = left (value, n + 1) + right (f(x), n + 1)
            cache.put (x, (value, n + 1))
        } // avg
        value
    } // apply

    private def left (fx: Double, n: Int)  = (fx / n) * (n - 1)
    private def right (fx: Double, n: Int) = (fx / n) 

    //::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    /** Clears the cache.
     */
    def clear () = cache.clear ()

} // CachedFunctionV2S

//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
/** Tests for the `CachedFunctionV2S` class.
 */
object CachedFunctionV2STest extends App {

    def f (x: VectorD) = {
        println ("evaluating")
        x(0) + 2
    } // f

    val cf = new CachedFunctionV2S (f)

    println("Test 1: Should only print \"evaluating\" once...")

    cf(new VectorD(1))
    cf(new VectorD(1))

    println("Test 2: Clearing the cache; should only print \"evaluating\" once...")

    cf.clear ()
    cf(new VectorD(1))
    cf(new VectorD(1))

    println("Test 3: Testing averaging...")

    import scalation.random.Normal
    val nv = Normal ()

    def sf (x: VectorD) = x(0) + nv.gen
    val csf = new CachedFunctionV2S (sf, true)

    for (i <- 0 until 100000) csf(new VectorD (1))
    println("eval = " + csf(new VectorD (1)))

} // CachedFunctionV2STest
