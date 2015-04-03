
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
/** @author  John Miller
 *  @version 1.1
 *  @date    Mon Sep  7 15:05:06 EDT 2009
 *  @see     LICENSE (MIT style license file).
 */

package scalation.process

import java.util.concurrent.ConcurrentLinkedQueue 

import collection.mutable.ListBuffer
import util.control.Breaks.{break, breakable}

import scalation.animation.CommandType._
import scalation.random.Variate
import scalation.scala2d.Ellipse
import scalation.scala2d.Colors._
import scalation.util.Monitor.trace

//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
/** The `Source` class is used to periodically inject entities (`SimActors`) into a
 *  running simulation model.  May act as an arrival generator.  Source is both
 *  a simulation `Component` and special `SimActor` and therefore runs in own thread.
 *  @param name          the name of the source
 *  @param director      the director controlling the model
 *  @param makeEntity    the function to make entities of a specified type
 *  @param subtype       indicator of the subtype of the entities to me made
 *  @param units         the number of entities to make
 *  @param iArrivalTime  the inter-arrival time distribution
 *  @param at            the location of the source (x, y, w, h)
 */
class Source (name: String, director: Model, makeEntity: () => SimActor, subtype: Int, units: Int,
              iArrivalTime: Variate, at: Array [Double])
      extends SimActor (name, director) with Component
{
    initStats (name)
    setAt (at)

    /** List holding SimActors that have left via a Sink.
     */
    private val zombies = new ConcurrentLinkedQueue [SimActor] ()

    //::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    /** Auxiliary constructor that uses defaults for width (w) and height (h).
     *  @param name          the name of the source
     *  @param director      the director controlling the model
     *  @param makeEntity    the function to make entities of a specified type
     *  @param units         the number of entities to make
     *  @param iArrivalTime  the inter-arrival time distribution
     *  @param xy            the (x, y) coordinates for the top-left corner of the source.
     */
    def this (name: String, director: Model, makeEntity: () => SimActor, subtype: Int, units: Int,
              iArrivalTime: Variate, xy: Tuple2 [Double, Double])
    {
        this (name, director, makeEntity, subtype, units, iArrivalTime, Array (xy._1, xy._2, 20.0, 20.0))
    } // constructor

    //::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    /** Add a zombie actor to add to queue of zombies.
     *  @param z  the zombie actor to add to queue
     */
    def addZombie (z: SimActor) { z.zombify; zombies.add (z) }

    //::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    /** Display this source as a node on the animation canvas.
     */
    def display ()
    {
        director.animate (this, CreateNode, limegreen, Ellipse (), at)
    } // display

    //::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    /** The Sources as a special SimActor will act over time to make entities
     *  (other SimActors).
     */
    def act ()
    {
        for (rep <- 1 to director.reps) {                                // major loop - replications
            actTime = director.clock                                     // set to model start time

            breakable { for (i <- 1 to units) {                                       // minor loop - make actors
                if (director.stopped) break
                val (actor, isZombie) = if (zombies.isEmpty) (makeEntity (), false)   // make new actor
                                        else (zombies.poll (), true)                  // reuse zombie actor

                if (isZombie) actor.setArrivalT (director.clock)         // update zombie's re-creation time
                actor.setMySource (this)                                 // actor's source
                actor.setSubtype (subtype)                               // set the subtype
                trace (this, "generates", actor, director.clock)
                director.animate (actor, CreateToken, randomColor (actor.id), Ellipse (),
                         Array (at(0) + at(2) + RAD / 2.0, at(1) + at(3) / 2.0 - RAD))
                actor.schedule (0.0)

                if (i < units) {
                    val duration = iArrivalTime.gen
                    tally (duration)
                    schedule (duration)
                    yieldToDirector ()                       // yield and wait duration time units
                } // if
            }} // for
            if (rep == director.reps) zombify ()
            yieldToDirector ()                              // yield and wait for next replication

        } // for
        // model must terminate this source actor
    } // act

    //::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    /** Clean-up all of the zombie actors for this source.  The model should
     *  call this method after all replications have been run.
     */
    def cleanup ()
    {
        while (! zombies.isEmpty) {
            val z = zombies.poll ()
            z.yyield (null, true)                        // yield and terminate zombie actor
            z.interrupt ()                               // must interrupt zombie actor to terminate
        } // while
    } // cleanup

} // Source class


//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
/** The `Source` companion object provides a builder method for sources.
 */
object Source
{
    //::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    /** Create a source using defaults for width (w) and heigth (h).
     *  @param name          the name of the source
     *  @param director      the director controlling the model
     *  @param makeEntity    the function to make entities of a specified type
     *  @param units         the number of entities to make
     *  @param iArrivalTime  the inter-arrival time distribution
     *  @param xy            the (x, y) coordinates for the top-left corner of the source.
     */
    def apply (name: String, director: Model, makeEntity: () => SimActor, subtype: Int, units: Int,
              iArrivalTime: Variate, xy: Tuple2 [Int, Int]): Source =
    {
        new Source (name, director, makeEntity, subtype, units, iArrivalTime,
                    Array (xy._1.toDouble, xy._2.toDouble, 20.0, 20.0))
    } // apply

    //::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    /** Create a group of related sources using defaults for width (w) and heigth (h).
     *  @param director      the director controlling the model
     *  @param makeEntity    the function to make entities of a specified type
     *  @param units         the number of entities to make
     *  @param xy            the (x, y) coordinates for the top-left corner of the reference source.
     *  @param src           repeated source specific info: name, subtype, distribution, offset
     
    def group (director: Model, makeEntity: () => SimActor, units: Int, xy: Tuple2 [Int, Int],
               src: Tuple4 [String, Int, Tuple2 [Variate, Variate], Tuple2 [Int, Int]]*): List [Source] =
    {
        val sourceGroup = new ListBuffer [Source] ()
        var main = true
        for (s <- src) { 
            if (main) {
                sourceGroup += Source (s._1, director, makeEntity, s._2, units, s._3._1,
                                      (xy._1 + s._4._1, xy._2 + s._4._2))
            } else {
                sourceGroup += Source (s._1, director, makeEntity, s._2, units, s._3._2,
                                      (xy._1 + s._4._1, xy._2 + s._4._2))
            } // if
            main = ! main
        } // for
            
        sourceGroup.toList
    } // group
*/
    //::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    /** Create a group of related sources using defaults for width (w) and heigth (h).
     *  @param director      the director controlling the model
     *  @param makeEntity    the function to make entities of a specified type
     *  @param units         the number of entities to make
     *  @param xy            the (x, y) coordinates for the top-left corner of the reference source.
     *  @param src           repeated source specific info: name, subtype, distribution, offset
     */
    def group (director: Model, makeEntity: () => SimActor, units: Int, xy: Tuple2 [Int, Int],
               src: Tuple4 [String, Int, Variate, Tuple2 [Int, Int]]*): List [Source] =
    {
        val sourceGroup = new ListBuffer [Source] ()
        for (s <- src) sourceGroup += Source (s._1, director, makeEntity, s._2, units, s._3,
                                             (xy._1 + s._4._1, xy._2 + s._4._2))
        sourceGroup.toList
    } // group    

} // Source object

