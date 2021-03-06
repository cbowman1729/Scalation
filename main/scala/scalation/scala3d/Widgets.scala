
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
/** @author  John Miller
 *  @version 1.1
 *  @date    Tue May 13 16:18:42 EDT 2014
 *  @see     LICENSE (MIT style license file).
 *
 *  Intended to make switching GUI's easier.
 */

package scalation.scala3d

import javafx.application.Application
import javafx.event.{Event, EventHandler}
import javafx.scene.Group
import javafx.scene.input.MouseEvent
import javafx.scene.paint.{Color, PhongMaterial}
import javafx.scene.transform.{Rotate, Scale, Translate}

//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
/** `AppHelper` is a convenience class for `Application`.
 */
case class AppHelper() extends javafx.application.Application()
{
    val root = new Group()
    val camera = new javafx.scene.PerspectiveCamera(true)
    var anchorX: Double = _
    var anchorY: Double = _
    var anchorAngle: Double = _

    def mkEventHandler[E <: Event](f: E => Unit) = new EventHandler[E] { def handle(e: E) = f(e) }

    override def start (primaryStage: javafx.stage.Stage)
    {
        println("start")

        buildCamera()
        //Create Box (aka teller)
        val c = new Box(20,20,20)
        val boxMaterial = new PhongMaterial()
        boxMaterial.setDiffuseColor(Color.GREEN)
        boxMaterial.setSpecularColor(Color.WHITESMOKE)
        c.setMaterial(boxMaterial)
        val translateBox = new Translate(20,0,0)
        c.getTransforms().add(translateBox)

        //Create Sphere (aka source/sink)
        val sphere = new Sphere(10.0)
        val sphereMaterial = new PhongMaterial()
        sphereMaterial.setDiffuseColor(Color.GREEN)
        sphereMaterial.setSpecularColor(Color.WHITESMOKE)
        sphere.setMaterial(sphereMaterial)
        val translateSphere = new Translate(-100,0,0)
        sphere.getTransforms().add(translateSphere)

        //Create Sphere (aka source/sink)
        val sphere2 = new Sphere(10.0)
        val sphereMaterial2 = new PhongMaterial()
        sphereMaterial2.setDiffuseColor(Color.RED)
        sphereMaterial2.setSpecularColor(Color.WHITESMOKE)
        sphere2.setMaterial(sphereMaterial2)
        val translateSphere2 = new Translate(100,0,0)
        sphere2.getTransforms().add(translateSphere2)

        //Create Cylinder (aka tellerQ)
        val cylinder = new Cylinder(10, 20)
        val cylinderMaterial = new PhongMaterial()
        cylinderMaterial.setDiffuseColor(Color.ORANGE)
        cylinderMaterial.setSpecularColor(Color.WHITESMOKE)
        cylinder.setMaterial(cylinderMaterial)
        cylinder.getTransforms().add(new Rotate(90,0,0))

        val group1 = new Group(sphere)
        val group2 = new Group(sphere2)
        val group3 = new Group(cylinder)
        val group4 = new Group(c)
        root.getChildren().addAll(group1, group2, group3, group4)

        try {
            val scene = new Scene(root, 1024, 768, true)

            scene.setOnMousePressed (mkEventHandler ((event: MouseEvent) =>
            {
                anchorX = event.getSceneX ()
                anchorY = event.getSceneY ()
                anchorAngle = camera.getRotate ()
            }))

            scene.setOnMouseDragged(mkEventHandler((event: MouseEvent) =>
            {
                println ("Rotate: " + (anchorAngle + anchorX - event.getSceneX ()))
                camera.setRotate (anchorAngle + anchorX - event.getSceneX ())
            }))

            scene.setFill (javafx.scene.paint.Color.GREY)
            primaryStage.setTitle ("Scala3d sample")
            primaryStage.setScene (scene)
            primaryStage.show ()
            scene.setCamera (camera)
        } catch {
            case e: Exception => e.printStackTrace ()
        } // try
    } // start
    
    def buildCamera()
    {
        val translateZ = new Translate (0,0,-1000)
        camera.getTransforms ().add (translateZ)
        camera.setNearClip (0.0001)
        camera.setFarClip (1000.0)
    } // buildCamera

} // AppHelper


//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
/** `Stage` is a convenience class for `Stage`.
 */
class Stage extends javafx.stage.Stage

//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
/** `BorderPane` is a convenience class for `BorderPane`.
 */
class BorderPane extends javafx.scene.layout.BorderPane


//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
/** `Scene` is a convenience class for `Scene`.
 */
class Scene(root: javafx.scene.Parent, width: Int, height: Int, depthField: Boolean)
      extends javafx.scene.Scene (root, width, height, depthField)


//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
/** The `AppTest` is used to test the `AppHelper` class.
 */
object AppTest
{
    def main(args: Array [String]) {
        javafx.application.Application.launch (classOf [AppHelper], args: _*)
    } // main

} // AppTest

