
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
/** @author  John Miller, Jason Kaine
 *  @version 1.1
 *  @date    Sat Dec 20 13:06:20 EST 2014
 *  @see     LICENSE (MIT style license file).
 */

package scalation.scala3d

import javafx.scene.paint.Color

import util.Random
 
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
/** The `Colors` convenience object defines numerous common colors.
 *  Source of colors:  http://www.webmoments.com/colorchart.htm
 */
object Colors
{
    type Color = javafx.scene.paint.Color

    val antiquewhite     = Color.ANTIQUEWHITE
    val aqua             = Color.AQUA
    val aquamarine       = Color.AQUAMARINE
    val azure            = Color.AZURE
    val beige            = Color.BEIGE
    val bisque           = Color.BISQUE
    val black            = Color.BLACK
    val blanchedalmond   = Color.BLANCHEDALMOND
    val blue             = Color.BLUE
    val blueviolet       = Color.BLUEVIOLET
    val brown            = Color.BROWN
    val burlywood        = Color.BURLYWOOD
    val cadetblue        = Color.CADETBLUE
    val chartreuse       = Color.CHARTREUSE
    val chocolate        = Color.CHOCOLATE
    val coral            = Color.CORAL
    val cornflowerblue   = Color.CORNFLOWERBLUE
    val cornsilk         = Color.CORNSILK
    val crimson          = Color.CRIMSON
    val cyan             = Color.CYAN
    val darkblue         = Color.DARKBLUE
    val darkcyan         = Color.DARKCYAN
    val darkgoldenrod    = Color.DARKGOLDENROD
    val darkgray         = Color.DARKGRAY
    val darkgreen        = Color.DARKGREEN
    val darkkhaki        = Color.DARKKHAKI
    val darkmagenta      = Color.DARKMAGENTA
    val darkolivegreen   = Color.DARKOLIVEGREEN
    val darkorange       = Color.DARKORANGE
    val darkorchid       = Color.DARKORCHID
    val darkred          = Color.DARKRED
    val darksalmon       = Color.DARKSALMON
    val darkseagreen     = Color.DARKSEAGREEN
    val darkslateblue    = Color.DARKSLATEBLUE
    val darkslategray    = Color.DARKSLATEGRAY
    val darkturquoise    = Color.DARKTURQUOISE
    val darkviolet       = Color.DARKVIOLET
    val deeppink         = Color.DEEPPINK
    val deepskyblue      = Color.DEEPSKYBLUE
    val dimgray          = Color.DIMGRAY
    val dodgerblue       = Color.DODGERBLUE
    val firebrick        = Color.FIREBRICK
    val floralwhite      = Color.FLORALWHITE
    val forestgreen      = Color.FORESTGREEN
    val fuchsia          = Color.FUCHSIA
    val gainsboro        = Color.GAINSBORO
    val ghostwhite       = Color.GHOSTWHITE
    val gold             = Color.GOLD
    val goldenrod        = Color.GOLDENROD
    val gray             = Color.GRAY
    val green            = Color.GREEN
    val greenyellow      = Color.GREENYELLOW
    val honeydew         = Color.HONEYDEW
    val hotpink          = Color.HOTPINK
    val indianred        = Color.INDIANRED
    val indigo           = Color.INDIGO
    val ivory            = Color.IVORY
    val khaki            = Color.KHAKI
    val lavender         = Color.LAVENDER
    val lavenderblush    = Color.LAVENDERBLUSH
    val lawngreen        = Color.LAWNGREEN
    val lemonchiffon     = Color.LEMONCHIFFON
    val lightblue        = Color.LIGHTBLUE
    val lightcoral       = Color.LIGHTCORAL
    val lightcyan        = Color.LIGHTCYAN
    val lightgoldenrodyellow = Color.LIGHTGOLDENRODYELLOW
    val lightgreen       = Color.LIGHTGREEN
    val lightgrey        = Color.LIGHTGREY
    val lightpink        = Color.LIGHTPINK
    val lightsalmon      = Color.LIGHTSALMON
    val lightseagreen    = Color.LIGHTSEAGREEN
    val lightskyblue     = Color.LIGHTSKYBLUE
    val lightslategray   = Color.LIGHTSLATEGRAY
    val lightsteelblue   = Color.LIGHTSTEELBLUE
    val lightyellow      = Color.LIGHTYELLOW
    val lime             = Color.LIME
    val limegreen        = Color.LIMEGREEN
    val linen            = Color.LINEN
    val magenta          = Color.MAGENTA
    val maroon           = Color.MAROON
    val mediumaquamarine = Color.MEDIUMAQUAMARINE
    val mediumblue       = Color.MEDIUMBLUE
    val mediumorchid     = Color.MEDIUMORCHID
    val mediumpurple     = Color.MEDIUMPURPLE
    val mediumseagreen   = Color.MEDIUMSEAGREEN
    val mediumslateblue  = Color.MEDIUMSLATEBLUE
    val mediumspringgreen = Color.MEDIUMSPRINGGREEN
    val mediumturquoise  = Color.MEDIUMTURQUOISE
    val mediumvioletred  = Color.MEDIUMVIOLETRED
    val midnightblue     = Color.MIDNIGHTBLUE
    val mintcream        = Color.MINTCREAM
    val mistyrose        = Color.MISTYROSE
    val moccasin         = Color.MOCCASIN
    val navajowhite      = Color.NAVAJOWHITE
    val navy             = Color.NAVY
    val oldlace          = Color.OLDLACE
    val olive            = Color.OLIVE
    val olivedrab        = Color.OLIVEDRAB
    val orange           = Color.ORANGE
    val orangered        = Color.ORANGERED
    val orchid           = Color.ORCHID
    val palegoldenrod    = Color.PALEGOLDENROD
    val palegreen        = Color.PALEGREEN
    val paleturquoise    = Color.PALETURQUOISE
    val palevioletred    = Color.PALEVIOLETRED
    val papayawhip       = Color.PAPAYAWHIP
    val peachpuff        = Color.PEACHPUFF
    val peru             = Color.PERU
    val pink             = Color.PINK
    val plum             = Color.PLUM
    val powderblue       = Color.POWDERBLUE
    val purple           = Color.PURPLE
    val red              = Color.RED
    val rosybrown        = Color.ROSYBROWN
    val royalblue        = Color.ROYALBLUE
    val saddlebrown      = Color.SADDLEBROWN
    val salmon           = Color.SALMON
    val sandybrown       = Color.SANDYBROWN
    val seagreen         = Color.SEAGREEN
    val seashell         = Color.SEASHELL
    val sienna           = Color.SIENNA
    val silver           = Color.SILVER
    val skyblue          = Color.SKYBLUE
    val slateblue        = Color.SLATEBLUE
    val slategray        = Color.SLATEGRAY
    val snow             = Color.SNOW
    val springgreen      = Color.SPRINGGREEN
    val steelblue        = Color.STEELBLUE
    val tan              = Color.TAN
    val teal             = Color.TEAL
    val thistle          = Color.THISTLE
    val tomato           = Color.TOMATO
    val turquoise        = Color.TURQUOISE
    val violet           = Color.VIOLET
    val wheat            = Color.WHEAT
    val white            = Color.WHITE
    val whitesmoke       = Color.WHITESMOKE
    val yellow           = Color.YELLOW
    val yellowgreen      = Color.YELLOWGREEN

    val hi = 0xE0        // high intensity
    val md = 0x80        // mid intensity

    //::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    /** The `Randi` inner case class provides an implementation of the 'igen'
     *  method eliminating any dependency of the 'random' package.
     *  @param a  the lower bound (inclusive)
     *  @param b  the upper bound (inclusive)
     */
    case class Randi (a: Int, b: Int)
    {
         private val rng = new Random ()
         val diff = b - a + 1
         def igen: Int = a + rng.nextInt (diff)
    } // Randi class

    //::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    /** Generate (somewhat randomly) an RGB color, based on the object's id.
     *  It randomizes off of 27/25 base colors that include all the 16 colors
     *  specified by W3C, except silver.
     *  @see http://www.w3.org/TR/css3-color
     *  @param id    the identifier/index for some object
     *  @param full  whether to use all 27 color bases or exclude black and white
     */
    def randomColor (id: Int, full: Boolean = true): Color =
    {
        val r     = Randi (0x0, 0x1F)
        val bases = if (full) 27 else 25
        id % bases match {
            case  0 => new Color ((hi+r.igen)/255, r.igen/255, r.igen/255, 1)             // RED
            case  1 => new Color ((hi+r.igen)/255, (md+r.igen)/255, r.igen/255, 1)        // orange
            case  2 => new Color ((hi+r.igen)/255, (md+r.igen)/255, (md+r.igen)/255, 1)   // peach

            case  3 => new Color ((hi+r.igen)/255, (hi+r.igen)/255, r.igen/255, 1)        // YELLOW
            case  4 => new Color ((hi+r.igen)/255, (hi+r.igen)/255, (md+r.igen)/255, 1)   // light yellow
            case  5 => new Color ((md+r.igen)/255, (md+r.igen)/255, r.igen/255, 1)        // OLIVE

            case  6 => new Color ((md+r.igen)/255, (hi+r.igen)/255, r.igen/255, 1)        // Lime green
            case  7 => new Color ((md+r.igen)/255, (hi+r.igen)/255, (md+r.igen)/255, 1)   // light green
            case  8 => new Color (r.igen/255,      (hi+r.igen)/255, r.igen/255, 1)        // LIME

            case  9 => new Color (r.igen/255,      (hi+r.igen)/255, (md+r.igen)/255, 1)   // sea green
            case 10 => new Color (r.igen/255,      (md+r.igen)/255, r.igen/255, 1)        // GREEN
            case 11 => new Color (r.igen/255,      (md+r.igen)/255, (md+r.igen)/255, 1)   // TEAL

            case 12 => new Color ((md+r.igen)/255, (hi+r.igen)/255, (hi+r.igen)/255, 1)   // AQUA
            case 13 => new Color (r.igen/255,      (hi+r.igen)/255, (hi+r.igen)/255, 1)   // cyan
            case 14 => new Color (r.igen/255,      (md+r.igen)/255, (hi+r.igen)/255, 1)   // blue green

            case 15 => new Color ((md+r.igen)/255, (md+r.igen)/255, (hi+r.igen)/255, 1)   // periwinkle
            case 16 => new Color (r.igen/255,      r.igen/255,    (hi+r.igen)/255, 1)     // BLUE
            case 17 => new Color (r.igen/255,      r.igen/255,    (md+r.igen)/255, 1)     // NAVY

            case 18 => new Color ((md+r.igen)/255, r.igen/255,    (hi+r.igen)/255, 1)     // violet
            case 19 => new Color ((md+r.igen)/255, r.igen/255,    (md+r.igen)/255, 1)     // PURPLE
            case 20 => new Color ((hi+r.igen)/255, r.igen/255,    (hi+r.igen)/255, 1)     // FUCHSIA

            case 21 => new Color ((hi+r.igen)/255, (md+r.igen)/255, (hi+r.igen)/255, 1)   // pink
            case 22 => new Color ((md+r.igen)/255, r.igen/255,    r.igen/255, 1)          // MAROON
            case 23 => new Color ((hi+r.igen)/255, r.igen/255,    (md+r.igen)/255, 1)     // magenta

            case 24 => new Color ((md+r.igen)/255, (md+r.igen)/255, (md+r.igen)/255, 1)   // GRAY
            case 25 => new Color (r.igen/255,      r.igen/255,    r.igen/255, 1)          // BLACK
            case 26 => new Color ((hi+r.igen)/255, (hi+r.igen)/255, (hi+r.igen)/255, 1)   // WHITE
        } // match
    } // randomColor

    //::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    /*  Based upon an object's id, generate (somewhat randomly) a color.
     *  @param id  the identifier/index for some object
     *
    def randomColor (id: Int) =
    {
        rc = (hi+r.igen)/255, (0x1F, 0xF7)
        id % 3 match {
            case 0 => new Color (255, (31 * id) % 224, 0, 1)
            case 1 => new Color (255, 0, (31 * id) % 256, 1)
            case _ => new Color (32 + (31 * id) % 224, 0, 255, 1)
        } // match
    } // randomColor
     */

} // Colors object


//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
/** The `ColorsTest` object is used to test the `Colors` object.
 */
object ColorsTest extends App
{
    import Colors.randomColor
    for (i <- 0 until 20) println ("color " + i + " = " + randomColor (i))

} // ColorsTest object

