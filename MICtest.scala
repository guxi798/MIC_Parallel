import scalation.analytics.Probability._

import scalation.linalgebra.MatrixD
import scalation.linalgebra.MatrixC
import scala.collection.mutable.MutableList
import scalation.linalgebra.VectorD
import scalation.linalgebra_gen.Matrices.MatrixI
import scalation.linalgebra_gen.MatrixN
import scala.math._
import scala.io.Source
import java.io._
import Array._

object MICtest extends App{
val beginning = System.currentTimeMillis();
val filename = "soybean.tab"

    val numOfRow = Source.fromFile(filename).getLines().size
    var col_count = 1
    var element = Source.fromFile(filename).getLines().next
    //println("element = "+ element)

      for(temp <- element)
      {
        if (temp == '\t'){
          col_count += 1
        }//if
      }//for temp <- element
 //println("col_count = " + col_count)
col_count = (col_count/100).toInt
 var mymatrix = ofDim[Double](numOfRow, col_count)
 //println("mymatrix = " + mymatrix.deep)

 //assignment
 var row = 0
    for (line <- Source.fromFile(filename).getLines()) {
      //println("line" + line)
        for(i <- 0 until col_count){
          mymatrix(row)(i) = (line.stripSuffix("\n").split("\t").toArray.deep)(i).toString.toDouble
        //mymatrix(row)(i+1) = (line.stripSuffix("\n").split(",").toArray.deep)(i+1).toString.toDouble
    }//for column
        row += 1
      }//for
 //println("mymatrix = " + mymatrix.deep)

  var D = new MatrixD(mymatrix)
  //println("D = " + D)
  //println("D(0) = " + D.col(0))
  //println("D(1) = " + D.col(1))
  var dcount = 0
  val writer = new PrintWriter(new File("test.txt" ))
  for(i <- 0 until col_count)//outer
  {
    for(j <- 0 until col_count)//inner
    {
println("Column (" + i + "," + j+")")
      if(i==j)
      {
        writer.write("1")
        if(j != col_count-1){
        writer.write("\t")
        }
        }//if(i==j)
      
      else{
        dcount+=1
        if(j<i)
        {
        
        }//if(j<i)
        else{
          var d1 = new MatrixD(Array(D.col(i), D.col(j))).t
  //println("d1 = "+d1)
  val y = 2
  //println(y)
  var m = new MIC
  writer.write((m.ApproxMaxMI(d1, 2, y, 1*3)(2)).toString)
  if(j != col_count-1){
  writer.write("\t")
        }//else
  }
  //println(m.ApproxMaxMI(d1, 2, y, 1*3)(2))
      }
    }//inner for
    writer.write("\n")
  }//for
  val finish = System.currentTimeMillis();
writer.close()
println("Running time = " + (beginning-finish))
}
