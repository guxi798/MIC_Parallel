import scalation.analytics.Probability._

import scalation.linalgebra.MatrixD
import scalation.linalgebra.MatrixC
import scala.collection.mutable.MutableList
import scalation.linalgebra.VectorD
import scalation.linalgebra_gen.Matrices.MatrixI
import scalation.linalgebra_gen.MatrixN
import scala.math._
import scala.collection.mutable.HashMap

class MIC {//(D: MatrixD, x:Int , y:Int, k_hat:Int){
  def ApproxMaxMI(D: MatrixD, x:Int , y:Int, k_hat:Int): Array[Double] = {
    
    var Dy = EquipartitionYAxis(D, y)
    var Dx = SortInIncreasingOrderByValue(Dy, 0)
    var Q = Dx.col(2).apply.map(_.toInt)
   
    return OptimizeXAxis(Dx.selectCols(Array(0,1)), Q, x, k_hat)
  }
  
  def EquipartitionYAxis(D: MatrixD, y: Int): MatrixD = 
  {
    var D_Y = SortInIncreasingOrderByValue(D, 1) //sorted increasing by Y
    var i:Int = 0
    var currRow:Int = 1
    var desiredRowSize:Int = D_Y.dim1/y
    var no = 0
    var Q = new Array[Int](D.dim1)
    
    while(i < D_Y.dim1)
    { var s = 0
      
     var temp = D_Y(i,1)
     for(j <- 0 until D_Y.dim1)
     {
       if(D_Y(j,1)==temp)
       {
         s +=1
       }//if(t==temp)
     }//for
 
    if(no!=0 && ((no + s -desiredRowSize).abs>= (no - desiredRowSize).abs) ){
      currRow += 1
      no = 0
      desiredRowSize = (D_Y.dim1-i)/(y-currRow+1)
    }//if
    
    for(x <- i until i+s)
    {
      Q(x) = currRow
    }//for
    i+=s
    no +=s
    }//while
  
    var i_t = 0
	var Djoin = Array[VectorD]() //join D_Y+Q
	
	  for(e <- D_Y){
	    var t =  new VectorD(e :+ Q(i_t).toDouble)
	    i_t += 1
	    Djoin = Djoin :+ t
	  }
	  //println("EquipartitionYAxis end")
	
	 new MatrixD(Djoin)
  }//EquipartitionYAxis

  
  def SortInIncreasingOrderByValue(D:MatrixD, i: Int): MatrixD = {
    // get sorted D by X or Y (0,1)
    var a = Array[VectorD]()
    for(e <- D){
      var t =  new VectorD(e)
      a = a :+ t
      }//for
    
    var Dsort = new MatrixD(a.sortBy(_(i)))
    Dsort
  }
  
  def GetClumpsPartition (D: MatrixD, Q: Array[Int]): MatrixD= {
	  /*
	   * returns the map P defining the clumps for the set D
	   * constrain that points with the same x-value in the same clump
	   */
    var Q1:Array[Int]=Q.clone
	  var i = 0
	  var s = 0
	  while(i < D.dim1){
	    for(j <- i+1 until D.dim1){
	      if(D(i,0) == D(j,0)){
	        s += 1
	        if(Q1(i) != Q1(j)){
	          Q1(j) = Q1(i)
	        }
	      } //if
	    } //for j
	    i += (s+1)
	  } //while i
	  
	  var P = new Array[Int](D.dim1) 
	  i = 1
	  P(0) = i
	  for(j <- 1 until D.dim1){
	    if(Q1(j) != Q1(j-1)){
	      i += 1
	    }
	      P(j) = i
	  }
	  
	  var i_t = 0
	  var Djoin = Array[VectorD]() //join D_Y+Q 
	
	  for(e <- D){
	    var t =  new VectorD(e :+ P(i_t).toDouble)
	    i_t += 1
	    Djoin = Djoin :+ t
	  }

	  //println("GetClimpsPartition end")
	  
	new MatrixD(Djoin)
	  
	}//GetClumpsPartition

  
  def OptimizeXAxis(D0: MatrixD, Q: Array[Int], x: Int, k_hat:Int): Array[Double] = {

    var Djoin = GetClumpsPartition(D0, Q)
    var D = Djoin.selectCols(Array(0,1))
    var P = Djoin.col(2).apply.map(_.toInt)
    //println("Plength = "+P.length)
    var C = Array(0)
    for(i <- 1 until P.length){
      if(P(i) != P(i-1)){
        C = C :+ i
      }
    }
    C = C :+ (P.length-1)
    //println("C = "+C.deep)
    var k = C.length-1
    
    //println("x = "+x+"k = "+k)
 
    var Pmat = Array.ofDim[Int](k+1, x+1, k) 
    var Imat = Array.ofDim[Double](k+1, x+1)
   
    for(t <- 2 to k){
      /* Find the optimal partitions of size 2 */
      var Fmax = -1.
      for(s <- 1 until t){
        var Fpart = HP(D, Array(C(0), C(s), C(t))) - HPQ(D, Array(C(0), C(s), C(t)), Q)
        if(Fpart > Fmax){
          Fmax = Fpart
          Pmat(t)(2) = Array(C(0), C(s), C(t))
          Imat(t)(2) = Fmax + HQ(D, Q)
        }// if max
      }// for s
    }// for t
    
    for(l <- 3 to x; t <- l to k){
      var Fmax = -1.0
      for(s <- l-1 until t){
        var Fpart = (s/t) * (Imat(s)(l-1) - HQ(D, Q)) - ((t-s)/t * HPQ(D, Array(C(0), C(s), C(t)), Q))
          HP(D, Array(C(0), C(s), C(t))) - HPQ(D, Array(C(0), C(s), C(t)), Q)
          if(Fpart > Fmax){
            Fmax = Fpart
            Pmat(t)(l) = Insert(Pmat(s)(l-1), C(t))
            Imat(t)(l) = HQ(D, Q) + HP(D, Pmat(t)(l)) - HPQ(D, Pmat(t)(l), Q)
        }// if max
      }// for s
    } // for l and t
    /*
    for(i <- (k+1) until x){
      Pmat(k)(i) = Pmat(k)(k)
      Imat(k)(i) = Imat(k)(k)
    }
    * 
    */
   Imat(k)
  }
  
  def Insert(P: Array[Int], t:Int): Array[Int] = {
    for(i <- 0 until P.length){
      if(P(i) == t){
        return P
      }else if(P(i) > t){
        return Array.concat(P.slice(0, i), t +: P.slice(i, P.length))
      }
    }
    P :+ t
  }
  
  def HP(D: MatrixD, P: Array[Int]): Double = {
    var px = new VectorD(P.length-1) //, Array[Double]())
    for(i <- 2 until P.length){
      px(i-1) = (P(i).toDouble - P(i-1).toDouble)/D.dim1.toDouble 
    }
    var H = entropy(px)
    //println("HP = " + H)
    H
    
  }
  
  def HQ(D: MatrixD, Q: Array[Int]): Double = {
    var py = new VectorD(Q.reduce(max)) //, Array[Double]())
    for(j <- 0 until Q.length){
      py(Q(j)-1) += 1.0
    }
    var py1 = py / D.dim1.toDouble
    var H = entropy(py1)
    //println("HQ = " + H)
    H
  }
  
  def HPQ(D: MatrixD, P: Array[Int], Q: Array[Int]): Double = {
    var pxy = new MatrixD(P.length-1, Q.reduce(max)) //, Array[Array[Double]]())
    //println("P = "+P.deep)
    //println("Q = "+Q.deep)
    for(i <- 1 until P.length){
      for(j <- P(i-1) until P(i)){
          pxy(i-1, Q(j)-1) += 1.0
     }//for j
    } //for i
    //println("pxy = "+pxy)
    pxy = pxy / (D.dim1 * D.dim2).toDouble
    //println("pxy = "+pxy)
    var HPQ = entropy(pxy)
    //println("HPQ = " + HPQ)
    HPQ
  }
  
    /*
  def ApproxCharacteristicMatrix(D:MatrixD, B:Int, c:Int): MatrixD = {
    var Dswap = new MatrixD(D.dim1, D.dim2)
    Dswap.col(0) += D.col(1)
    Dswap.col(1) += D.col(0)
    
    var Ixy = new MatrixD(B/2+1, B/2+1)
    var Mxy = new MatrixD(B/2+1, B/2+1)
    
    for(y <- 2 to B/2){
      var x = B/y
      var I = ApproxMaxMI(D, x, y, c*x)
      var Iswap = ApproxMaxMI(Dswap, x, y, c*x)
      for(i <- 0 until I.length){
        Ixy(x)(y) = max(I(i), Iswap(i))
        Mxy(x)(y) = Ixy(x)(y) / min(log(x), log(y))
      }
    }
    Mxy
  }
  * 
  */
} // class MIC