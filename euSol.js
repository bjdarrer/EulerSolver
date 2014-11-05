/**
 * Simulating and visualizing reaction-diffusion systems
 * based on Subquantum Kinetics Model G
 * Programmed by David Jonsson, Athens January 2014
  *                             Stockholm March, April 2014
 * @constructor
 */
function EuSol(anchor) {
  if (anchor) {
  }
  else
    anchor = window;
  var sqkMgthis = this;

  /**
   * Size setup parameters.
   */
/** @const */  this.size = {
    'rows' : 500,
    'columns' : 500
  }; 
  
  /**
   * Reaction parameters.
   * @const
   * @ctype {object}
   */
/** @const */  this.k = {
    "1":  1,
    "_1": 1,
    "2":  1,
    "_2": 1,
    "3":  1,
    "_3": 0,
    "4":  1,
    "_4": 0,
    "5":  1,
    "_5": 0
  };
  
  /**
   * Diffusion parameters.
   * @const
   * @ctype {object}
   */
/** @const */  this.diffusion = {
    "DG" : 1,
    "DX" : 1,
    "DY" : 1
  };
  
  /**
   * Initial static concentrations.
   * @const
   * @ctype {object}
   */
/** @const */  var concentrations = /*new Float64Array(8);*/ {
    "A": 1,
    "B": 1,
    "G": 1,
    "X": 1,
    "Y": 1,
    "Z": 1,
    "Omega": 1
  };
  /**
   * Initial static concentrations.
   * @const
   * @ctype {object}
   */
/** @const */  this.dt = 0.5;

//Constants during a run (usually) 
//  const x = 0;

  // Configurations, user changeable  
  // Canvas id
  this.plotDOMid = "c";
  this.statusDOMclass = "simulationStatus";
  this.runSimulation = true;
  this.mask = "none";
  this.maskParam = 1;
  this.calcStepIterationDOM = document.getElementsByClassName('calcStepIteration')[0];
  this.calcStepDurationDOM = document.getElementsByClassName('calcStepDuration')[0];

  // End of configurations, user changeable
  
  this.initializeArrays = function() {
    //Initial dynamic concentrations. Format GXYGXY ... 64 bit double array might be needed
    this.GXYarr = new Array(2);
    if (window.Float64Array) {
      this.GXYarr = [new Float64Array(3 * this.size.rows * this.size.columns),
      new Float64Array(3 * this.size.rows * this.size.columns)];
    } else {
      this.GXYarr = [new Array(3 * this.size.rows * this.size.columns),
      new Array(3 * this.size.rows * this.size.columns)];
    }
    for (var index = 0; index < this.GXYarr[0].length; index += 3) {
      this.GXYarr[0][index] = concentrations.G;
      this.GXYarr[1][index] = concentrations.G;
      this.GXYarr[0][index + 1] = concentrations.X;
      this.GXYarr[1][index + 1] = concentrations.X;
      this.GXYarr[0][index + 2] = concentrations.Y;
      this.GXYarr[1][index + 2] = concentrations.Y;
    }
    this.iteration = 0;
    this.shownIteration = 0;
  };
  
  this.initializeArrays();
  this.GXYbuffer = 0;
    
  // Laplacian numerial operator âˆ‡^2, wheight 1 on adjacents and Â½ on diagonals 
  var laplace2dDiagonal = function(V, index) {
    return -6*V[index] +
            V[index + 3] + V[index - 3] + V[index + 3 * sqkMgthis.size["columns"]] + V[index - 3 * sqkMgthis.size["columns"]] + //closest neighbours
            0.5*(V[index + 3 * sqkMgthis.size.columns + 3] + V[index + 3 * sqkMgthis.size.columns - 3] + V[index - 3 * sqkMgthis.size.columns + 3] + V[index - 3 * sqkMgthis.size.columns - 3]); //diagonal neighbours
  };
  var laplace2dAdjacentOnly = function(V, index) {
    return -4*V[index] +
            V[index + 3] + V[index - 3] + V[index + 3 * sqkMgthis.size["columns"]] + V[index - 3 * sqkMgthis.size["columns"]];  //closest neighbours
  };
//  this.laplace = laplace2dAdjacentOnly;
  this.laplace = laplace2dDiagonal;

//Masking, boundary conditions
  this.masking = function(i, j, mask) {
    var im = (i / this.size.columns - 1/2);
    var jm = (j / 3 / this.size.rows - 1/2);
    switch (mask) {
      case "none":
        return true;
      case "elliptic":
        return im * im + jm * jm < sqkMgthis.maskParam;
      case "vertical":
        return im < sqkMgthis.maskParam;
      case "horizontal":
        return jm < sqkMgthis.maskParam;
      case "rhombic":
        return im + jm < sqkMgthis.maskParam;
      case "triangle":
        return im - jm > sqkMgthis.maskParam;
      case "innerElliptic":
        return im * im + jm * jm > sqkMgthis.maskParam;
      default:
        return true;
    }
  }  

//Concentrations, not complete
  this.gxyConcentrations = function(arr3, readIndex, readBuffer, writeBuffer) {
    var G = readBuffer[3 * this.size.columns * i + j];
    var X = readBuffer[3 * this.size.columns * i + j + 1];
    var Y = readBuffer[3 * this.size.columns * i + j + 2];
    G += (this.diffusion['DG'] * this.laplace(arr3[readIndex], 3 * this.size.columns * i + j) - (sqkMgthis.k["_1"] + sqkMgthis.k["2"]) * G + sqkMgthis.k["_2"] * X + sqkMgthis.k["1"] * concentrations['A']) * sqkMgthis.dt;

    X += (this.diffusion['DX'] * this.laplace(arr3[readIndex], 3 * this.size.columns * i + j + 1) + sqkMgthis.k["2"] * G - (sqkMgthis.k["_2"] + sqkMgthis.k["3"] * concentrations['B'] + sqkMgthis.k["5"]) * X +
      sqkMgthis.k["_3"] * concentrations['Z'] * Y - sqkMgthis.k["_4"] * X * X * X + sqkMgthis.k["4"] * X * X * Y + sqkMgthis.k["_5"] * concentrations['Î©']) * sqkMgthis.dt;

    Y += (this.diffusion['DY'] * this.laplace(arr3[readIndex], 3 * this.size.columns * i + j + 2) + sqkMgthis.k["3"] * concentrations['B'] * X - sqkMgthis.k["_3"]*concentrations['Z']*Y + 
      sqkMgthis.k["_4"] * X * X * X - sqkMgthis.k["4"]*X*X * Y ) * sqkMgthis.dt;
    return [G, X, Y];
  }  

  //Nondimensionalized potentials
//  âˆ‚G
//  âˆ‚t = âˆ‡2G âˆ’ qG + gX + a (7a)
//  âˆ‚X
//  âˆ‚t = dxâˆ‡  2X + pG âˆ’ (1 + b)X + uY âˆ’ sX3 + X
//  2Y + w (7b)
//  âˆ‚Y
//  âˆ‚t = dyâˆ‡
//  2Y + bX âˆ’ uY + sX3 âˆ’ X
//  2Y (7c)
  
  this.gxyPotentials = function(arr3, readIndex) {
    var G = arr3[readIndex][3 * this.size.columns * i + j];
    var X = arr3[readIndex][3 * this.size.columns * i + j + 1];
    var Y = arr3[readIndex][3 * this.size.columns * i + j + 2];
    G += (this.laplace(arr3[readIndex], 3 * this.size.columns * i + j) - (sqkMgthis.k["_1"] + sqkMgthis.k["2"]) / (sqkMgthis.k["_2"] + sqkMgthis.k["5"]) * G + 1/(1-sqkMgthis.k["5"]/sqkMgthis.k["_2"]) * X + sqkMgthis.k["1"] * Math.sqrt(sqkMgthis.k["4"])*Math.pow(sqkMgthis.k["_2"] +sqkMgthis.k["5"], 1.5) * concentrations['A'] ) * dt;

    X += (this.diffusion['DX'] / this.diffusion['DG'] * this.laplace(arr3[readIndex], 3 * this.size.columns * i + j + 1) + sqkMgthis.k["2"]/(sqkMgthis.k["_2"] - sqkMgthis.k["5"]) * G - 
      + (1 + sqkMgthis.k["_2"] / (sqkMgthis.k["_2"] + sqkMgthis.k["5"]) ) * concentrations['B'] * X + //Done until here
      sqkMgthis.k["_3"] * concentrations['Z'] * Y - sqkMgthis.k["_4"] * X * X * X + sqkMgthis.k["4"] * X * X * Y + sqkMgthis.k["_5"] * concentrations['Î©']) * dt;

    Y += (this.diffusion['DY'] * this.laplace(arr3[readIndex], 3 * this.size.columns * i + j + 2) + sqkMgthis.k["3"] * concentrations['B'] * X - sqkMgthis.k["_3"]*concentrations['Z']*Y + 
      sqkMgthis.k["_4"] * X * X * X - sqkMgthis.k["4"]*X*X * Y ) * dt;
    return [G, X, Y];
  }  

  //The 3 valued return function in the iteration
  this.gxy=this.gxyConcentrations;

//  Reaction diffusion functions   
//  Iterate, leave the boundary with i,j=1 ...
  this.gxyStep = function(arr3) {
    var functionStart = window.performance.now();
    var readIndex, writeIndex;
    var conc3 = [];
    this.extremesGXY = {
        "gMin": Infinity,
        "gMax": -Infinity,
        "xMin": Infinity,
        "xMax": -Infinity,
        "yMin": Infinity,
        "yMax": -Infinity
      };
    if (this.GXYbuffer === 0) {
      this.GXYbuffer = 1;
      readIndex = 0;
      writeIndex = 1;
    } else {
      this.GXYbuffer = 0;
      readIndex = 1;
      writeIndex = 0;      
    }
    
    var readBuffer = arr3[readIndex];
    var writeBuffer = arr3[writeIndex];

    var  o = 0;
    for (j = 3; j < 3 * this.size.rows - 3; j += 3 ) {
      for (i = 1; i < this.size.columns - 2; i++ ) {
        if (this.masking(i, j, this.mask)) {
//          var G = arr3[readIndex][3 * this.size.columns * i + j];
//          var X = arr3[readIndex][3 * this.size.columns * i + j + 1];
//          var Y = arr3[readIndex][3 * this.size.columns * i + j + 2];

//          G += (this.diffusion['DG'] * this.laplace(arr3[readIndex], 3 * this.size.columns * i + j) - (sqkMgthis.k["_1"] + sqkMgthis.k["2"]) * G + sqkMgthis.k["_2"] * X + sqkMgthis.k["1"] * concentrations['A']) * dt;
//
//          X += (this.diffusion['DX'] * this.laplace(arr3[readIndex], 3 * this.size.columns * i + j + 1) + sqkMgthis.k["2"] * G - (sqkMgthis.k["_2"] + sqkMgthis.k["3"] * concentrations['B'] + sqkMgthis.k["5"]) * X +
//            sqkMgthis.k["_3"] * concentrations['Z'] * Y - sqkMgthis.k["_4"] * X * X * X + sqkMgthis.k["4"] * X * X * Y + sqkMgthis.k["_5"] * concentrations['Î©']) * dt;
//
//          Y += (this.diffusion['DY'] * this.laplace(arr3[readIndex], 3 * this.size.columns * i + j + 2) + sqkMgthis.k["3"] * concentrations['B'] * X - sqkMgthis.k["_3"]*concentrations['Z']*Y + 
//            sqkMgthis.k["_4"] * X * X * X - sqkMgthis.k["4"]*X*X * Y ) * dt;
          conc3 = sqkMgthis.gxy(arr3, readIndex, readBuffer);

//          arr3[writeIndex][3 * this.size.columns * i + j] = G;
//          arr3[writeIndex][3 * this.size.columns * i + j + 1] = X;
//          arr3[writeIndex][3 * this.size.columns * i + j + 2] = Y;

          writeBuffer[3 * this.size.columns * i + j] = conc3[0];
          writeBuffer[3 * this.size.columns * i + j + 1] = conc3[1];
          writeBuffer[3 * this.size.columns * i + j + 2] = conc3[2];

          sqkMgthis.extremesGXY.gMin = Math.min(sqkMgthis.extremesGXY.gMin , conc3[0]);
          sqkMgthis.extremesGXY.gMax = Math.max(sqkMgthis.extremesGXY.gMax , conc3[0]);
          sqkMgthis.extremesGXY.xMin = Math.min(sqkMgthis.extremesGXY.xMin , conc3[1]);
          sqkMgthis.extremesGXY.xMax = Math.max(sqkMgthis.extremesGXY.xMax , conc3[1]);
          sqkMgthis.extremesGXY.yMin = Math.min(sqkMgthis.extremesGXY.yMin , conc3[2]);
          sqkMgthis.extremesGXY.yMax = Math.max(sqkMgthis.extremesGXY.yMax , conc3[2]);
          
           /*
          if (sqkMgthis.extremesGXY.gMin > conc3[0]) {
            sqkMgthis.extremesGXY.gMin = conc3[0];
          } else if (sqkMgthis.extremesGXY.gMax < conc3[0]) {
            sqkMgthis.extremesGXY.gMax = conc3[0];
          }
          if (sqkMgthis.extremesGXY.xMin > conc3[1]) {
            sqkMgthis.extremesGXY.xMin = conc3[1];
          } else if (sqkMgthis.extremesGXY.xMax < conc3[1]) {
            sqkMgthis.extremesGXY.xMax = conc3[1];
          }
          if (sqkMgthis.extremesGXY.yMin > conc3[2]) {
            sqkMgthis.extremesGXY.yMin = conc3[2];
          } else if (sqkMgthis.extremesGXY.yMax < conc3[2]) {
            sqkMgthis.extremesGXY.yMax = conc3[2];
          }*/
        }
      };
    };    
    this.iteration++;
    this.calcStepDuration = (performance.now() - functionStart).toFixed(3);
//    console.log("iteration ", this.iteration, " done in ", performance.now() - functionStart, " milliseconds");
  };
  
  //Plot the data in a 2D Canvas, (G, X, Y) mapped to (R, G, B)
  this.plot2D = function(arr3) {
    var c = document.getElementById(sqkMgthis.plotDOMid);
    var ctx = c.getContext("2d");
    var imgData = ctx.getImageData(0, 0, c.width, c.height);
    var j = 0; // Image pixel color index
    for (var i = 0; i < arr3[sqkMgthis.GXYbuffer].length; i += 3) {
      imgData.data[j++] = (sqkMgthis.extremesGXY.gMax === sqkMgthis.extremesGXY.gMin) ? 0 : 255 * (arr3[sqkMgthis.GXYbuffer][i]     - sqkMgthis.extremesGXY.gMin) / (sqkMgthis.extremesGXY.gMax - sqkMgthis.extremesGXY.gMin);
      imgData.data[j++] = (sqkMgthis.extremesGXY.xMax === sqkMgthis.extremesGXY.xMin) ? 0 : 255 * (arr3[sqkMgthis.GXYbuffer][i + 1] - sqkMgthis.extremesGXY.xMin) / (sqkMgthis.extremesGXY.xMax - sqkMgthis.extremesGXY.xMin);
      imgData.data[j++] = (sqkMgthis.extremesGXY.yMax === sqkMgthis.extremesGXY.yMin) ? 0 : 255 * (arr3[sqkMgthis.GXYbuffer][i + 2] - sqkMgthis.extremesGXY.yMin) / (sqkMgthis.extremesGXY.yMax - sqkMgthis.extremesGXY.yMin);
      imgData.data[j++] = 255;
    }
    ctx.putImageData(imgData, 0, 0);
    
    this.calcStepIterationDOM.innerHTML = this.iteration;
    this.calcStepDurationDOM.innerHTML = this.calcStepDuration;
  };
  
  this.loopSimulation = function () {
    if (sqkMgthis.runSimulation) {
      sqkMgthis.gxyStep(sqkMgthis.GXYarr);
      setTimeout(sqkMgthis.loopSimulation, 0);
    }
  };
  
  this.loopVisualization = function () {
    if (sqkMgthis.shownIteration < sqkMgthis.iteration) {
      sqkMgthis.plot2D(sqkMgthis.GXYarr, sqkMgthis.plotDOMid);
      sqkMgthis.shownIteration = sqkMgthis.iteration;
    }
    sqkMgthis.loopRequestAnimationFrameID = requestAnimationFrame(sqkMgthis.loopVisualization);
  };
  
  this.loadNewParamsFromView = function() {
  //Get values from view "changeable values"
    var newParamsClass = document.getElementsByClassName("newGparams")[0];
    concentrations.A = parseFloat((newParamsClass.getElementsByClassName("A")[0]).value);
    concentrations.B = parseFloat((newParamsClass.getElementsByClassName("B")[0]).value);
    concentrations.G = parseFloat((newParamsClass.getElementsByClassName("G")[0]).value);
    concentrations.X = parseFloat((newParamsClass.getElementsByClassName("X")[0]).value);
    concentrations.Y = parseFloat((newParamsClass.getElementsByClassName("Y")[0]).value);
    concentrations.Z = parseFloat((newParamsClass.getElementsByClassName("Z")[0]).value);
    concentrations.Omega = parseFloat((newParamsClass.getElementsByClassName("Omega")[0]).value);

    this.diffusion.DG = parseFloat((newParamsClass.getElementsByClassName("DG")[0]).value);
    this.diffusion.DX = parseFloat((newParamsClass.getElementsByClassName("DX")[0]).value);
    this.diffusion.DY = parseFloat((newParamsClass.getElementsByClassName("DY")[0]).value);

    this.k["1"]  = parseFloat((newParamsClass.getElementsByClassName("k1")[0]).value);
    this.k["_1"] = parseFloat((newParamsClass.getElementsByClassName("k_1")[0]).value);
    this.k["2"]  = parseFloat((newParamsClass.getElementsByClassName("k2")[0]).value);
    this.k["_2"] = parseFloat((newParamsClass.getElementsByClassName("k_2")[0]).value);
    this.k["3"]  = parseFloat((newParamsClass.getElementsByClassName("k3")[0]).value);
    this.k["_3"] = parseFloat((newParamsClass.getElementsByClassName("k_3")[0]).value);
    this.k["4"]  = parseFloat((newParamsClass.getElementsByClassName("k4")[0]).value);
    this.k["_4"] = parseFloat((newParamsClass.getElementsByClassName("k_4")[0]).value);
    this.k["5"]  = parseFloat((newParamsClass.getElementsByClassName("k5")[0]).value);
    this.k["_5"] = parseFloat((newParamsClass.getElementsByClassName("k_5")[0]).value);

    this.dt = parseFloat((newParamsClass.getElementsByClassName("dt")[0]).value);

    this.size.rows = parseFloat((newParamsClass.getElementsByClassName("height")[0]).value);
    this.size.columns = parseFloat((newParamsClass.getElementsByClassName("width")[0]).value);

    this.updateParamsInView("loadedGparams");
  };
    
  this.updateParamsInView = function(paramGroup) {    
  //Update the web view list of params        
    var updateClass = document.getElementsByClassName(paramGroup)[0]; 
    (updateClass.getElementsByClassName("A")[0]).setAttribute("value", concentrations.A);
    (updateClass.getElementsByClassName("B")[0]).setAttribute("value", concentrations.B);
    (updateClass.getElementsByClassName("G")[0]).setAttribute("value", concentrations.G);
    (updateClass.getElementsByClassName("X")[0]).setAttribute("value", concentrations.X);
    (updateClass.getElementsByClassName("Y")[0]).setAttribute("value", concentrations.Y);
    (updateClass.getElementsByClassName("Z")[0]).setAttribute("value", concentrations.Z);
    (updateClass.getElementsByClassName("Omega")[0]).setAttribute("value", concentrations.Omega);
    
    (updateClass.getElementsByClassName("DG")[0]).setAttribute("value", this.diffusion.DG);
    (updateClass.getElementsByClassName("DX")[0]).setAttribute("value", this.diffusion.DX);
    (updateClass.getElementsByClassName("DY")[0]).setAttribute("value", this.diffusion.DY);

    (updateClass.getElementsByClassName("k1")[0]).setAttribute("value", this.k["1"]);
    (updateClass.getElementsByClassName("k_1")[0]).setAttribute("value", this.k["_1"]);
    (updateClass.getElementsByClassName("k2")[0]).setAttribute("value", this.k["2"]);
    (updateClass.getElementsByClassName("k_2")[0]).setAttribute("value", this.k["_2"]);
    (updateClass.getElementsByClassName("k3")[0]).setAttribute("value", this.k["3"]);
    (updateClass.getElementsByClassName("k_3")[0]).setAttribute("value", this.k["_3"]);
    (updateClass.getElementsByClassName("k4")[0]).setAttribute("value", this.k["4"]);
    (updateClass.getElementsByClassName("k_4")[0]).setAttribute("value", this.k["_4"]);
    ((updateClass).getElementsByClassName("k5")[0]).setAttribute("value", this.k["5"]);
    ((updateClass).getElementsByClassName("k_5")[0]).setAttribute("value", this.k["_5"]);

    ((document.getElementsByClassName(paramGroup)[0]).getElementsByClassName("dt")[0]).setAttribute("value", this.dt);

    (updateClass.getElementsByClassName("width")[0]).setAttribute("value", this.size["columns"]);    
    (updateClass.getElementsByClassName("height")[0]).setAttribute("value", this.size["rows"]);
    var c = document.getElementById(sqkMgthis.plotDOMid);
    c.height = sqkMgthis.size["rows"];
    c.width = sqkMgthis.size["columns"];
  };
  this.updateParamsInView("loadedGparams");
  this.updateParamsInView("newGparams");
  
  this.showStatus = function(text, bgColor) {
    var status = document.getElementsByClassName(sqkMgthis.statusDOMclass);
    status[0].style.backgroundColor = bgColor;
    status[0].innerHTML = text;
  };
};
window['EuSol'] = EuSol;


