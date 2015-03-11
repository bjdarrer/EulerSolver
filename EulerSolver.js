/**
 * Simulating and visualizing reaction-diffusion systems
 * Programmed by David Jonsson, Athens January 2014
  *                             Stockholm March, April, November 2014, 2015
 * @constructor
 */
function EulerSolver() {
  var sqkMgthis = this;

  /**
   * Size of simulation and canvas
   * @type {number}
   */
  this.size = {
    'rows' : 500,
    'columns' : 500
  }; 
  
  /**
   * Reaction parameters.
   * @type {number}
   */
  var k = {
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
   * @type {number}
   */
  var diffusion = {
    "DG" : 1,
    "DX" : 1,
    "DY" : 1
  };
  
  /**
   * Initial static concentrations.
   * @type {number}
   */
  var concentrations = {
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
   * @type {number}
   */
  this.dt = 0.5;

  /**
   * User changeable settings.
   */
  this.plotDOMid = "c"; // Canvas id
  this.statusDOMclass = "simulationStatus";
  this.runSimulation = true;
  this.mask = "none";
  this.maskParam = 1;
  this.calcStepIterationDOM = document.getElementsByClassName('calcStepIteration')[0];
  this.calcStepDurationDOM = document.getElementsByClassName('calcStepDuration')[0];

  /**
   * Initial static concentrations. Format GXYGXY...
   */
  this.initializeArrays = function() {
    //Initial dynamic concentrations. Format GXYGXY ... 64 bit double array might be needed
    this.GXYarr = new Array(2);
    if (window.Float64Array) {
      this.GXYarr = [new Float64Array(3 * this.size.rows * this.size.columns),
          new Float64Array(3 * this.size.rows * this.size.columns)];
    } else {
      this.GXYarr = [new Array(3 * this.size.rows * this.size.columns),
          new Array(3 * this.size.rows * this.size.columns)];
    };
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
    
  /**
   * Laplacian numerical operator nabla^2, wheight 1 on adjacents and ½ on diagonals.
   */
  var laplace2dDiagonal = function(V, index) {
    return -6*V[index] +
            V[index + 3] + V[index - 3] + V[index + 3 * sqkMgthis.size["columns"]] + V[index - 3 * sqkMgthis.size["columns"]] + //closest neighbours
            0.5*(V[index + 3 * sqkMgthis.size.columns + 3] + V[index + 3 * sqkMgthis.size.columns - 3] + V[index - 3 * sqkMgthis.size.columns + 3] + V[index - 3 * sqkMgthis.size.columns - 3]); //diagonal neighbours
  };

  /**
   * Laplacian numerical operator nabla^2, wheight 1 on adjacents and 0 on diagonals.
   */
  var laplace2dAdjacentOnly = function(V, index) {
    return -4*V[index] +
            V[index + 3] + V[index - 3] + V[index + 3 * sqkMgthis.size["columns"]] + V[index - 3 * sqkMgthis.size["columns"]];  //closest neighbours
  };
//  this.laplace = laplace2dAdjacentOnly;
  this.laplace = laplace2dDiagonal;

  /**
   * Masking, boundary conditions.
   */
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

  /**
   * Name and TeX code for equations.
   * @type {Array}
   */
  this.equations = [
    {name: "model_G_concentrations", tex: '\begin{align}'+"\n" +
                                       '\\ \frac{\partial G}{\partial t} & = k_1A -k_2G + k_{-2}X + D_g \nabla^2 G'+"\n" +
                                       '\\ \frac{\partial X}{\partial t} & = k_2G +k_4X^2Y - k_{3}BX - k_5X + D_x \nabla^2 X'+"\n" +
                                       '\\ \frac{\partial Y}{\partial t} & = k_3BX -k_4X^2Y + D_y \nabla^2 Y'+"\n" +
                                       '\end{align}'},
    {name: "model_G_nondimensionalized", tex: '\begin{align}'+"\n" +
                                       '\\ \frac{\partial G}{\partial t} & = \nabla^2 G -qG + gX +a'+"\n" +
                                       '\\ \frac{\partial X}{\partial t} & = d_x \nabla^2 X + pG - (1+b)X +uY -sX^3 +x^2Y + w'+"\n" +
                                       '\\ \frac{\partial Y}{\partial t} & = d_y \nabla^2 Y +bX -uY +sX^3 -X^2Y'+"\n" +
                                       '\end{align}'}
  ]

  /**
   * iteration step function for concentrations
   * @param {Array} arr3 - Double buffer simulation data.
   * @param {number} readIndex - Index of the double buffers.
   * @param {string} readBuffer - The read buffer of the two buffers.
   * @param {string} writeBuffer - The write buffer of the two buffers.
   */
  this.gxyConcentrations = function(arr3, readIndex, readBuffer, writeBuffer) {
    var G = readBuffer[3 * this.size.columns * i + j];
    var X = readBuffer[3 * this.size.columns * i + j + 1];
    var Y = readBuffer[3 * this.size.columns * i + j + 2];
    G += (diffusion['DG'] * this.laplace(arr3[readIndex], 3 * this.size.columns * i + j) - (k["_1"] + k["2"]) * G + k["_2"] * X + k["1"] * concentrations['A']) * sqkMgthis.dt;
    X += (diffusion['DX']*this.laplace(arr3[readIndex], 3*this.size.columns*i + j + 1) + k["2"]*G - (k["_2"] + k["3"]*concentrations['B'] + k["5"])*X +
        k["_3"]*concentrations['Z']*Y - k["_4"]*X*X*X + k["4"] * X * X * Y + k["_5"] * concentrations['Omega'])*sqkMgthis.dt;
    Y += (diffusion['DY']*this.laplace(arr3[readIndex], 3*this.size.columns*i + j + 2) + k["3"]*concentrations['B']*X - k["_3"]*concentrations['Z']*Y +
        k["_4"]*X*X*X - k["4"]*X*X*Y )*sqkMgthis.dt;
    return [G, X, Y];
  }  

  /**
   *  iteration step function for nondimensionalized potentials
   *  ∂G
   *  ∂t = âˆ‡2G âˆ’ qG + gX + a (7a)
   *  ∂X
   *  ∂t = dxâˆ‡  2X + pG âˆ’ (1 + b)X + uY âˆ’ sX^3 + X^2Y + w (7b)
   *  ∂Y
   *  ∂t = dyâˆ‡ 2Y + bX âˆ’ uY + sX3 âˆ’ X^2Y (7c)
   *
   * @param {Array} arr3 - Double buffer simulation data.
   * @param {number} readIndex - Index of the double buffers.
   */
  this.gxyPotentials = function(arr3, readIndex) {
    var G = arr3[readIndex][3 * this.size.columns * i + j];
    var X = arr3[readIndex][3 * this.size.columns * i + j + 1];
    var Y = arr3[readIndex][3 * this.size.columns * i + j + 2];
    G += (this.laplace(arr3[readIndex], 3 * this.size.columns * i + j) - (k["_1"] + k["2"]) / (k["_2"] + k["5"]) * G + 1/(1-k["5"]/k["_2"]) * X + k["1"] * Math.sqrt(k["4"])*Math.pow(k["_2"] +k["5"], 1.5) * concentrations['A'] ) * dt;
    X += (diffusion['DX'] / diffusion['DG'] * this.laplace(arr3[readIndex], 3 * this.size.columns * i + j + 1) + k["2"]/(k["_2"] - k["5"]) * G -
        + (1 + k["_2"] / (k["_2"] + k["5"]) ) * concentrations['B'] * X + //Done until here
        k["_3"] * concentrations['Z'] * Y - k["_4"] * X * X * X + k["4"] * X * X * Y + k["_5"] * concentrations['Î©']) * dt;
    Y += (diffusion['DY'] * this.laplace(arr3[readIndex], 3 * this.size.columns * i + j + 2) + k["3"] * concentrations['B'] * X - k["_3"]*concentrations['Z']*Y +
        k["_4"] * X * X * X - k["4"]*X*X * Y ) * dt;
    return [G, X, Y];
  }  

 /**
  * Choosing the set of differential equations for returning 3 valued GXY in the iteration.
  * @type {Object}
  */
  this.gxy=this.gxyConcentrations;

  /**
   *  Reaction diffusion functions. Iterating, special treating the boundary.
   *
   * @param {Array} arr3 - Double buffer simulation data.
   */
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
          conc3 = sqkMgthis.gxy(arr3, readIndex, readBuffer);

          writeBuffer[3 * this.size.columns * i + j]     = conc3[0];
          writeBuffer[3 * this.size.columns * i + j + 1] = conc3[1];
          writeBuffer[3 * this.size.columns * i + j + 2] = conc3[2];

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
          }
        }
      };
    };    
    this.iteration++;
    this.calcStepDuration = (performance.now() - functionStart).toFixed(3);
  };
  
  /**
   *  Plot the data in a 2D Canvas, (G, X, Y) mapped to (R, G, B).
   *
   * @param {Array} arr3 - Double buffer simulation data.
   */
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
  
  /**
   *  Loop the simulation with setTimeout(...,0).
   */
  this.loopSimulation = function () {
    if (sqkMgthis.runSimulation) {
      sqkMgthis.gxyStep(sqkMgthis.GXYarr);
      setTimeout(sqkMgthis.loopSimulation, 0);
    }
  };
  
  /**
   *  Loop the visualization with requestAnimationFrame or new data.
   */
  this.loopVisualization = function () {
    if (sqkMgthis.shownIteration < sqkMgthis.iteration) {
      sqkMgthis.plot2D(sqkMgthis.GXYarr, sqkMgthis.plotDOMid);
      sqkMgthis.shownIteration = sqkMgthis.iteration;
    }
    sqkMgthis.loopRequestAnimationFrameID = requestAnimationFrame(sqkMgthis.loopVisualization);
  };

  /**
   *  Fill variables and objects from the view's "changeable values".
   */
  this.loadNewParamsFromView = function() {
    var newParamsClass = document.getElementsByClassName("newGparams")[0];
    concentrations.A = parseFloat((newParamsClass.getElementsByClassName("A")[0]).value);
    concentrations.B = parseFloat((newParamsClass.getElementsByClassName("B")[0]).value);
    concentrations.G = parseFloat((newParamsClass.getElementsByClassName("G")[0]).value);
    concentrations.X = parseFloat((newParamsClass.getElementsByClassName("X")[0]).value);
    concentrations.Y = parseFloat((newParamsClass.getElementsByClassName("Y")[0]).value);
    concentrations.Z = parseFloat((newParamsClass.getElementsByClassName("Z")[0]).value);
    concentrations.Omega = parseFloat((newParamsClass.getElementsByClassName("Omega")[0]).value);

    diffusion.DG = parseFloat((newParamsClass.getElementsByClassName("DG")[0]).value);
    diffusion.DX = parseFloat((newParamsClass.getElementsByClassName("DX")[0]).value);
    diffusion.DY = parseFloat((newParamsClass.getElementsByClassName("DY")[0]).value);

    k["1"]  = parseFloat((newParamsClass.getElementsByClassName("k1")[0]).value);
    k["_1"] = parseFloat((newParamsClass.getElementsByClassName("k_1")[0]).value);
    k["2"]  = parseFloat((newParamsClass.getElementsByClassName("k2")[0]).value);
    k["_2"] = parseFloat((newParamsClass.getElementsByClassName("k_2")[0]).value);
    k["3"]  = parseFloat((newParamsClass.getElementsByClassName("k3")[0]).value);
    k["_3"] = parseFloat((newParamsClass.getElementsByClassName("k_3")[0]).value);
    k["4"]  = parseFloat((newParamsClass.getElementsByClassName("k4")[0]).value);
    k["_4"] = parseFloat((newParamsClass.getElementsByClassName("k_4")[0]).value);
    k["5"]  = parseFloat((newParamsClass.getElementsByClassName("k5")[0]).value);
    k["_5"] = parseFloat((newParamsClass.getElementsByClassName("k_5")[0]).value);

    this.dt = parseFloat((newParamsClass.getElementsByClassName("dt")[0]).value);

    this.size.rows = parseFloat((newParamsClass.getElementsByClassName("height")[0]).value);
    this.size.columns = parseFloat((newParamsClass.getElementsByClassName("width")[0]).value);

    this.updateParamsInView("loadedGparams");
  };

  /**
   *  Fill variables and objects from the view's "changeable values".
   *
   * @param {String} paramGroup - class of parameters in view to update.
   */
  this.updateParamsInView = function(paramGroup) {    
    var updateClass = document.getElementsByClassName(paramGroup)[0];
    (updateClass.getElementsByClassName("A")[0]).setAttribute("value", concentrations.A);
    (updateClass.getElementsByClassName("B")[0]).setAttribute("value", concentrations.B);
    (updateClass.getElementsByClassName("G")[0]).setAttribute("value", concentrations.G);
    (updateClass.getElementsByClassName("X")[0]).setAttribute("value", concentrations.X);
    (updateClass.getElementsByClassName("Y")[0]).setAttribute("value", concentrations.Y);
    (updateClass.getElementsByClassName("Z")[0]).setAttribute("value", concentrations.Z);
    (updateClass.getElementsByClassName("Omega")[0]).setAttribute("value", concentrations.Omega);
    
    (updateClass.getElementsByClassName("DG")[0]).setAttribute("value", diffusion.DG);
    (updateClass.getElementsByClassName("DX")[0]).setAttribute("value", diffusion.DX);
    (updateClass.getElementsByClassName("DY")[0]).setAttribute("value", diffusion.DY);

    (updateClass.getElementsByClassName("k1")[0]).setAttribute("value", k["1"]);
    (updateClass.getElementsByClassName("k_1")[0]).setAttribute("value", k["_1"]);
    (updateClass.getElementsByClassName("k2")[0]).setAttribute("value", k["2"]);
    (updateClass.getElementsByClassName("k_2")[0]).setAttribute("value", k["_2"]);
    (updateClass.getElementsByClassName("k3")[0]).setAttribute("value", k["3"]);
    (updateClass.getElementsByClassName("k_3")[0]).setAttribute("value", k["_3"]);
    (updateClass.getElementsByClassName("k4")[0]).setAttribute("value", k["4"]);
    (updateClass.getElementsByClassName("k_4")[0]).setAttribute("value", k["_4"]);
    ((updateClass).getElementsByClassName("k5")[0]).setAttribute("value", k["5"]);
    ((updateClass).getElementsByClassName("k_5")[0]).setAttribute("value", k["_5"]);

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

