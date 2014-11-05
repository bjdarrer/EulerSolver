/* 
 * Making stuff more compatible between browsers
 */

window.requestAnimationFrame = window.requestAnimationFrame || window.mozRequestAnimationFrame ||
  window.webkitRequestAnimationFrame || window.msRequestAnimationFrame ||
  function( callback ){
    window.setTimeout(callback, 1000 / 60);
  };

//Fix for Safari and Android Stock  
if (!window.performance) {
  window.performance = {};
  window.performance.now = Date.now; 
}
