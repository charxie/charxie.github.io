//image layer manipulation

var middleX, middleY, pos;

function start1() {
  document.layers["imgLayer"].document.map.src="pix/contact_map.gif";
  start();
}

function start2() {
  document.layers["imgLayer"].document.map.src="pix/gp.jpeg";
  start();
}

function start() {
  // get size of image
  //var width= document.layers["imgLayer"].document.map.width;
  //var height= document.layers["imgLayer"].document.map.height;

  var width = document.width;
  var height = 200;
  
  // calculate pixel in the middle of image
  middleX= Math.round(width/2);
  middleY= Math.round(height/2);

  // starting position
  pos= 0;

  // start it!
  show();
}

function show() {

  // increase size of clipping area  
  pos+= 2; // step size
  document.layers["imgLayer"].clip.left= middleX- pos;
  document.layers["imgLayer"].clip.top= middleY- pos;
  document.layers["imgLayer"].clip.right= middleX+ pos;
  document.layers["imgLayer"].clip.bottom= middleY+ pos;

  // check if the whole image has been displayed
  if (!((pos > middleX) && (pos > middleY))) 
    setTimeout("show()", 20);  

}

