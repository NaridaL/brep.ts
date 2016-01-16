"use strict";
var log = Math.log, ceil = Math.ceil;
function setContextFont() {
}
window.onload = function () {
	var ctx = $("myCanvas").getContext('2d');

	var str = "400mm";
	var font = "64px monospace";
	
	ctx.font = font ;
	var canvasWidth = 1 << ceil(log(ctx.measureText(str).width)/log(2));
	$("myCanvas").width = canvasWidth;
	
	
	ctx.font = font;
	ctx.fillStyle = "#000000";
	ctx.textAlign = "left";	// This determines the alignment of text, e.g. left, center, right
	ctx.textBaseline = "top";	// This determines the baseline of the text, e.g. top, middle, bottom
	
	ctx.fillText(str, 0, 0);

	
}