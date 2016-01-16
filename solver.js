window.onload = function () {
	var f = math.parse("a + b == 2");
	var g = math.parse("a* b == 4");
	console.log(f.toString(), g);
}