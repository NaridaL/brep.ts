/*
 * lightgl.js
 * http://github.com/evanw/lightgl.js/
 *
 * Copyright 2011 Evan Wallace
 * Released under the MIT license
 */
"use strict"
var GL = (function() {

// src/main.js
// The internal `gl` variable holds the current WebGL context.
	var gl;
	var V3 = NLA.Vector3, M4 = NLA.Matrix4x4

	var GL = {
		// ### Initialization
		//
		// `GL.create()` creates a new WebGL context and augments it with more
		// methods. The alpha channel is disabled by default because it usually causes
		// unintended transparencies in the canvas.
		/**
		 *
		 * @param options
		 * @returns {WebGLRenderingContext}
		 */
		create: function(options) {
			options = options || {};
			var canvas = options.canvas || document.createElement('canvas');
			if (!options.canvas) {
				canvas.width = 800;
				canvas.height = 600;
			}
			if (!('alpha' in options)) options.alpha = false;
			try { gl = canvas.getContext('webgl', options); } catch (e) {}
			try { gl = gl || canvas.getContext('experimental-webgl', options); } catch (e) {}
			if (!gl) throw new Error('WebGL not supported');
			gl.HALF_FLOAT_OES = 0x8D61;
			addMatrixStack();
			addImmediateMode();
			addEventListeners();
			addOtherMethods();
			return gl;
		},

		// `GL.keys` contains a mapping of key codes to booleans indicating whether
		// that key is currently pressed.
		keys: {},

		// Export all external classes.
		Indexer: Indexer,
		Buffer: Buffer,
			Mesh: Mesh,
		HitTest: HitTest,
		Raytracer: Raytracer,
		Shader: Shader,
		Texture: Texture,
	};

// ### Matrix stack
//
// Implement the OpenGL modelview and projection matrix stacks, along with some
// other useful GLU matrix functions.

	function addMatrixStack() {
		gl.MODELVIEW = ENUM | 1;
		gl.PROJECTION = ENUM | 2;
		var tempMatrix = M4();
		var resultMatrix = M4();
		gl.modelviewMatrix = M4();
		gl.projectionMatrix = M4();
		var modelviewStack = [];
		var projectionStack = [];
		var matrix, stack;
		gl.matrixMode = function(mode) {
			switch (mode) {
				case gl.MODELVIEW:
					matrix = 'modelviewMatrix';
					stack = modelviewStack;
					break;
				case gl.PROJECTION:
					matrix = 'projectionMatrix';
					stack = projectionStack;
					break;
				default:
					throw new Error('invalid matrix mode ' + mode);
			}
		};
		gl.loadIdentity = function() {
			M4.identity(gl[matrix]);
		};
		gl.loadMatrix = function(m4) {
			M4.copy(m4, gl[matrix])
		}
		gl.multMatrix = function(m) {
			M4.multiply(gl[matrix], m, resultMatrix);
			var temp = resultMatrix
			resultMatrix = gl[matrix]
			gl[matrix] = temp
		};
		gl.perspective = function(fov, aspect, near, far) {
			gl.multMatrix(M4.perspective(fov, aspect, near, far, tempMatrix));
		};
		gl.frustum = function(l, r, b, t, n, f) {
			gl.multMatrix(M4.frustum(l, r, b, t, n, f, tempMatrix));
		};
		gl.ortho = function(l, r, b, t, n, f) {
			gl.multMatrix(M4.ortho(l, r, b, t, n, f, tempMatrix));
		};
		gl.scale = function(x, y, z) {
			gl.multMatrix(M4.scaling(x, y, z, tempMatrix));
		};
		gl.translate = function(x, y, z) {
			if (undefined !== y) {
				gl.multMatrix(M4.translation(x, y, z, tempMatrix));
			} else {
				gl.multMatrix(M4.translation(x, tempMatrix));
			}
		};
		gl.rotate = function(a, x, y, z) {
			gl.multMatrix(M4.rotation(a, x, y, z, tempMatrix));
		};
		gl.lookAt = function(eye, center, up) {
			gl.multMatrix(M4.lookAt(eye, center, up, tempMatrix));
		};
		/**
		 * Test
		 */
		gl.pushMatrix = function() {
			stack.push(M4.copy(gl[matrix]));
		};
		gl.popMatrix = function() {
			gl[matrix] = stack.pop()
		};
		gl.project = function(objX, objY, objZ, modelview, projection, viewport) {
			modelview = modelview || gl.modelviewMatrix;
			projection = projection || gl.projectionMatrix;
			viewport = viewport || gl.getParameter(gl.VIEWPORT);
			var point = projection.transformPoint(modelview.transformPoint(V3(objX, objY, objZ)));
			return V3(
				viewport[0] + viewport[2] * (point.x * 0.5 + 0.5),
				viewport[1] + viewport[3] * (point.y * 0.5 + 0.5),
				point.z * 0.5 + 0.5
			);
		};
		gl.unProject = function(winX, winY, winZ, modelview, projection, viewport) {
			modelview = modelview || gl.modelviewMatrix;
			projection = projection || gl.projectionMatrix;
			viewport = viewport || gl.getParameter(gl.VIEWPORT);
			var point = V3(
				(winX - viewport[0]) / viewport[2] * 2 - 1,
				(winY - viewport[1]) / viewport[3] * 2 - 1,
				winZ * 2 - 1
			);
			return M4.inverse(M4.multiply(projection, modelview, tempMatrix), resultMatrix).transformPoint(point);
		};
		gl.matrixMode(gl.MODELVIEW);
	}

// ### Immediate mode
//
// Provide an implementation of OpenGL's deprecated immediate mode. This is
// depricated for a reason: constantly re-specifying the geometry is a bad
// idea for performance. You should use a `GL.Mesh` instead, which specifies
// the geometry once and caches it on the graphics card. Still, nothing
// beats a quick `gl.begin(gl.POINTS); gl.vertex(1, 2, 3); gl.end();` for
// debugging. This intentionally doesn't implement fixed-function lighting
// because it's only meant for quick debugging tasks.

	function addImmediateMode() {
		var immediateMode = {
			mesh: new Mesh({ coords: true, colors: true, triangles: false }),
			mode: -1,
			coord: [0, 0, 0, 0],
			color: [1, 1, 1, 1],
			pointSize: 1,
			shader: new Shader(`
uniform float pointSize;
varying vec4 color;
varying vec4 coord;
void main() {
	color = gl_Color;
	coord = gl_TexCoord;
	gl_Position = gl_ModelViewProjectionMatrix * gl_Vertex;
	gl_PointSize = pointSize;
}`, `
uniform sampler2D texture;
uniform float pointSize;
uniform bool useTexture;
varying vec4 color;
varying vec4 coord;
void main() {
	gl_FragColor = color;
	if (useTexture) gl_FragColor *= texture2D(texture, coord.xy);
}`)
		};
		gl.pointSize = function(pointSize) {
			immediateMode.shader.uniforms({ pointSize: pointSize });
		};
		gl.begin = function(mode) {
			if (immediateMode.mode != -1) throw new Error('mismatched gl.begin() and gl.end() calls');
			immediateMode.mode = mode;
			immediateMode.mesh.colors = [];
			immediateMode.mesh.coords = [];
			immediateMode.mesh.vertices = [];
		};
		gl.color = function(r, g, b, a) {
			immediateMode.color = (arguments.length == 1) ? r.toArray().concat(1) : [r, g, b, a || 1];
		};
		gl.texCoord = function(s, t) {
			immediateMode.coord = (arguments.length == 1) ? s.toArray(2) : [s, t];
		};
		gl.vertex = function(x, y, z) {
			immediateMode.mesh.colors.push(immediateMode.color);
			immediateMode.mesh.coords.push(immediateMode.coord);
			immediateMode.mesh.vertices.push(arguments.length == 1 ? x.toArray() : [x, y, z]);
		};
		gl.end = function() {
			if (immediateMode.mode == -1) throw new Error('mismatched gl.begin() and gl.end() calls');
			console.log(immediateMode.mesh.toSource())
			immediateMode.mesh.compile();
			immediateMode.shader.uniforms({
				useTexture: !!gl.getParameter(gl.TEXTURE_BINDING_2D)
			}).draw(immediateMode.mesh, immediateMode.mode);
			immediateMode.mode = -1;
		};
	}

// ### Improved mouse events
//
// This adds event listeners on the `gl.canvas` element that call
// `gl.onmousedown()`, `gl.onmousemove()`, and `gl.onmouseup()` with an
// augmented event object. The event object also has the properties `x`, `y`,
// `deltaX`, `deltaY`, and `dragging`.
	function addEventListeners() {
		var context = gl, oldX = 0, oldY = 0, buttons = {}, hasOld = false;
		function isDragging() {
			for (var b in buttons) {
				if (b in buttons && buttons[b]) return true;
			}
			return false;
		}
		function augment(original) {
			// Make a copy of original, a native `MouseEvent`, so we can overwrite
			// WebKit's non-standard read-only `x` and `y` properties (which are just
			// duplicates of `pageX` and `pageY`). We can't just use
			// `Object.create(original)` because some `MouseEvent` functions must be
			// called in the context of the original event object.
			var e = {};
			for (var name in original) {
				if (typeof original[name] == 'function') {
					e[name] = (function(callback) {
						return function() {
							callback.apply(original, arguments);
						};
					})(original[name]);
				} else {
					e[name] = original[name];
				}
			}
			e.original = original;
			e.x = e.pageX;
			e.y = e.pageY;
			for (var obj = gl.canvas; obj; obj = obj.offsetParent) {
				e.x -= obj.offsetLeft;
				e.y -= obj.offsetTop;
			}
			if (hasOld) {
				e.deltaX = e.x - oldX;
				e.deltaY = e.y - oldY;
			} else {
				e.deltaX = 0;
				e.deltaY = 0;
				hasOld = true;
			}
			oldX = e.x;
			oldY = e.y;
			e.dragging = isDragging();
			e.preventDefault = function (){}//e.original.preventDefault
			e.stopPropagation = e.original.stopPropagation
			return e;
		}
		function mousedown(e) {
			gl = context;
			if (!isDragging()) {
				// Expand the event handlers to the document to handle dragging off canvas.
				on(document, 'mousemove', mousemove);
				on(document, 'mouseup', mouseup);
				off(gl.canvas, 'mousemove', mousemove);
				off(gl.canvas, 'mouseup', mouseup);
			}
			buttons[e.which] = true;
			e = augment(e);
			if (gl.onmousedown) gl.onmousedown(e);
			e.preventDefault();
		}
		function mousemove(e) {
			gl = context;
			e = augment(e);
			if (gl.onmousemove) gl.onmousemove(e);
			e.preventDefault();
		}
		function mouseup(e) {
			gl = context;
			buttons[e.which] = false;
			if (!isDragging()) {
				// Shrink the event handlers back to the canvas when dragging ends.
				off(document, 'mousemove', mousemove);
				off(document, 'mouseup', mouseup);
				on(gl.canvas, 'mousemove', mousemove);
				on(gl.canvas, 'mouseup', mouseup);
			}
			e = augment(e);
			if (gl.onmouseup) gl.onmouseup(e);
			e.preventDefault();
		}
		function reset() {
			hasOld = false;
		}
		function resetAll() {
			buttons = {};
			hasOld = false;
		}
		on(gl.canvas, 'mousedown', mousedown);
		on(gl.canvas, 'mousemove', mousemove);
		on(gl.canvas, 'mouseup', mouseup);
		on(gl.canvas, 'mouseover', reset);
		on(gl.canvas, 'mouseout', reset);
		on(document, 'contextmenu', resetAll);
	}

// ### Automatic keyboard state
//
// The current keyboard state is stored in `GL.keys`, a map of integer key
// codes to booleans indicating whether that key is currently pressed. Certain
// keys also have named identifiers that can be used directly, such as
// `GL.keys.SPACE`. Values in `GL.keys` are initially undefined until that
// key is pressed for the first time. If you need a boolean value, you can
// cast the value to boolean by applying the not operator twice (as in
// `!!GL.keys.SPACE`).

	function mapKeyCode(code) {
		var named = {
			8: 'BACKSPACE',
			9: 'TAB',
			13: 'ENTER',
			16: 'SHIFT',
			27: 'ESCAPE',
			32: 'SPACE',
			37: 'LEFT',
			38: 'UP',
			39: 'RIGHT',
			40: 'DOWN'
		};
		return named[code] || (code >= 65 && code <= 90 ? String.fromCharCode(code) : null);
	}

	function on(element, name, callback) {
		element.addEventListener(name, callback);
	}

	function off(element, name, callback) {
		element.removeEventListener(name, callback);
	}

	on(document, 'keydown', function(e) {
		if (!e.altKey && !e.ctrlKey && !e.metaKey) {
			var key = mapKeyCode(e.keyCode);
			if (key) GL.keys[key] = true;
			GL.keys[e.keyCode] = true;
		}
	});

	on(document, 'keyup', function(e) {
		if (!e.altKey && !e.ctrlKey && !e.metaKey) {
			var key = mapKeyCode(e.keyCode);
			if (key) GL.keys[key] = false;
			GL.keys[e.keyCode] = false;
		}
	});

	function addOtherMethods() {
		// ### Multiple contexts
		//
		// When using multiple contexts in one web page, `gl.makeCurrent()` must be
		// called before issuing commands to a different context.
		(function(context) {
			gl.makeCurrent = function() {
				gl = context;
			};
		})(gl);

		// ### Animation
		//
		// Call `gl.animate()` to provide an animation loop that repeatedly calls
		// `gl.onupdate()` and `gl.ondraw()`.
		gl.animate = function() {
			var post =
				window.requestAnimationFrame ||
				window.mozRequestAnimationFrame ||
				window.webkitRequestAnimationFrame ||
				function(callback) { setTimeout(callback, 1000 / 60); };
			var time = new Date().getTime();
			var context = gl;
			function update() {
				gl = context;
				var now = new Date().getTime();
				if (gl.onupdate) gl.onupdate((now - time) / 1000);
				if (gl.ondraw) gl.ondraw();
				post(update);
				time = now;
			}
			update();
		};

		// ### Fullscreen
		//
		// Provide an easy way to get a fullscreen app running, including an
		// automatic 3D perspective projection matrix by default. This should be
		// called once.
		//
		// Just fullscreen, no automatic camera:
		//
		//     gl.fullscreen({ camera: false });
		//
		// Adjusting field of view, near plane distance, and far plane distance:
		//
		//     gl.fullscreen({ fov: 45, near: 0.1, far: 1000 });
		//
		// Adding padding from the edge of the window:
		//
		//     gl.fullscreen({ paddingLeft: 250, paddingBottom: 60 });
		//
		gl.fullscreen = function(options) {
			options = options || {};
			var top = options.paddingTop || 0;
			var left = options.paddingLeft || 0;
			var right = options.paddingRight || 0;
			var bottom = options.paddingBottom || 0;
			if (!document.body) {
				throw new Error('document.body doesn\'t exist yet (call gl.fullscreen() from ' +
					'window.onload() or from inside the <body> tag)');
			}
			document.body.appendChild(gl.canvas);
			document.body.style.overflow = 'hidden';
			gl.canvas.style.position = 'absolute';
			gl.canvas.style.left = left + 'px';
			gl.canvas.style.top = top + 'px';
			function resize() {
				gl.canvas.width = window.innerWidth - left - right;
				gl.canvas.height = window.innerHeight - top - bottom;
				gl.viewport(0, 0, gl.canvas.width, gl.canvas.height);
				if (options.camera || !('camera' in options)) {
					gl.matrixMode(gl.PROJECTION);
					gl.loadIdentity();
					gl.perspective(options.fov || 45, gl.canvas.width / gl.canvas.height,
						options.near || 0.1, options.far || 1000);
					gl.matrixMode(gl.MODELVIEW);
				}
				if (gl.ondraw) gl.ondraw();
			}
			on(window, 'resize', resize);
			resize();
		};
	}

// A value to bitwise-or with new enums to make them distinguishable from the
// standard WebGL enums.
	var ENUM = 0x12340000;

// src/matrix.js
// Represents a 4x4 matrix stored in row-major order that uses Float32Arrays
// when available. Matrix operations can either be done using convenient
// methods that return a new matrix for the result or optimized methods
// that store the result in an existing matrix to avoid generating garbage.

	var hasFloat32Array = (typeof Float32Array != 'undefined');

// ### new GL.Matrix([elements])
//
// src/mesh.js
// Represents indexed triangle geometry with arbitrary additional attributes.
// You need a shader to draw a mesh; meshes can't draw themselves.
//
// A mesh is a collection of `GL.Buffer` objects which are either vertex buffers
// (holding per-vertex attributes) or index buffers (holding the order in which
// vertices are rendered). By default, a mesh has a position vertex buffer called
// `vertices` and a triangle index buffer called `triangles`. New buffers can be
// added using `addVertexBuffer()` and `addIndexBuffer()`. Two strings are
// required when adding a new vertex buffer, the name of the data array on the
// mesh instance and the name of the GLSL attribute in the vertex shader.
//
// Example usage:
//
//     var mesh = new GL.Mesh({ coords: true, lines: true });
//
//     // Default attribute "vertices", available as "gl_Vertex" in
//     // the vertex shader
//     mesh.vertices = [[0, 0, 0], [1, 0, 0], [0, 1, 0], [1, 1, 0]];
//
//     // Optional attribute "coords" enabled in constructor,
//     // available as "gl_TexCoord" in the vertex shader
//     mesh.coords = [[0, 0], [1, 0], [0, 1], [1, 1]];
//
//     // Custom attribute "weights", available as "weight" in the
//     // vertex shader
//     mesh.addVertexBuffer('weights', 'weight');
//     mesh.weights = [1, 0, 0, 1];
//
//     // Default index buffer "triangles"
//     mesh.triangles = [[0, 1, 2], [2, 1, 3]];
//
//     // Optional index buffer "lines" enabled in constructor
//     mesh.lines = [[0, 1], [0, 2], [1, 3], [2, 3]];
//
//     // Upload provided data to GPU memory
//     mesh.compile();

// ### new GL.Indexer()
//
// Generates indices into a list of unique objects from a stream of objects
// that may contain duplicates. This is useful for generating compact indexed
// meshes from unindexed data.
	function Indexer() {
		this.unique = [];
		this.indices = [];
		this.map = {};
	}

	Indexer.prototype = {
		// ### .add(v)
		//
		// Adds the object `obj` to `unique` if it hasn't already been added. Returns
		// the index of `obj` in `unique`.
		add: function(obj) {
			var key = JSON.stringify(obj);
			if (!(key in this.map)) {
				this.map[key] = this.unique.length;
				this.unique.push(obj);
			}
			return this.map[key];
		}
	};

// ### new GL.Buffer(target, type)
//
// Provides a simple method of uploading data to a GPU buffer. Example usage:
//
//     var vertices = new GL.Buffer(gl.ARRAY_BUFFER, Float32Array);
//     var indices = new GL.Buffer(gl.ELEMENT_ARRAY_BUFFER, Uint16Array);
//     vertices.data = [[0, 0, 0], [1, 0, 0], [0, 1, 0], [1, 1, 0]];
//     indices.data = [[0, 1, 2], [2, 1, 3]];
//     vertices.compile();
//     indices.compile();
//
	function Buffer(target, type) {
		this.buffer = null;
		this.target = target;
		this.type = type;
		this.data = [];
	}

	Buffer.prototype = {
		// ### .compile(type)
		//
		// Upload the contents of `data` to the GPU in preparation for rendering. The
		// data must be a list of lists where each inner list has the same length. For
		// example, each element of data for vertex normals would be a list of length three.
		// This will remember the data length and element length for later use by shaders.
		// The type can be either `gl.STATIC_DRAW` or `gl.DYNAMIC_DRAW`, and defaults to
		// `gl.STATIC_DRAW`.
		//
		// This could have used `[].concat.apply([], this.data)` to flatten
		// the array but Google Chrome has a maximum number of arguments so the
		// concatenations are chunked to avoid that limit.
		compile: function(type) {
			this.buffer = this.buffer || gl.createBuffer();
			if (this.data.length == 0) {
				console.warn("empty buffer " + this.name)
				//console.trace()
			}
			if (this.data.length == 0 || this.data[0] instanceof V3) {
				var buffer = V3.flattenV3Array(this.data) // asserts that all elements are V3s
				this.buffer.length = this.data.length
				this.buffer.spacing = 3
				this.buffer.count = this.data.length
			} else {
				var data = [];
				for (var i = 0, chunk = 10000; i < this.data.length; i += chunk) {
					data = Array.prototype.concat.apply(data, this.data.slice(i, i + chunk));
				}
				var buffer = new this.type(data)
				var spacing = this.data.length ? data.length / this.data.length : 0;
				assert(spacing % 1 == 0, `buffer ${this.name}elements not of consistent size, average size is ` + spacing)
				assert(data.every(v => 'number' == typeof v), () => "data.every(v => 'number' == typeof v)" + data.toSource())
				this.buffer.length = data.length
				if (NLA.DEBUG) { this.maxValue = data.max() }
				this.buffer.spacing = spacing
				this.buffer.count = data.length / spacing
			}
			gl.bindBuffer(this.target, this.buffer)
			gl.bufferData(this.target, buffer, type || gl.STATIC_DRAW);
		}
	};

// ### new GL.Mesh([options])
//
// Represents a collection of vertex buffers and index buffers. Each vertex
// buffer maps to one attribute in GLSL and has a corresponding property set
// on the Mesh instance. There is one vertex buffer by default: `vertices`,
// which maps to `gl_Vertex`. The `coords`, `normals`, and `colors` vertex
// buffers map to `gl_TexCoord`, `gl_Normal`, and `gl_Color` respectively,
// and can be enabled by setting the corresponding options to true. There are
// two index buffers, `triangles` and `lines`, which are used for rendering
// `gl.TRIANGLES` and `gl.LINES`, respectively. Only `triangles` is enabled by
// default, although `computeWireframe()` will add a normal buffer if it wasn't
// initially enabled.
	/**
	 *
	 * @param options
	 * @constructor
	 */
	function Mesh(options) {
		options = options || {};
		this.hasBeenCompiled = false
		this.vertexBuffers = {}
		this.indexBuffers = {}
		this.addVertexBuffer('vertices', 'gl_Vertex');
		if (options.coords) this.addVertexBuffer('coords', 'gl_TexCoord')
		if (options.normals) this.addVertexBuffer('normals', 'gl_Normal')
		if (options.colors) this.addVertexBuffer('colors', 'gl_Color')
		if (!('triangles' in options) || options.triangles) this.addIndexBuffer('TRIANGLES')
		if (options.lines) this.addIndexBuffer('LINES')
	}

	Mesh.prototype = {
		// ### .addVertexBuffer(name, attribute)
		//
		// Add a new vertex buffer with a list as a property called `name` on this object
		// and map it to the attribute called `attribute` in all shaders that draw this mesh.
		addVertexBuffer: function(name, attribute) {
			assert(!this.vertexBuffers[attribute])
			assert(!this[name])
			this.hasBeenCompiled = false
			assert('string' == typeof name)
			assert('string' == typeof attribute)
			var buffer = this.vertexBuffers[attribute] = new Buffer(gl.ARRAY_BUFFER, Float32Array);
			buffer.name = name;
			this[name] = [];
		},

		// ### .addIndexBuffer(name)
		//
		// Add a new index buffer with a list as a property called `name` on this object.
		addIndexBuffer: function(name) {
			this.hasBeenCompiled = false
			var buffer = this.indexBuffers[name] = new Buffer(gl.ELEMENT_ARRAY_BUFFER, Uint16Array);
			buffer.name = name;
			this[name.toLowerCase()] = [];
		},

		// ### .compile()
		//
		// Upload all attached buffers to the GPU in preparation for rendering. This
		// doesn't need to be called every frame, only needs to be done when the data
		// changes.
		compile: function() {
			// TODO: ensure element buffer indexes are in bounds
			var minVSize = Infinity
			var minBufferName
			for (var attribute in this.vertexBuffers) {
				var buffer = this.vertexBuffers[attribute]
				buffer.data = this[buffer.name]
				buffer.compile()
				if (this[buffer.name].length < minVSize) {
					minBufferName = attribute
					minVSize = this[buffer.name].length
				}
			}

			for (var name in this.indexBuffers) {
				var buffer = this.indexBuffers[name];
				buffer.data = this[name.toLowerCase()];
				buffer.compile();
				if (NLA.DEBUG && buffer.maxValue >= minVSize) {
					throw new Error(`max index value for buffer ${name}
					is too large ${buffer.maxValue} min Vbuffer size: ${minVSize} ${minBufferName}`)
				}
			}
			this.hasBeenCompiled = true
		},

		// ### .transform(matrix)
		//
		// Transform all vertices by `matrix` and all normals by the inverse transpose
		// of `matrix`.
		transform: function(m4) {
			var mesh = new GL.Mesh({vertices: this.vertices, normals: this.normals})
			mesh.vertices = m4.transformedPoints(this.vertices)
			if (this.normals) {
				var invTrans = m4.inverse().transpose();
				this.normals = this.normals.map(n => invTrans.transformVector(n).normalized())
			}
			mesh.triangles = this.triangles
			mesh.compile()
			return mesh
		},

		// ### .computeNormals()
		//
		// Computes a new normal for each vertex from the average normal of the
		// neighboring triangles. This means adjacent triangles must share vertices
		// for the resulting normals to be smooth.
		computeNormals: function() {
			if (!this.normals) this.addVertexBuffer('normals', 'gl_Normal')
			this.normals = Array.fromFunction(this.vertices.length, i => V3.ZERO)

			for (var i = 0; i < this.triangles.length; i++) {
				var t = this.triangles[i]
				var a = this.vertices[t[0]]
				var b = this.vertices[t[1]]
				var c = this.vertices[t[2]]
				var normal = b.minus(a).cross(c.minus(a)).normalized();
				this.normals[t[0]] = this.normals[t[0]].plus(normal);
				this.normals[t[1]] = this.normals[t[1]].plus(normal);
				this.normals[t[2]] = this.normals[t[2]].plus(normal);
			}
			for (var i = 0; i < this.vertices.length; i++) {
				this.normals[i] = this.normals[i].normalized()
			}
			return this
		},

		// ### .computeWireframe()
		//
		// Populate the `lines` index buffer from the `triangles` index buffer.
		computeWireframe: function() {
			assert(false)
			var canonEdges = new Set()
			function canonEdge(i0, i1) {
				var iMin = min(i0, i1), iMax = max(i0, i1)
				return (iMin << 16) | iMax
			}
			function uncanonEdge(key) {
				return [key >> 16, key & 0xffff]
			}
			for (var i = 0; i < this.triangles.length; i++) {
				var t = this.triangles[i]
				canonEdges.add(canonEdge(t[0], t[1]))
				canonEdges.add(canonEdge(t[1], t[2]))
				canonEdges.add(canonEdge(t[2], t[0]))
			}
			if (!this.lines) this.addIndexBuffer('LINES');
			canonEdges.forEach(val => this.lines.push(val >> 16, val & 0xffff))
			return this
		},
		computeWireframeFromFlatTriangles: function() {
			var canonEdges = new Set()
			function canonEdge(i0, i1) {
				var iMin = min(i0, i1), iMax = max(i0, i1)
				return (iMin << 16) | iMax
			}
			function uncanonEdge(key) {
				return [key >> 16, key & 0xffff]
			}
			var t = this.triangles
			for (var i = 0; i < this.triangles.length; i += 3) {
				canonEdges.add(canonEdge(t[i + 0], t[i + 1]))
				canonEdges.add(canonEdge(t[i + 1], t[i + 2]))
				canonEdges.add(canonEdge(t[i + 2], t[i + 0]))
			}
			if (!this.lines) this.addIndexBuffer('LINES');
			//this.lines = new Array(canonEdges.size)
			canonEdges.forEach(val => this.lines.push(val >> 16, val & 0xffff))
			return this
		},
		computeNormalLines: function (length) {
			length = length || 1
			var vs = this.vertices, si = this.vertices.length
			if (!this.lines) this.addIndexBuffer('LINES');

			for (var i = 0; i < this.normals.length; i++) {
				vs[si + i] = vs[i].plus(this.normals[i].times(length))
				this.lines.push(si + i, i)
			}
		},

		// ### .getAABB()
		//
		// Computes the axis-aligned bounding box, which is an object whose `min` and
		// `max` properties contain the minimum and maximum coordinates of all vertices.
		getAABB: function() {
			return new AABB().addPoints(this.vertices)
		},

		// ### .getBoundingSphere()
		//
		// Computes a sphere that contains all vertices (not necessarily the smallest
		// sphere). The returned object has two properties, `center` and `radius`.
		getBoundingSphere: function() {
			var aabb = this.getAABB();
			var sphere = { center: aabb.getCenter(), radius: 0 };
			for (var i = 0; i < this.vertices.length; i++) {
				sphere.radius = Math.max(sphere.radius, this.vertices[i].minus(sphere.center).length())
			}
			return sphere;
		}
	};

// ### GL.Mesh.plane([options])
//
// Generates a square 2x2 mesh the xy plane centered at the origin. The
// `options` argument specifies options to pass to the mesh constructor.
// Additional options include `detailX` and `detailY`, which set the tesselation
// in x and y, and `detail`, which sets both `detailX` and `detailY` at once.
// Two triangles are generated by default.
// Example usage:
//
//     var mesh1 = GL.Mesh.plane();
//     var mesh2 = GL.Mesh.plane({ detail: 5 });
//     var mesh3 = GL.Mesh.plane({ detailX: 20, detailY: 40 });
//
	Mesh.plane = function(options) {
		options = options || {};
		var mesh = new Mesh(options);
		var detailX = options.detailX || options.detail || 1;
		var detailY = options.detailY || options.detail || 1;

		for (var y = 0; y <= detailY; y++) {
			var t = y / detailY;
			for (var x = 0; x <= detailX; x++) {
				var s = x / detailX;
				mesh.vertices.push([2 * s - 1, 2 * t - 1, 0]);
				if (mesh.coords) mesh.coords.push([s, t]);
				if (mesh.normals) mesh.normals.push([0, 0, 1]);
				if (x < detailX && y < detailY) {
					var i = x + y * (detailX + 1);
					mesh.triangles.push([i, i + 1, i + detailX + 1]);
					mesh.triangles.push([i + detailX + 1, i + 1, i + detailX + 2]);
				}
			}
		}

		mesh.compile();
		return mesh;
	};

	var cubeData = [
		[0, 4, 2, 6, -1, 0, 0], // -x
		[1, 3, 5, 7, +1, 0, 0], // +x
		[0, 1, 4, 5, 0, -1, 0], // -y
		[2, 6, 3, 7, 0, +1, 0], // +y
		[0, 2, 1, 3, 0, 0, -1], // -z
		[4, 5, 6, 7, 0, 0, +1]  // +z
	];

	function pickOctant(i) {
		return V3((i & 1) * 2 - 1, (i & 2) - 1, (i & 4) / 2 - 1);
	}

// ### GL.Mesh.cube([options])
//
// Generates a 2x2x2 box centered at the origin. The `options` argument
// specifies options to pass to the mesh constructor.
	Mesh.cube = function(options) {
		var mesh = new Mesh(options);

		for (var i = 0; i < cubeData.length; i++) {
			var data = cubeData[i], v = i * 4;
			for (var j = 0; j < 4; j++) {
				var d = data[j];
				mesh.vertices.push(pickOctant(d).toArray());
				if (mesh.coords) mesh.coords.push([j & 1, (j & 2) / 2]);
				if (mesh.normals) mesh.normals.push(data.slice(4, 7));
			}
			mesh.triangles.push([v, v + 1, v + 2]);
			mesh.triangles.push([v + 2, v + 1, v + 3]);
		}

		mesh.compile();
		return mesh;
	};

// ### GL.Mesh.sphere([options])
//
// Generates a geodesic sphere of radius 1. The `options` argument specifies
// options to pass to the mesh constructor in addition to the `detail` option,
// which controls the tesselation level. The detail is `6` by default.
// Example usage:
//
//     var mesh1 = GL.Mesh.sphere();
//     var mesh2 = GL.Mesh.sphere({ detail: 2 });
//
	Mesh.sphere = function(options) {
		function tri(a, b, c) { return flip ? [a, c, b] : [a, b, c]; }
		function fix(x) { return x + (x - x * x) / 2; }
		options = options || {};
		var mesh = new Mesh(options);
		var indexer = new Indexer();
		var detail = options.detail || 6;

		for (var octant = 0; octant < 8; octant++) {
			var scale = pickOctant(octant);
			var flip = scale.x * scale.y * scale.z > 0;
			var data = [];
			for (var i = 0; i <= detail; i++) {
				// Generate a row of vertices on the surface of the sphere
				// using barycentric coordinates.
				for (var j = 0; i + j <= detail; j++) {
					var a = i / detail;
					var b = j / detail;
					var c = (detail - i - j) / detail;
					var vertex = { vertex: V3(fix(a), fix(b), fix(c)).normalized().multiply(scale).toArray() };
					if (mesh.coords) vertex.coord = scale.y > 0 ? [1 - a, c] : [c, 1 - a];
					data.push(indexer.add(vertex));
				}

				// Generate triangles from this row and the previous row.
				if (i > 0) {
					for (var j = 0; i + j <= detail; j++) {
						var a = (i - 1) * (detail + 1) + ((i - 1) - (i - 1) * (i - 1)) / 2 + j;
						var b = i * (detail + 1) + (i - i * i) / 2 + j;
						mesh.triangles.push(tri(data[a], data[a + 1], data[b]));
						if (i + j < detail) {
							mesh.triangles.push(tri(data[b], data[a + 1], data[b + 1]));
						}
					}
				}
			}
		}

		// Reconstruct the geometry from the indexer.
		mesh.vertices = indexer.unique.map(function(v) { return v.vertex; });
		if (mesh.coords) mesh.coords = indexer.unique.map(function(v) { return v.coord; });
		if (mesh.normals) mesh.normals = mesh.vertices;
		mesh.compile();
		return mesh;
	};
	/**
	 * Returns a sphere mesh with radius 1 created by subdividing the faces of a dodecahedron (20-sided)
	 * The sphere is positioned at the origin
	 * @param subdivisions
	 *      How many recursive divisions to do. A subdivision divides a triangle into 4,
	 *      so the total number of triangles is 20 * 4^subdivisions
	 * @returns {Mesh}
	 *      Contains vertex and normal buffers and index buffers for triangles and lines
	 */
	Mesh.sphere2 = function (subdivisions) {
		var golden = (1 + Math.sqrt(5)) / 2, u = V3(1, golden, 0).normalized(), s = u.x, t = u.y
		// base vertices of dodecahedron
		var vertices = [
			V3(-s,  t,  0),
			V3( s,  t,  0),
			V3(-s, -t,  0),
			V3( s, -t,  0),

			V3( 0, -s,  t),
			V3( 0,  s,  t),
			V3( 0, -s, -t),
			V3( 0,  s, -t),

			V3( t,  0, -s),
			V3( t,  0,  s),
			V3(-t,  0, -s),
			V3(-t,  0,  s),]
		// base triangles of dodecahedron
		var triangles = [
			0, 11, 5,
			0, 5, 1,
			0, 1, 7,
			0, 7, 10,
			0, 10, 11,

			// 5 adjacent faces
			1, 5, 9,
			5, 11, 4,
			11, 10, 2,
			10, 7, 6,
			7, 1, 8,

			// 5 faces around point 3
			3, 9, 4,
			3, 4, 2,
			3, 2, 6,
			3, 6, 8,
			3, 8, 9,

			// 5 adjacent faces
			4, 9, 5,
			2, 4, 11,
			6, 2, 10,
			8, 6, 7,
			9, 8, 1,
		]

		/**
		 * Tesselates triangle a b c
		 * a b c must already be in vertices with the indexes ia ib ic
		 * res is the number of subdivisions to do. 0 just results in triangle and line indexes being added to the
		 * respective buffers.
		 * @param a
		 * @param b
		 * @param c
		 * @param res
		 * @param vertices
		 * @param triangles
		 * @param ia
		 * @param ib
		 * @param ic
		 * @param lines
		 */
		function tesselateRecursively(a, b, c, res, vertices, triangles, ia, ib, ic, lines) {
			if (0 == res) {
				triangles.push([ia, ib, ic])
				if (ia < ib) lines.push(ia, ib)
				if (ib < ic) lines.push(ib, ic)
				if (ic < ia) lines.push(ic, ia)
			} else {
				// subdivide the triangle abc into 4 by adding a vertex (with the correct distance from the origin)
				// between each segment ab, bc and cd, then calling the function recursively
				var abMid1 = a.plus(b).toLength(1), bcMid1 = b.plus(c).toLength(1), caMid1 = c.plus(a).toLength(1)
				var iabm = vertices.length, ibcm = iabm + 1, icam = iabm + 2
				vertices.push(abMid1, bcMid1, caMid1)
				tesselateRecursively(abMid1, bcMid1, caMid1, res - 1, vertices, triangles, iabm, ibcm, icam, lines)
				tesselateRecursively(a, abMid1, caMid1, res - 1, vertices, triangles, ia, iabm, icam, lines)
				tesselateRecursively(b, bcMid1, abMid1, res - 1, vertices, triangles, ib, ibcm, iabm, lines)
				tesselateRecursively(c, caMid1, bcMid1, res - 1, vertices, triangles, ic, icam, ibcm, lines)
			}
		}

		var mesh = new GL.Mesh({normals: true, colors: false, lines: true});

		mesh.vertices.pushAll(vertices)
		for (var i = 0; i < 20; i++) {
			var [ia, ic, ib] = triangles.slice(i * 3, i * 3 + 3)
			tesselateRecursively(vertices[ia], vertices[ic], vertices[ib], subdivisions, mesh.vertices, mesh.triangles, ia, ic, ib, mesh.lines)
		}

		mesh.normals = mesh.vertices
		mesh.compile()
		return mesh

	}

// ### GL.Mesh.load(json[, options])
//
// Creates a mesh from the JSON generated by the `convert/convert.py` script.
// Example usage:
//
//     var data = {
//       vertices: [[0, 0, 0], [1, 0, 0], [0, 1, 0]],
//       triangles: [[0, 1, 2]]
//     };
//     var mesh = GL.Mesh.load(data);
//
	Mesh.load = function(json, options) {
		options = options || {};
		if (!('coords' in options)) options.coords = !!json.coords;
		if (!('normals' in options)) options.normals = !!json.normals;
		if (!('colors' in options)) options.colors = !!json.colors;
		if (!('triangles' in options)) options.triangles = !!json.triangles;
		if (!('lines' in options)) options.lines = !!json.lines;
		var mesh = new Mesh(options);
		mesh.vertices = json.vertices;
		if (mesh.coords) mesh.coords = json.coords;
		if (mesh.normals) mesh.normals = json.normals;
		if (mesh.colors) mesh.colors = json.colors;
		if (mesh.triangles) mesh.triangles = json.triangles;
		if (mesh.lines) mesh.lines = json.lines;
		mesh.compile();
		return mesh;
	};

// src/raytracer.js
// Provides a convenient raytracing interface.

// ### new GL.HitTest([t, hit, normal])
//
// This is the object used to return hit test results. If there are no
// arguments, the constructed argument represents a hit infinitely far
// away.
	function HitTest(t, hit, normal) {
		this.t = arguments.length ? t : Number.MAX_VALUE;
		this.hit = hit;
		this.normal = normal;
	}

// ### .mergeWith(other)
//
// Changes this object to be the closer of the two hit test results.
	HitTest.prototype = {
		mergeWith: function(other) {
			if (other.t > 0 && other.t < this.t) {
				this.t = other.t;
				this.hit = other.hit;
				this.normal = other.normal;
			}
		}
	};

// ### new GL.Raytracer()
//
// This will read the current modelview matrix, projection matrix, and viewport,
// reconstruct the eye position, and store enough information to later generate
// per-pixel rays using `getRayForPixel()`.
//
// Example usage:
//
//     var tracer = new GL.Raytracer();
//     var ray = tracer.getRayForPixel(
//       gl.canvas.width / 2,
//       gl.canvas.height / 2);
//     var result = GL.Raytracer.hitTestSphere(
//       tracer.eye, ray, new GL.Vector(0, 0, 0), 1);
	function Raytracer() {
		var v = gl.getParameter(gl.VIEWPORT);
		var m = gl.modelviewMatrix.m;

		var axisX = V3(m[0], m[4], m[8]);
		var axisY = V3(m[1], m[5], m[9]);
		var axisZ = V3(m[2], m[6], m[10]);
		var offset = V3(m[3], m[7], m[11]);
		this.eye = V3(-offset.dot(axisX), -offset.dot(axisY), -offset.dot(axisZ));

		var minX = v[0], maxX = minX + v[2];
		var minY = v[1], maxY = minY + v[3];
		this.ray00 = gl.unProject(minX, minY, 1).minus(this.eye);
		this.ray10 = gl.unProject(maxX, minY, 1).minus(this.eye);
		this.ray01 = gl.unProject(minX, maxY, 1).minus(this.eye);
		this.ray11 = gl.unProject(maxX, maxY, 1).minus(this.eye);
		this.viewport = v;
	}

	Raytracer.prototype = {
		// ### .getRayForPixel(x, y)
		//
		// Returns the ray originating from the camera and traveling through the pixel `x, y`.
		getRayForPixel: function(x, y) {
			x = (x - this.viewport[0]) / this.viewport[2];
			y = 1 - (y - this.viewport[1]) / this.viewport[3];
			var ray0 = Vector.lerp(this.ray00, this.ray10, x);
			var ray1 = Vector.lerp(this.ray01, this.ray11, x);
			return Vector.lerp(ray0, ray1, y).normalized();
		}
	};

// ### GL.Raytracer.hitTestBox(origin, ray, min, max)
//
// Traces the ray starting from `origin` along `ray` against the axis-aligned box
// whose coordinates extend from `min` to `max`. Returns a `HitTest` with the
// information or `null` for no intersection.
//
// This implementation uses the [slab intersection method](http://www.siggraph.org/education/materials/HyperGraph/raytrace/rtinter3.htm).
	Raytracer.hitTestBox = function(origin, ray, min, max) {
		var tMin = min.minus(origin).divide(ray);
		var tMax = max.minus(origin).divide(ray);
		var t1 = Vector.min(tMin, tMax);
		var t2 = Vector.max(tMin, tMax);
		var tNear = t1.max();
		var tFar = t2.min();

		if (tNear > 0 && tNear < tFar) {
			var epsilon = 1.0e-6, hit = origin.add(ray.multiply(tNear));
			min = min.add(epsilon);
			max = max.minus(epsilon);
			return new HitTest(tNear, hit, V3(
				(hit.x > max.x) - (hit.x < min.x),
				(hit.y > max.y) - (hit.y < min.y),
				(hit.z > max.z) - (hit.z < min.z)
			));
		}

		return null;
	};

// ### GL.Raytracer.hitTestSphere(origin, ray, center, radius)
//
// Traces the ray starting from `origin` along `ray` against the sphere defined
// by `center` and `radius`. Returns a `HitTest` with the information or `null`
// for no intersection.
	Raytracer.hitTestSphere = function(origin, ray, center, radius) {
		var offset = origin.minus(center);
		var a = ray.dot(ray);
		var b = 2 * ray.dot(offset);
		var c = offset.dot(offset) - radius * radius;
		var discriminant = b * b - 4 * a * c;

		if (discriminant > 0) {
			var t = (-b - Math.sqrt(discriminant)) / (2 * a), hit = origin.add(ray.multiply(t));
			return new HitTest(t, hit, hit.minus(center).divide(radius));
		}

		return null;
	};

// ### GL.Raytracer.hitTestTriangle(origin, ray, a, b, c)
//
// Traces the ray starting from `origin` along `ray` against the triangle defined
// by the points `a`, `b`, and `c`. Returns a `HitTest` with the information or
// `null` for no intersection.
	Raytracer.hitTestTriangle = function(origin, ray, a, b, c) {
		var ab = b.minus(a);
		var ac = c.minus(a);
		var normal = ab.cross(ac).normalized();
		var t = normal.dot(a.minus(origin)) / normal.dot(ray);

		if (t > 0) {
			var hit = origin.add(ray.multiply(t));
			var toHit = hit.minus(a);
			var dot00 = ac.dot(ac);
			var dot01 = ac.dot(ab);
			var dot02 = ac.dot(toHit);
			var dot11 = ab.dot(ab);
			var dot12 = ab.dot(toHit);
			var divide = dot00 * dot11 - dot01 * dot01;
			var u = (dot11 * dot02 - dot01 * dot12) / divide;
			var v = (dot00 * dot12 - dot01 * dot02) / divide;
			if (u >= 0 && v >= 0 && u + v <= 1) return new HitTest(t, hit, normal);
		}

		return null;
	};

// src/shader.js
// Provides a convenient wrapper for WebGL shaders. A few uniforms and attributes,
// prefixed with `gl_`, are automatically added to all shader sources to make
// simple shaders easier to write.
//
// Example usage:
//
//     var shader = new GL.Shader('\
//       void main() {\
//         gl_Position = gl_ModelViewProjectionMatrix * gl_Vertex;\
//       }\
//     ', '\
//       uniform vec4 color;\
//       void main() {\
//         gl_FragColor = color;\
//       }\
//     ');
//
//     shader.uniforms({
//       color: [1, 0, 0, 1]
//     }).draw(mesh);

	function regexMap(regex, text, callback) {
		var result;
		while ((result = regex.exec(text)) != null) {
			callback(result);
		}
	}

// Non-standard names beginning with `gl_` must be mangled because they will
// otherwise cause a compiler error.
	var LIGHTGL_PREFIX = 'LIGHTGL';

// ### new GL.Shader(vertexSource, fragmentSource)
//
// Compiles a shader program using the provided vertex and fragment shaders.
	function Shader(vertexSource, fragmentSource) {
		// Allow passing in the id of an HTML script tag with the source
		function followScriptTagById(id) {
			var element = document.getElementById(id);
			return element ? element.text : id;
		}
		vertexSource = followScriptTagById(vertexSource);
		fragmentSource = followScriptTagById(fragmentSource);

		// Headers are prepended to the sources to provide some automatic functionality.
		var header = `
	uniform mat3 gl_NormalMatrix;
	uniform mat4 gl_ModelViewMatrix;
	uniform mat4 gl_ProjectionMatrix;
	uniform mat4 gl_ModelViewProjectionMatrix;
	uniform mat4 gl_ModelViewMatrixInverse;
	uniform mat4 gl_ProjectionMatrixInverse;
	uniform mat4 gl_ModelViewProjectionMatrixInverse;
`
		var vertexHeader = header + `
	attribute vec4 gl_Vertex;
	attribute vec4 gl_TexCoord;
	attribute vec3 gl_Normal;
	attribute vec4 gl_Color;
	vec4 ftransform() {
		return gl_ModelViewProjectionMatrix * gl_Vertex;
	}
`
		var fragmentHeader = `  precision highp float;` + header;

		// Check for the use of built-in matrices that require expensive matrix
		// multiplications to compute, and record these in `usedMatrices`.
		var source = vertexSource + fragmentSource;
		var usedMatrices = {};
		regexMap(/\b(gl_[^;]*)\b;/g, header, function(groups) {
			var name = groups[1];
			if (source.indexOf(name) != -1) {
				var capitalLetters = name.replace(/[a-z_]/g, '');
				usedMatrices[capitalLetters] = LIGHTGL_PREFIX + name;
			}
		});
		if (source.indexOf('ftransform') != -1) usedMatrices.MVPM = LIGHTGL_PREFIX + 'gl_ModelViewProjectionMatrix';
		this.usedMatrices = usedMatrices;

		// The `gl_` prefix must be substituted for something else to avoid compile
		// errors, since it's a reserved prefix. This prefixes all reserved names with
		// `_`. The header is inserted after any extensions, since those must come
		// first.
		function fix(header, source) {
			var replaced = {};
			var match = /^((\s*\/\/.*\n|\s*#extension.*\n)+)[^]*$/.exec(source);
			source = match ? match[1] + header + source.substr(match[1].length) : header + source;
			regexMap(/\bgl_\w+\b/g, header, function(result) {
				if (!(result in replaced)) {
					source = source.replace(new RegExp('\\b' + result + '\\b', 'g'), LIGHTGL_PREFIX + result);
					replaced[result] = true;
				}
			});
			return source;
		}
		vertexSource = fix(vertexHeader, vertexSource);
		fragmentSource = fix(fragmentHeader, fragmentSource);

		// Compile and link errors are thrown as strings.
		function compileSource(type, source) {
			var shader = gl.createShader(type);
			gl.shaderSource(shader, source);
			gl.compileShader(shader);
			if (!gl.getShaderParameter(shader, gl.COMPILE_STATUS)) {
				throw new Error('compile error: ' + gl.getShaderInfoLog(shader));
			}
			return shader;
		}
		this.program = gl.createProgram();
		gl.attachShader(this.program, compileSource(gl.VERTEX_SHADER, vertexSource));
		gl.attachShader(this.program, compileSource(gl.FRAGMENT_SHADER, fragmentSource));
		gl.linkProgram(this.program);
		if (!gl.getProgramParameter(this.program, gl.LINK_STATUS)) {
			throw new Error('link error: ' + gl.getProgramInfoLog(this.program));
		}
		this.attributes = {};
		this.uniformLocations = {};

		// Sampler uniforms need to be uploaded using `gl.uniform1i()` instead of `gl.uniform1f()`.
		// To do this automatically, we detect and remember all uniform samplers in the source code.
		var isSampler = {};
		regexMap(/uniform\s+sampler(1D|2D|3D|Cube)\s+(\w+)\s*;/g, vertexSource + fragmentSource, function(groups) {
			isSampler[groups[2]] = 1;
		});
		this.isSampler = isSampler;

		if (NLA.DEBUG) {
			this.uniformInfo = {}
			for (var i = gl.getProgramParameter(this.program, gl.ACTIVE_UNIFORMS); i-- > 0;) {
				var info = gl.getActiveUniform(this.program, i)
				this.uniformInfo[info.name] = info
			}
		}
	}

	function isArray(obj) {
		var type = typeof obj
		//console.log(type, type.toString(), obj.toString(), obj, Object.prototype.toString.apply(obj))
		return obj.constructor == Array || Float32Array == obj.constructor || Float64Array == obj.constructor
	}

	function isNumber(obj) {
		var str = Object.prototype.toString.call(obj);
		return str == '[object Number]' || str == '[object Boolean]';
	}

	var tempMatrix = M4();
	var resultMatrix = M4();

	Shader.prototype = {
		// ### .uniforms(uniforms)
		//
		// Set a uniform for each property of `uniforms`. The correct `gl.uniform*()` method is
		// inferred from the value types and from the stored uniform sampler flags.
		uniforms: function(uniforms) {
			gl.useProgram(this.program);

			for (var name in uniforms) {
				var types = ["FLOAT", "FLOAT_MAT2", "FLOAT_MAT3", "FLOAT_MAT4", "FLOAT_VEC2", "FLOAT_VEC3", "FLOAT_VEC4", "INT", "INT_VEC2", "INT_VEC3", "INT_VEC4", "UNSIGNED_INT"]
				var location = this.uniformLocations[name] || gl.getUniformLocation(this.program, name);
				assert(!!location, name)
				if (!location) continue;
				this.uniformLocations[name] = location;
				var value = uniforms[name];
				if (NLA.DEBUG) {
					var info = this.uniformInfo[name]
					assert(info.type != gl.FLOAT || "number" == typeof value)
					assert(info.type != gl.INT || "number" == typeof value && value % 1 == 0)
					assert(info.type != gl.FLOAT_VEC3 || value instanceof NLA.Vector3)
					assert(info.type != gl.FLOAT_MAT4 || value instanceof NLA.Matrix4x4, () => value.toSource())
					assert(info.type != gl.FLOAT_MAT3 || value.length == 9)
				}
				if (value instanceof V3) {
					value = value.toArray()
				} else if (value instanceof M4) {
					value = value.m;
				}
				if (isArray(value)) {
					switch (value.length) {
						case 1: gl.uniform1fv(location, new Float32Array(value)); break;
						case 2: gl.uniform2fv(location, new Float32Array(value)); break;
						case 3: gl.uniform3fv(location, new Float32Array(value)); break;
						case 4: gl.uniform4fv(location, new Float32Array(value)); break;
						// Matrices are automatically transposed, since WebGL uses column-major
						// indices instead of row-major indices.
						case 9: gl.uniformMatrix3fv(location, false, new Float32Array([
							value[0], value[3], value[6],
							value[1], value[4], value[7],
							value[2], value[5], value[8]
						])); break;
						case 16: gl.uniformMatrix4fv(location, false, new Float32Array([
							value[0], value[4], value[8], value[12],
							value[1], value[5], value[9], value[13],
							value[2], value[6], value[10], value[14],
							value[3], value[7], value[11], value[15]
						])); break;
						default: throw new Error('don\'t know how to load uniform "' + name + '" of length ' + value.length);
					}
				} else if (isNumber(value)) {
					(this.isSampler[name] ? gl.uniform1i : gl.uniform1f).call(gl, location, value);
				} else {
					throw new Error('attempted to set uniform "' + name + '" to invalid value ' + value);
				}
			}

			return this;
		},

		// ### .draw(mesh[, mode])
		//
		// Sets all uniform matrix attributes, binds all relevant buffers, and draws the
		// mesh geometry as indexed triangles or indexed lines. Set `mode` to `gl.LINES`
		// (and either add indices to `lines` or call `computeWireframe()`) to draw the
		// mesh in wireframe.
		draw: function(mesh, mode, start, count) {
			assert(mesh.hasBeenCompiled, 'mesh.hasBeenCompiled')
			mode = mode || 'TRIANGLES'
			assert(GL.drawModes.includes(mode), 'GL.drawModes.includes(mode) ' + mode)
			this.drawBuffers(mesh.vertexBuffers, mesh.indexBuffers[mode], gl[mode], start, count);
		},

		// ### .drawBuffers(vertexBuffers, indexBuffer, mode)
		//
		// Sets all uniform matrix attributes, binds all relevant buffers, and draws the
		// indexed mesh geometry. The `vertexBuffers` argument is a map from attribute
		// names to `Buffer` objects of type `gl.ARRAY_BUFFER`, `indexBuffer` is a `Buffer`
		// object of type `gl.ELEMENT_ARRAY_BUFFER`, and `mode` is a WebGL primitive mode
		// like `gl.TRIANGLES` or `gl.LINES`. This method automatically creates and caches
		// vertex attribute pointers for attributes as needed.
		drawBuffers: function(vertexBuffers, indexBuffer, mode, start, count) {
			assert(GL.drawModes.map(m => gl[m]).includes(mode), 'GL.drawModes.map(m => gl[m]).includes(mode) ' + mode)
			// Only construct up the built-in matrices we need for this shader.
			var used = this.usedMatrices;
			var MVM = gl.modelviewMatrix;
			var PM = gl.projectionMatrix;
			var MVMI = (used.MVMI || used.NM) ? MVM.inversed() : null;
			var PMI = (used.PMI) ? PM.inversed() : null;
			var MVPM = (used.MVPM || used.MVPMI) ? PM.times(MVM) : null;
			var matrices = {};
			if (used.MVM) matrices[used.MVM] = MVM;
			if (used.MVMI) matrices[used.MVMI] = MVMI;
			if (used.PM) matrices[used.PM] = PM;
			if (used.PMI) matrices[used.PMI] = PMI;
			if (used.MVPM) matrices[used.MVPM] = MVPM;
			if (used.MVPMI) matrices[used.MVPMI] = MVPM.inversed();
			if (used.NM) {
				var m = MVMI.m;
				// transpose normal matrix
				matrices[used.NM] = [m[0], m[4], m[8], m[1], m[5], m[9], m[2], m[6], m[10]];
			}
			this.uniforms(matrices);

			// Create and enable attribute pointers as necessary.
			var length = 0;
			for (var attribute in vertexBuffers) {
				var buffer = vertexBuffers[attribute];
				var location = this.attributes[attribute] ||
					gl.getAttribLocation(this.program, attribute.replace(/^(gl_.*)$/, LIGHTGL_PREFIX + '$1'));
				if (location == -1 || !buffer.buffer) continue;
				this.attributes[attribute] = location;
				gl.bindBuffer(gl.ARRAY_BUFFER, buffer.buffer, buffer.buffer.spacing, buffer.buffer.length);
				gl.enableVertexAttribArray(location);
				//console.log(attribute)
				try {
					gl.vertexAttribPointer(location, buffer.buffer.spacing, gl.FLOAT, false, 0, 0);
				} catch (e) {
					console.log("caught")
				}

				// TODO: uuuh
				length = buffer.buffer.count
			}

			// Disable unused attribute pointers.
			for (var attribute in this.attributes) {
				if (!(attribute in vertexBuffers)) {
					gl.disableVertexAttribArray(this.attributes[attribute]);
				}
			}

			// Draw the geometry.
			if (length && (!indexBuffer || indexBuffer.buffer)) {
				if (indexBuffer) {
					if (start + count > indexBuffer.buffer.length) {
						throw new Error("Buffer not long enough for passed parameters start/length/buffer length" +" "+ start +" "+ count + " "+ indexBuffer.buffer.length)
					}
					gl.bindBuffer(gl.ELEMENT_ARRAY_BUFFER, indexBuffer.buffer);
					gl.drawElements(mode, count || indexBuffer.buffer.length, gl.UNSIGNED_SHORT, 2 * (start || 0));
				} else {
					if (start + count > length) {
						throw new Error("invalid")
					}
					//console.log("drawArrays", mode, start || 0, count || length);
					gl.drawArrays(mode, start || 0, count || length)
				}
			}

			return this;
		}
	};
	GL.drawModes = ['POINTS', 'LINES', 'LINE_STRIP', 'LINE_LOOP', 'TRIANGLES', 'TRIANGLE_STRIP', 'TRIANGLE_FAN']
// src/texture.js
// Provides a simple wrapper around WebGL textures that supports render-to-texture.

// ### new GL.Texture(width, height[, options])
//
// The arguments `width` and `height` give the size of the texture in texels.
// WebGL texture dimensions must be powers of two unless `filter` is set to
// either `gl.NEAREST` or `gl.LINEAR` and `wrap` is set to `gl.CLAMP_TO_EDGE`
// (which they are by default).
//
// Texture parameters can be passed in via the `options` argument.
// Example usage:
//
//     var t = new GL.Texture(256, 256, {
//       // Defaults to gl.LINEAR, set both at once with "filter"
//       magFilter: gl.NEAREST,
//       minFilter: gl.LINEAR,
//
//       // Defaults to gl.CLAMP_TO_EDGE, set both at once with "wrap"
//       wrapS: gl.REPEAT,
//       wrapT: gl.REPEAT,
//
//       format: gl.RGB, // Defaults to gl.RGBA
//       type: gl.FLOAT // Defaults to gl.UNSIGNED_BYTE
//     });
	function Texture(width, height, options) {
		options = options || {};
		this.id = gl.createTexture();
		this.width = width;
		this.height = height;
		this.format = options.format || gl.RGBA;
		this.type = options.type || gl.UNSIGNED_BYTE;
		var magFilter = options.filter || options.magFilter || gl.LINEAR;
		var minFilter = options.filter || options.minFilter || gl.LINEAR;
		if (this.type === gl.FLOAT) {
			if (!Texture.canUseFloatingPointTextures()) {
				throw new Error('OES_texture_float is required but not supported');
			}
			if ((minFilter !== gl.NEAREST || magFilter !== gl.NEAREST) &&
				!Texture.canUseFloatingPointLinearFiltering()) {
				throw new Error('OES_texture_float_linear is required but not supported');
			}
		} else if (this.type === gl.HALF_FLOAT_OES) {
			if (!Texture.canUseHalfFloatingPointTextures()) {
				throw new Error('OES_texture_half_float is required but not supported');
			}
			if ((minFilter !== gl.NEAREST || magFilter !== gl.NEAREST) &&
				!Texture.canUseHalfFloatingPointLinearFiltering()) {
				throw new Error('OES_texture_half_float_linear is required but not supported');
			}
		}
		gl.bindTexture(gl.TEXTURE_2D, this.id);
		gl.pixelStorei(gl.UNPACK_FLIP_Y_WEBGL, 1);
		gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_MAG_FILTER, magFilter);
		gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_MIN_FILTER, minFilter);
		gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_WRAP_S, options.wrap || options.wrapS || gl.CLAMP_TO_EDGE);
		gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_WRAP_T, options.wrap || options.wrapT || gl.CLAMP_TO_EDGE);
		gl.texImage2D(gl.TEXTURE_2D, 0, this.format, width, height, 0, this.format, this.type, null);
	}

	var framebuffer;
	var renderbuffer;
	var checkerboardCanvas;

	Texture.prototype = {
		// ### .bind([unit])
		//
		// Bind this texture to the given texture unit (0-7, defaults to 0).
		bind: function(unit) {
			gl.activeTexture(gl.TEXTURE0 + (unit || 0));
			gl.bindTexture(gl.TEXTURE_2D, this.id);
		},

		// ### .unbind([unit])
		//
		// Clear the given texture unit (0-7, defaults to 0).
		unbind: function(unit) {
			gl.activeTexture(gl.TEXTURE0 + (unit || 0));
			gl.bindTexture(gl.TEXTURE_2D, null);
		},

		// ### .canDrawTo()
		//
		// Check if rendering to this texture is supported. It may not be supported
		// for floating-point textures on some configurations.
		canDrawTo: function() {
			framebuffer = framebuffer || gl.createFramebuffer();
			gl.bindFramebuffer(gl.FRAMEBUFFER, framebuffer);
			gl.framebufferTexture2D(gl.FRAMEBUFFER, gl.COLOR_ATTACHMENT0, gl.TEXTURE_2D, this.id, 0);
			var result = gl.checkFramebufferStatus(gl.FRAMEBUFFER) == gl.FRAMEBUFFER_COMPLETE;
			gl.bindFramebuffer(gl.FRAMEBUFFER, null);
			return result;
		},

		// ### .drawTo(callback)
		//
		// Render all draw calls in `callback` to this texture. This method sets up
		// a framebuffer with this texture as the color attachment and a renderbuffer
		// as the depth attachment. It also temporarily changes the viewport to the
		// size of the texture.
		//
		// Example usage:
		//
		//     texture.drawTo(function() {
		//       gl.clearColor(1, 0, 0, 1);
		//       gl.clear(gl.COLOR_BUFFER_BIT);
		//     });
		drawTo: function(callback) {
			var v = gl.getParameter(gl.VIEWPORT);
			framebuffer = framebuffer || gl.createFramebuffer();
			renderbuffer = renderbuffer || gl.createRenderbuffer();
			gl.bindFramebuffer(gl.FRAMEBUFFER, framebuffer);
			gl.bindRenderbuffer(gl.RENDERBUFFER, renderbuffer);
			if (this.width != renderbuffer.width || this.height != renderbuffer.height) {
				renderbuffer.width = this.width;
				renderbuffer.height = this.height;
				gl.renderbufferStorage(gl.RENDERBUFFER, gl.DEPTH_COMPONENT16, this.width, this.height);
			}
			gl.framebufferTexture2D(gl.FRAMEBUFFER, gl.COLOR_ATTACHMENT0, gl.TEXTURE_2D, this.id, 0);
			gl.framebufferRenderbuffer(gl.FRAMEBUFFER, gl.DEPTH_ATTACHMENT, gl.RENDERBUFFER, renderbuffer);
			if (gl.checkFramebufferStatus(gl.FRAMEBUFFER) != gl.FRAMEBUFFER_COMPLETE) {
				throw new Error('Rendering to this texture is not supported (incomplete framebuffer)');
			}
			gl.viewport(0, 0, this.width, this.height);

			callback();

			gl.bindFramebuffer(gl.FRAMEBUFFER, null);
			gl.bindRenderbuffer(gl.RENDERBUFFER, null);
			gl.viewport(v[0], v[1], v[2], v[3]);
		},

		// ### .swapWith(other)
		//
		// Switch this texture with `other`, useful for the ping-pong rendering
		// technique used in multi-stage rendering.
		swapWith: function(other) {
			var temp;
			temp = other.id; other.id = this.id; this.id = temp;
			temp = other.width; other.width = this.width; this.width = temp;
			temp = other.height; other.height = this.height; this.height = temp;
		}
	};

// ### GL.Texture.fromImage(image[, options])
//
// Return a new image created from `image`, an `<img>` tag.
	Texture.fromImage = function(image, options) {
		options = options || {};
		var texture = new Texture(image.width, image.height, options);
		try {
			gl.texImage2D(gl.TEXTURE_2D, 0, texture.format, texture.format, texture.type, image);
		} catch (e) {
			if (location.protocol == 'file:') {
				throw new Error('image not loaded for security reasons (serve this page over "http://" instead)');
			} else {
				throw new Error('image not loaded for security reasons (image must originate from the same ' +
					'domain as this page or use Cross-Origin Resource Sharing)');
			}
		}
		if (options.minFilter && options.minFilter != gl.NEAREST && options.minFilter != gl.LINEAR) {
			gl.generateMipmap(gl.TEXTURE_2D);
		}
		return texture;
	};

// ### GL.Texture.fromURL(url[, options])
//
// Returns a checkerboard texture that will switch to the correct texture when
// it loads.
	Texture.fromURL = function(url, options) {
		checkerboardCanvas = checkerboardCanvas || (function() {
				var c = document.createElement('canvas').getContext('2d');
				c.canvas.width = c.canvas.height = 128;
				for (var y = 0; y < c.canvas.height; y += 16) {
					for (var x = 0; x < c.canvas.width; x += 16) {
						c.fillStyle = (x ^ y) & 16 ? '#FFF' : '#DDD';
						c.fillRect(x, y, 16, 16);
					}
				}
				return c.canvas;
			})();
		var texture = Texture.fromImage(checkerboardCanvas, options);
		var image = new Image();
		var context = gl;
		image.onload = function() {
			context.makeCurrent();
			Texture.fromImage(image, options).swapWith(texture);
		};
		image.src = url;
		return texture;
	};

// ### GL.Texture.canUseFloatingPointTextures()
//
// Returns false if `gl.FLOAT` is not supported as a texture type. This is the
// `OES_texture_float` extension.
	Texture.canUseFloatingPointTextures = function() {
		return !!gl.getExtension('OES_texture_float');
	};

// ### GL.Texture.canUseFloatingPointLinearFiltering()
//
// Returns false if `gl.LINEAR` is not supported as a texture filter mode for
// textures of type `gl.FLOAT`. This is the `OES_texture_float_linear`
// extension.
	Texture.canUseFloatingPointLinearFiltering = function() {
		return !!gl.getExtension('OES_texture_float_linear');
	};

// ### GL.Texture.canUseFloatingPointTextures()
//
// Returns false if `gl.HALF_FLOAT_OES` is not supported as a texture type.
// This is the `OES_texture_half_float` extension.
	Texture.canUseHalfFloatingPointTextures = function() {
		return !!gl.getExtension('OES_texture_half_float');
	};

// ### GL.Texture.canUseFloatingPointLinearFiltering()
//
// Returns false if `gl.LINEAR` is not supported as a texture filter mode for
// textures of type `gl.HALF_FLOAT_OES`. This is the
// `OES_texture_half_float_linear` extension.
	Texture.canUseHalfFloatingPointLinearFiltering = function() {
		return !!gl.getExtension('OES_texture_half_float_linear');
	};


	return GL;
})();
function createRectangleMesh(x0, y0, x1, y1) {
	var mesh = new GL.Mesh({});
	mesh.vertices.push(
		V3.create(x0, y0, 0),
		V3.create(x1, y0, 0),
		V3.create(x0, y1, 0),
		V3.create(x1, y1, 0)
	);
	mesh.triangles.push(
		[0, 1, 2],
		[3, 2, 1]
	);
	mesh.compile();
	return mesh;
}
function rotationMesh(vertices, lineAxis, totalAngle, count, close, normals) {
	var mesh = new GL.Mesh({normals: !!normals})
	var vc = vertices.length, vTotal = vc * count

	for (var i = 0; i < count; i++) {
		var angle = totalAngle / count * i
		var m = M4()
		M4.rotationLine(lineAxis.anchor, lineAxis.dir1, angle, m)
		Array.prototype.push.apply(mesh.vertices, m.transformedPoints(vertices))
		normals && Array.prototype.push.apply(mesh.normals, m.transformedVectors(normals))

		// add triangles
		for (var j = 0; j < vc - 1; j++) {
			mesh.triangles.push([i * vc + j + 1, i * vc + j, (i + 1) * vc + j].map(x => x % vTotal))
			mesh.triangles.push([i * vc + j + 1, (i + 1) * vc + j, (i + 1) * vc + j + 1].map(x => x % vTotal))
		}
	}

	//TODO: make sure normals dont need to be adjusted
	mesh.compile()
	return mesh
}
/**
 *
 * @param {Array.<number>} triangles
 * @param {boolean} flipped
 * @param {number} a
 * @param {number} b
 * @param {number} c
 * @param {number} d
 */
function pushQuad(triangles, flipped, a, b, c, d) {
	if (flipped) {
		triangles.push(a, c, b,
			b, c, d)
	} else {
		triangles.push(a, b, c,
			b, d, c)
	}
}
/**
 *
 * @param {Array.<V3>} vertices
 * @param {V3} offset
 * @param {boolean} close
 * @param {Array.<V3>} normals
 */
function offsetMesh(vertices, offset, close, normals) {
	assertVectors.apply(undefined, vertices)
	assertVectors(offset)

	let mesh = new GL.Mesh({normals: !!normals})
	mesh.vertices = vertices.concat(vertices.map(v => v.plus(offset)))

	for (let i = 0; i < vertices.length - 1; i++) {
		pushQuad(mesh.triangles, false,
			i, i + 1,
			vertices.length + i, vertices.length + i + 1)
	}
	if (close) {
		pushQuad(mesh.triangles, false, vertices.length - 1, vertices.length, vertices.length * 2 - 1)
	}
	if (normals) {
		mesh.normals = normals.concat(normals)
	}
	return mesh
}