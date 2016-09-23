/*
 * lightgl.js
 * http://github.com/evanw/lightgl.js/
 *
 * Copyright 2011 Evan Wallace
 * Released under the MIT license
 */
"use strict"
var GL = (function() {


	/**
	 * never called, used so IDE type checker doesn't pick a fit
	 *
	 * @name LightGLContext
	 * @class
	 * @property matrixMode
	 * @augments WebGLRenderingContext
	 */
	function LightGLContext() {
		this.modelViewMatrix = M4()
		this.projectionMatrix = M4()
		this.onupdate = null
		this.ondraw = null
		this.onclick = null
		this.onmousemove = null
		this.onmouseup = null
		this.onmousedown = null
	}
	LightGLContext.MODELVIEW = {}
	LightGLContext.PROJECTION = {}

	var gl

	var _GL = {
		/**
		 * `GL.create()` creates a new WebGL context and augments it with more methods. The alpha channel is disabled
		 * by default because it usually causes unintended transparencies in the canvas.
		 *
		 * @param {Object} options
		 * @param {HTMLCanvasElement=} options.canvas Canvas to use. A new one will be created if undefined.
		 * @param {boolean=} options.alpha
		 * @return {LightGLContext}
		 */
		create: function(options) {
			options = options || {};
			var canvas = options.canvas || document.createElement('canvas');
			if (!options.canvas) {
				canvas.width = 800;
				canvas.height = 600;
			}
			if (!('alpha' in options)) options.alpha = false;
			try { gl = canvas.getContext('webgl', options); } catch (e) { console.log(e)}
			try { gl = gl || canvas.getContext('experimental-webgl', options); } catch (e) {}
			if (!gl) throw new Error('WebGL not supported');
			gl.HALF_FLOAT_OES = 0x8D61;
			addMatrixStack(gl);
			addImmediateMode(gl);
			addEventListeners(gl);
			addOtherMethods(gl);
			return gl;
		},

		// `GL.keys` contains a mapping of key codes to booleans indicating whether
		// that key is currently pressed.
		keys: {},

		// Export all external classes.
		Buffer: Buffer,
		Mesh: Mesh,
		Shader: Shader,
		Texture: Texture,
	};

	/**
	 * Implement the OpenGL modelview and projection matrix stacks, along with some other useful GLU matrix functions.
	 * @param {LightGLContext} gl
	 */
	function addMatrixStack(gl) {
		var tempMatrix = M4()
		var resultMatrix = M4();
		gl.modelViewMatrix = M4();
		gl.projectionMatrix = M4();
		var modelViewStack = [];
		var projectionStack = [];
		var matrix, stack;


		NLA.addOwnProperties(gl, /** @lends LightGLContext.prototype */ {
			MODELVIEW: LightGLContext.MODELVIEW,
			PROJECTION: LightGLContext.PROJECTION,
			matrixMode: function(mode) {
				switch (mode) {
					case this.MODELVIEW:
						matrix = 'modelViewMatrix';
						stack = modelViewStack;
						break;
					case this.PROJECTION:
						matrix = 'projectionMatrix';
						stack = projectionStack;
						break;
					default:
						throw new Error('invalid matrix mode ' + mode);
				}
			},
			modelViewMode: function () {
				Object.defineProperty(gl, 'currentMatrix', {
					get: function () { return this.modelViewMatrix },
					set: function (val) { this.modelViewMatrix = val},
					writable: true})
				matrix = 'modelViewMatrix'
				stack = modelViewStack
			},
			projectionMode: function () {
				matrix = 'projectionMatrix'
				stack = projectionStack
			},
			loadIdentity: function() {
				M4.identity(gl[matrix]);
			},
			loadMatrix: function(m4) {
				M4.copy(m4, gl[matrix])
			},
			multMatrix: function(m) {
				M4.multiply(gl[matrix], m, resultMatrix);
				var temp = resultMatrix
				resultMatrix = gl[matrix]
				gl[matrix] = temp
			},
			scaleVector: function(x, y, z) {
				this.multMatrix(M4.scaleVector(x, y, z));
			},
			perspective: function(fov, aspect, near, far) {
				this.multMatrix(M4.perspective(fov, aspect, near, far, tempMatrix));
			},
			frustum: function(l, r, b, t, n, f) {
				this.multMatrix(M4.frustum(l, r, b, t, n, f, tempMatrix));
			},
			ortho: function(l, r, b, t, n, f) {
				this.multMatrix(M4.ortho(l, r, b, t, n, f, tempMatrix));
			},
			scale: function(x, y, z) {
				this.multMatrix(M4.scaling(x, y, z, tempMatrix));
			},
			translate: function(x, y, z) {
				if (undefined !== y) {
					this.multMatrix(M4.translation(x, y, z, tempMatrix));
				} else {
					this.multMatrix(M4.translation(x, tempMatrix));
				}
			},
			rotate: function(a, x, y, z) {
				this.multMatrix(M4.rotation(a, x, y, z, tempMatrix));
			},
			lookAt: function(eye, center, up) {
				this.multMatrix(M4.lookAt(eye, center, up, tempMatrix));
			},
			pushMatrix: function() {
				stack.push(M4.copy(gl[matrix]));
			},
			popMatrix: function() {
				gl[matrix] = stack.pop()
			},
			project: function(objX, objY, objZ, modelview, projection, viewport) {
				modelview = modelview || this.modelViewMatrix;
				projection = projection || this.projectionMatrix;
				viewport = viewport || this.getParameter(this.VIEWPORT);
				var point = projection.transformPoint(modelview.transformPoint(V3(objX, objY, objZ)));
				return V3(
					viewport[0] + viewport[2] * (point.x * 0.5 + 0.5),
					viewport[1] + viewport[3] * (point.y * 0.5 + 0.5),
					point.z * 0.5 + 0.5
				);
			},
			unProject: function(winX, winY, winZ, modelview, projection, viewport) {
				modelview = modelview || this.modelViewMatrix;
				projection = projection || this.projectionMatrix;
				viewport = viewport || this.getParameter(this.VIEWPORT);
				var point = V3(
					(winX - viewport[0]) / viewport[2] * 2 - 1,
					(winY - viewport[1]) / viewport[3] * 2 - 1,
					winZ * 2 - 1
				);
				return M4.inverse(M4.multiply(projection, modelview, tempMatrix), resultMatrix).transformPoint(point);
			},
		})
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

	function addImmediateMode(gl) {
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
	color = LGL_Color;
	coord = LGL_TexCoord;
	gl_Position = LGL_ModelViewProjectionMatrix * LGL_Vertex;
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
		NLA.addOwnProperties(gl, /** @lends LightGLContext.prototype */ {
			pointSize: function(pointSize) {
				immediateMode.shader.uniforms({ pointSize: pointSize });
			},
			begin: function(mode) {
				if (immediateMode.mode != -1) throw new Error('mismatched gl.begin() and gl.end() calls');
				immediateMode.mode = mode;
				immediateMode.mesh.colors = [];
				immediateMode.mesh.coords = [];
				immediateMode.mesh.vertices = [];
			},
			color: function(r, g, b, a) {
				immediateMode.color = (arguments.length == 1) ? r.toArray().concat(1) : [r, g, b, a || 1];
			},
			texCoord: function(s, t) {
				immediateMode.coord = (arguments.length == 1) ? s.toArray(2) : [s, t];
			},
			vertex: function(x, y, z) {
				immediateMode.mesh.colors.push(immediateMode.color);
				immediateMode.mesh.coords.push(immediateMode.coord);
				immediateMode.mesh.vertices.push(arguments.length == 1 ? x.toArray() : [x, y, z]);
			},
			end: function() {
				if (immediateMode.mode == -1) throw new Error('mismatched gl.begin() and gl.end() calls');
				console.log(immediateMode.mesh.toSource())
				immediateMode.mesh.compile();
				immediateMode.shader.uniforms({
					useTexture: !!gl.getParameter(gl.TEXTURE_BINDING_2D)
				}).draw(immediateMode.mesh, immediateMode.mode);
				immediateMode.mode = -1;
			}
		})
	}

// ### Improved mouse events
//
// This adds event listeners on the `gl.canvas` element that call
// `gl.onmousedown()`, `gl.onmousemove()`, and `gl.onmouseup()` with an
// augmented event object. The event object also has the properties `x`, `y`,
// `deltaX`, `deltaY`, and `dragging`.
	function addEventListeners(gl) {
		var oldX = 0, oldY = 0, buttons = {}, hasOld = false;
		function isDragging() {
			for (var b in buttons) {
				if (b in buttons && buttons[b]) return true;
			}
			return false;
		}
		function augmented(original) {
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
		function mousemove(e) {
			e = augmented(e);
			if (gl.onmousemove) gl.onmousemove(e);
			e.preventDefault();
		}
		function mouseup(e) {
			buttons[e.which] = false;
			if (!isDragging()) {
				// Shrink the event handlers back to the canvas when dragging ends.
				document.removeEventListener('mousemove', mousemove);
				document.removeEventListener('mouseup', mouseup);
				gl.canvas.addEventListener('mousemove', mousemove);
				gl.canvas.addEventListener('mouseup', mouseup);
			}
			e = augmented(e);
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
		gl.canvas.addEventListener('mousedown', function (e) {
			if (!isDragging()) {
				// Expand the event handlers to the document to handle dragging off canvas.
				document.addEventListener('mousemove', mousemove)
				document.addEventListener('mouseup', mouseup)
				gl.canvas.removeEventListener('mousemove', mousemove)
				gl.canvas.removeEventListener('mouseup', mouseup)
			}
			buttons[e.which] = true;
			e = augmented(e);
			if (gl.onmousedown) gl.onmousedown(e);
			e.preventDefault();
		})
		gl.canvas.addEventListener('mousemove', mousemove)
		gl.canvas.addEventListener('mouseup', mouseup)
		gl.canvas.addEventListener('mouseover', reset)
		gl.canvas.addEventListener('mouseout', reset)
		document.addEventListener('contextmenu', resetAll)
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

	document.addEventListener('keydown', function (e) {
		if (!e.altKey && !e.ctrlKey && !e.metaKey) {
			var key = mapKeyCode(e.keyCode);
			if (key) _GL.keys[key] = true;
			_GL.keys[e.keyCode] = true;
		}
	});
	document.addEventListener('keyup', function (e) {
		if (!e.altKey && !e.ctrlKey && !e.metaKey) {
			var key = mapKeyCode(e.keyCode);
			if (key) _GL.keys[key] = false;
			_GL.keys[e.keyCode] = false;
		}
	});
	function addOtherMethods(_gl) {
		/**
		 * @memberOf LightGLContext.prototype
		 */
		_gl.makeCurrent = function () {
			gl = this
		}



		// ### Animation
		//
		// Call `gl.animate()` to provide an animation loop that repeatedly calls
		// `gl.onupdate()` and `gl.ondraw()`.
		/**
		 * @memberOf LightGLContext
		 */
		_gl.animate = function() {
			var raf =
				window.requestAnimationFrame ||
				window.mozRequestAnimationFrame ||
				window.webkitRequestAnimationFrame ||
				function(callback) { setTimeout(callback, 1000 / 60); };
			var time = new Date().getTime();
			function update() {
				var now = new Date().getTime();
				if (_gl.onupdate) _gl.onupdate((now - time) / 1000);
				if (_gl.ondraw) _gl.ondraw();
				raf(update);
				time = now;
			}
			update();
		};

		/**
		 * Provide an easy way to get a fullscreen app running, including an
		 * automatic 3D perspective projection matrix by default. This should be
		 * called once.
		 *
		 * Just fullscreen, no automatic camera:
		 *
		 *     gl.fullscreen({ camera: false });
		 *
		 * Adjusting field of view, near plane distance, and far plane distance:
		 *
		 *     gl.fullscreen({ fov: 45, near: 0.1, far: 1000 });
		 *
		 * Adding padding from the edge of the window:
		 *
		 *     gl.fullscreen({ paddingLeft: 250, paddingBottom: 60 });
		 *
		 * @param {object=} options
		 * @param {number=} options.paddingTop
		 * @param {number=} options.paddingLeft
		 * @param {number=} options.paddingRight
		 * @param {number=} options.paddingBottom
		 *
		 * @param {boolean=} options.camera
		 * @param {number=} options.fov
		 * @param {number=} options.near
		 * @param {number=} options.far
		 *
		 * @memberOf LightGLContext.prototype
		 * @this LightGLContext
		 *
		 */
		_gl.fullscreen = function(options) {
			options = options || {};
			var top = options.paddingTop || 0;
			var left = options.paddingLeft || 0;
			var right = options.paddingRight || 0;
			var bottom = options.paddingBottom || 0;
			if (!document.body) {
				throw new Error('document.body doesn\'t exist yet (call gl.fullscreen() from ' +
					'window.onload() or from inside the <body> tag)');
			}
			document.body.appendChild(this.canvas);
			document.body.style.overflow = 'hidden';
			this.canvas.style.position = 'absolute';
			this.canvas.style.left = left + 'px';
			this.canvas.style.top = top + 'px';
			const gl = this
			function windowOnResize() {
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
			window.addEventListener('resize', windowOnResize)
			windowOnResize();
		};
	}

	/**
	 * Provides a simple method of uploading data to a GPU buffer. Example usage:
	 *
	 *     var vertices = new GL.Buffer(gl.ARRAY_BUFFER, Float32Array);
	 *     vertices.data = [[0, 0, 0], [1, 0, 0], [0, 1, 0], [1, 1, 0]];
	 *     vertices.compile();
	 *
	 *     var indices = new GL.Buffer(gl.ELEMENT_ARRAY_BUFFER, Uint16Array);
	 *     indices.data = [[0, 1, 2], [2, 1, 3]];
	 *     indices.compile();
	 *
	 * @param {number} target
	 * Specifies the target to which the buffer object is bound.
	 * The symbolic constant must be GL_ARRAY_BUFFER or GL_ELEMENT_ARRAY_BUFFER.
	 * @param type Float32Array or Uint16Array
	 *
	 * @property {WebGLBuffer} buffer
	 * @constructor
	 */
	function Buffer(target, type) {
		assert(target == WebGLRenderingContext.ARRAY_BUFFER || target == WebGLRenderingContext.ELEMENT_ARRAY_BUFFER, 'target == gl.ARRAY_BUFFER || target == gl.ELEMENT_ARRAY_BUFFER')
		assert(type == Float32Array || type == Uint16Array, 'type == Float32Array || type == Uint16Array')
		this.buffer = null
		this.target = target
		this.type = type
		this.data = []
		/**
		 * Number of elements in buffer. 2 V3s is still 2, not 6.
		 * @type {number}
		 */
		this.count = 0
		/**
		 * Space between elements in buffer. 3 for V3s.
		 * @type {number}
		 */
		this.spacing = 0
		this.hasBeenCompiled = false
	}

	Buffer.prototype = {
		/**
		 * Upload the contents of `data` to the GPU in preparation for rendering. The data must be a list of lists
		 * where each inner list has the same length. For example, each element of data for vertex normals would be a
		 * list of length three. This will remember the data length and element length for later use by shaders.
		 *
		 * This could have used `[].concat.apply([], this.data)` to flatten the array but Google
		 * Chrome has a maximum number of arguments so the concatenations are chunked to avoid that limit.
		 *
		 * @param {number} type Either `gl.STATIC_DRAW` or `gl.DYNAMIC_DRAW`. Defaults to `gl.STATIC_DRAW`
		 */
		compile: function(type) {
			assert('undefined' == typeof type || gl.STATIC_DRAW == type || gl.DYNAMIC_DRAW == type, 'gl.STATIC_DRAW == type || gl.DYNAMIC_DRAW == type')
			this.buffer = this.buffer || gl.createBuffer();
			var buffer
			if (this.data.length == 0) {
				console.warn("empty buffer " + this.name)
				//console.trace()
			}
			if (this.data.length == 0 || this.data[0] instanceof V3) {
				buffer = V3.flattenV3Array(this.data) // asserts that all elements are V3s
				this.spacing = 3
				this.count = this.data.length
			} else {
				assert(Array != this.data[0].constructor, this.name + this.data[0])
				var data = [];
				for (var i = 0, chunk = 10000; i < this.data.length; i += chunk) {
					data = Array.prototype.concat.apply(data, this.data.slice(i, i + chunk));
				}
				buffer = new this.type(data)
				var spacing = this.data.length ? data.length / this.data.length : 0;
				assert(spacing % 1 == 0, `buffer ${this.name}elements not of consistent size, average size is ` + spacing)
				assert(data.every(v => 'number' == typeof v), () => "data.every(v => 'number' == typeof v)" + data.toSource())
				if (NLA_DEBUG) { this.maxValue = data.max() }
				this.spacing = spacing
				this.count = data.length
			}
			gl.bindBuffer(this.target, this.buffer)
			gl.bufferData(this.target, buffer, type || gl.STATIC_DRAW)
			this.hasBeenCompiled = true
		}
	};

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
	 * @param {Object=} options
	 * @param {boolean=} options.coords Texture coordinates.
	 * @param {boolean=} options.normals
	 * @param {boolean=} options.colors
	 * @param {boolean=} options.lines
	 * @param {boolean=} options.triangles
	 * @constructor
	 * @alias GL.Mesh
	 * @property {boolean} hasBeenCompiled
	 * @property {Object.<string, GL.Buffer>} vertexBuffers
	 * @property {Object.<string, GL.Buffer>} indexBuffers
	 * @property {number[]} triangles
	 * @property {number[]} lines
	 * @property {V3[]} normals
	 * @property {V3[]} vertices
	 */
	function Mesh(options) {
		options = options || {};
		this.hasBeenCompiled = false
		/** @dict */ this.vertexBuffers = {}
		/** @dict */ this.indexBuffers = {}
		this.addVertexBuffer('vertices', 'LGL_Vertex');
		if (options.coords) this.addVertexBuffer('coords', 'LGL_TexCoord')
		if (options.normals) this.addVertexBuffer('normals', 'LGL_Normal')
		if (options.colors) this.addVertexBuffer('colors', 'LGL_Color')
		if (!('triangles' in options) || options.triangles) this.addIndexBuffer('TRIANGLES')
		if (options.lines) this.addIndexBuffer('LINES')
	}


	Mesh.prototype = /** @lends GL.Mesh.prototype */ {
		/**
		 * Add a new vertex buffer with a list as a property called `name` on this object and map it to
		 * the attribute called `attribute` in all shaders that draw this mesh.
		 * @param {string} name
		 * @param {string} attribute
		 */
		addVertexBuffer: function(name, attribute) {
			assert(!this.vertexBuffers[attribute])
			//assert(!this[name])
			this.hasBeenCompiled = false
			assert('string' == typeof name)
			assert('string' == typeof attribute)
			var buffer = this.vertexBuffers[attribute] = new Buffer(WebGLRenderingContext.ARRAY_BUFFER, Float32Array);
			buffer.name = name;
			this[name] = [];
		},

		/**
		 * Add a new index buffer.
		 * @param {string} name
		 */
		addIndexBuffer: function(name) {
			this.hasBeenCompiled = false
			var buffer = this.indexBuffers[name] = new Buffer(gl.ELEMENT_ARRAY_BUFFER, Uint16Array);
			buffer.name = name.toLowerCase()
			this[name.toLowerCase()] = []
		},


		// /** @returns {number[]} */ get triangles() { return this.indexBuffers['TRIANGLES'].data },
		// /** @returns {number[]} */ get lines() { return this.indexBuffers['LINES'].data },
		// /** @returns {V3[]} */ get vertices() { return this.vertexBuffers['LGL_Vertex'].data },
		// /** @returns {V3[]} */ get normals() { return this.vertexBuffers['LGL_Normal'].data },

		/**
		 * Upload all attached buffers to the GPU in preparation for rendering. This doesn't need to be called every
		 * frame, only needs to be done when the data changes.
		 *
		 * Sets `this.hasBeenCompiled` to true.
		 */
		compile: function() {
			// figure out shortest vertex buffer to make sure indexBuffers are in bounds
			let minVertexBufferLength = Infinity, minBufferName
			Object.getOwnPropertyNames(this.vertexBuffers).forEach(attribute => {
				let buffer = this.vertexBuffers[attribute]
				buffer.data = this[buffer.name]
				buffer.compile()
				if (this[buffer.name].length < minVertexBufferLength) {
					minBufferName = attribute
					minVertexBufferLength = this[buffer.name].length
				}
			})

			for (var name in this.indexBuffers) {
				let buffer = this.indexBuffers[name];
				buffer.data = this[buffer.name]
				buffer.compile();
				// if (NLA_DEBUG && buffer.maxValue >= minVertexBufferLength) {
				// 	throw new Error(`max index value for buffer ${name}
				// 	is too large ${buffer.maxValue} min Vbuffer size: ${minVertexBufferLength} ${minBufferName}`)
				// }
			}
			this.hasBeenCompiled = true
		},

		/**
		 * Transform all vertices by `matrix` and all normals by the inverse transpose of `matrix`.
		 *
		 * Index buffer data is referenced.
		 *
		 * @param {M4} m4
		 * @returns {GL.Mesh}
		 */
		transform: function(m4) {
			var mesh = new _GL.Mesh({vertices: this.vertices, normals: this.normals})
			mesh.vertices = m4.transformedPoints(this.vertices)
			if (this.normals) {
				let invTrans = m4.as3x3().inversed().transposed().normalized();
				mesh.normals = this.normals.map(n => invTrans.transformVector(n))
			}
			for (let name in this.indexBuffers) {
				mesh.addIndexBuffer(name)
				mesh[name.toLowerCase()] = this[name.toLowerCase()]
			}
			mesh.compile()
			return mesh
		},

		/**
		 * Computes a new normal for each vertex from the average normal of the neighboring triangles. This means
		 * adjacent triangles must share vertices for the resulting normals to be smooth.
		 * @returns {GL.Mesh}
		 */
		computeNormals: function() {
			if (!this.normals) this.addVertexBuffer('normals', 'LGL_Normal')
			this.vertexBuffers['LGL_Normal'].data = NLA.arrayFromFunction(this.vertices.length, i => V3.ZERO)

			let {triangles, vertices, normals} = this
			for (let i = 0; i < triangles.length; i++) {
				let t = triangles[i]
				let a = vertices[t[0]]
				let b = vertices[t[1]]
				let c = vertices[t[2]]
				let normal = b.minus(a).cross(c.minus(a)).normalized();
				normals[t[0]] = normals[t[0]].plus(normal);
				normals[t[1]] = normals[t[1]].plus(normal);
				normals[t[2]] = normals[t[2]].plus(normal);
			}
			for (let i = 0; i < vertices.length; i++) {
				normals[i] = normals[i].normalized()
			}
			return this
		},

		/**
		 * Populate the `lines` index buffer from the `triangles` index buffer.
		 * @returns {GL.Mesh}
		 */
		computeWireframe: function() {
			var canonEdges = new Set()
			function canonEdge(i0, i1) {
				var iMin = min(i0, i1), iMax = max(i0, i1)
				return (iMin << 16) | iMax
			}
			// function uncanonEdge(key) {
			// 	return [key >> 16, key & 0xffff]
			// }
			for (var i = 0; i < this.triangles.length; i++) {
				var t = this.triangles[i]
				canonEdges.add(canonEdge(t[0], t[1]))
				canonEdges.add(canonEdge(t[1], t[2]))
				canonEdges.add(canonEdge(t[2], t[0]))
			}
			if (!this.lines) this.addIndexBuffer('LINES');
			canonEdges.forEach(val => this.lines.push(val >> 16, val & 0xffff))
			this.hasBeenCompiled = false
			return this
		},
		computeWireframeFromFlatTriangles: function() {
			var canonEdges = new Set()
			function canonEdge(i0, i1) {
				var iMin = min(i0, i1), iMax = max(i0, i1)
				return (iMin << 16) | iMax
			}
			// function uncanonEdge(key) {
			// 	return [key >> 16, key & 0xffff]
			// }
			var t = this.triangles
			for (var i = 0; i < this.triangles.length; i += 3) {
				canonEdges.add(canonEdge(t[i + 0], t[i + 1]))
				canonEdges.add(canonEdge(t[i + 1], t[i + 2]))
				canonEdges.add(canonEdge(t[i + 2], t[i + 0]))
			}
			if (!this.lines) this.addIndexBuffer('LINES');
			//this.lines = new Array(canonEdges.size)
			canonEdges.forEach(val => this.lines.push(val >> 16, val & 0xffff))
			this.hasBeenCompiled = false
			return this
		},
		computeWireframeFromFlatTrianglesClosedMesh: function() {
			if (!this.lines) this.addIndexBuffer('LINES');
			let tris = this.triangles
			let lines = this.lines
			for (let i = 0; i < this.triangles.length; i += 3) {
				if (tris[i + 0] < tris[i + 1]) lines.push(tris[i + 0], tris[i + 1])
				if (tris[i + 1] < tris[i + 2]) lines.push(tris[i + 1], tris[i + 2])
				if (tris[i + 2] < tris[i + 0]) lines.push(tris[i + 2], tris[i + 0])
			}
			this.hasBeenCompiled = false
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
			this.hasBeenCompiled = false
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
	/**
	 * Generates a square mesh in the XY plane.
	 * Texture coordinates (buffer "coords") are set to go from 0 to 1 in either direction.
	 *
	 *
	 * @param {Object=} options
	 * @param {number=} options.detail Defaults to 1
	 * @param {number=} options.detailX Defaults to options.detail. Number of subdivisions in X direction.
	 * @param {number=} options.detailY Defaults to options.detail. Number of subdivisions in Y direction.j
	 * @param {number=} options.width defaults to 1
	 * @param {number=} options.height defaults to 1
	 * @param {number=} options.startX defaults to 0
	 * @param {number=} options.startY defaults to 0
	 *
	 * @returns {GL.Mesh}
	 * @alias GL.Mesh.plane
	 */
	Mesh.plane = function(options) {
		options = options || {}
		let detailX = options.detailX || options.detail || 1
		let detailY = options.detailY || options.detail || 1
		let startX = options.startX || 0
		let startY = options.startY || 0
		let width = options.width || 1
		let height = options.height || 1
		let mesh = new Mesh({lines: true, normals: true, coords: true, triangles: true})

		for (let j = 0; j <= detailY; j++) {
			let t = j / detailY
			for (let i = 0; i <= detailX; i++) {
				let s = i / detailX
				mesh.vertices.push(V3.create(startX + s * width, startY + t * height, 0))
				mesh.coords.push(s, t)
				mesh.normals.push(V3.Z)
				if (i < detailX && j < detailY) {
					let offset = i + j * (detailX + 1)
					mesh.triangles.push(
						offset, offset + 1, offset + detailX + 1,
						offset + detailX + 1, offset + 1, offset + detailX + 2)
				}
			}
		}

		for (let i = 0; i < detailX; i++) {
			mesh.lines.push(i, i + 1)
			mesh.lines.push((detailX + 1) * detailY + i, (detailX + 1) * detailY + i + 1)
		}
		for (let j = 0; j < detailY; j++) {
			mesh.lines.push(detailX * j, detailX * (j + 1) + 1)
			mesh.lines.push(detailX * (j + 1), detailX * (j + 2) + 1)
		}

		mesh.compile()
		return mesh
	};

	// unique corners of a unit cube. Used by Mesh.cube to generate a cube mesh.
	const UNIT_CUBE_CORNERS = [
		V3.ZERO,
		V3(0, 0, 1),
		V3(0, 1, 0),
		V3(0, 1, 1),

		V3(1, 0, 0),
		V3(1, 0, 1),
		V3(1, 1, 0),
		V3.ONES
	]

	/**
	 * Generates a unit cube (1x1x1) starting at the origin and extending into the (+ + +) octant.
	 * I.e. box from V3.ZERO to V3(1,1,1)
	 * Creates line, triangle, vertex and normal buffers.
	 *
	 * @returns {GL.Mesh}
	 */
	Mesh.cube = function() {
		var mesh = new Mesh({lines: true, triangles: true, normals: true})

		// basically indexes for faces of the cube. vertices each need to be added 3 times,
		// as they have different normals depending on the face being rendered
		let VERTEX_CORNERS = [
			0, 1, 2, 3, // X = 0
			4, 5, 6, 7, // X = 1

			0, 4, 1, 5, // Y = 0
			2, 6, 3, 7, // Y = 1

			2, 6, 0, 4, // Z = 0
			3, 7, 1, 5, // Z = 1
		]
		mesh.vertices = VERTEX_CORNERS.map(i => UNIT_CUBE_CORNERS[i])
		mesh.normals = [V3.X.negated(), V3.X, V3.Y.negated(), V3.Y, V3.Z.negated(), V3.Z].map(v => [v, v, v, v]).concatenated()
		for (let i = 0; i < 6 * 4; i += 4) {
			pushQuad(/** @type {number[]} */ mesh.triangles, 0 != i % 8,
				VERTEX_CORNERS[i], VERTEX_CORNERS[i+1], VERTEX_CORNERS[i+2], VERTEX_CORNERS[i+3])
		}
		// indexes of lines relative to UNIT_CUBE_CORNERS. Mapped to VERTEX_CORNERS.indexOf
		// so they make sense in the context of the mesh
		mesh.lines = [
			0, 1,
			0, 2,
			1, 3,
			2, 3,

			0,4,
			1,5,
			2,6,
			3,7,

			4,5,
			4,6,
			5,7,
			6,7,
		].map(i => VERTEX_CORNERS.indexOf(i))

		mesh.compile()
		return mesh
	};

	Mesh.dodecahedron = function () {
		return Mesh.sphere(0)
	}

	/**
	 * Returns a sphere mesh with radius 1 created by subdividing the faces of a dodecahedron (20-sided) recursively
	 * The sphere is positioned at the origin
	 * @param {number=} subdivisions
	 *      How many recursive divisions to do. A subdivision divides a triangle into 4,
	 *      so the total number of triangles is 20 * 4^subdivisions
	 * @returns {GL.Mesh}
	 *      Contains vertex and normal buffers and index buffers for triangles and lines
	 */
	Mesh.sphere = function (subdivisions) {
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
				triangles.push(ia, ib, ic)
				if (ia < ib) lines.push(ia, ib)
				if (ib < ic) lines.push(ib, ic)
				if (ic < ia) lines.push(ic, ia)
			} else {
				// subdivide the triangle abc into 4 by adding a vertex (with the correct distance from the origin)
				// between each segment ab, bc and cd, then calling the function recursively
				var abMid1 = a.plus(b).toLength(1), bcMid1 = b.plus(c).toLength(1), caMid1 = c.plus(a).toLength(1)
				// indexes of new vertices:
				var iabm = vertices.length, ibcm = iabm + 1, icam = iabm + 2
				vertices.push(abMid1, bcMid1, caMid1)
				tesselateRecursively(abMid1, bcMid1, caMid1, res - 1, vertices, triangles, iabm, ibcm, icam, lines)
				tesselateRecursively(a, abMid1, caMid1, res - 1, vertices, triangles, ia, iabm, icam, lines)
				tesselateRecursively(b, bcMid1, abMid1, res - 1, vertices, triangles, ib, ibcm, iabm, lines)
				tesselateRecursively(c, caMid1, bcMid1, res - 1, vertices, triangles, ic, icam, ibcm, lines)
			}
		}

		var mesh = new _GL.Mesh({normals: true, colors: false, lines: true});
		mesh.vertices.pushAll(vertices)
		subdivisions = undefined == subdivisions ? 4 : subdivisions
		for (var i = 0; i < 20; i++) {
			var [ia, ic, ib] = triangles.slice(i * 3, i * 3 + 3)
			tesselateRecursively(vertices[ia], vertices[ic], vertices[ib], subdivisions, mesh.vertices, mesh.triangles, ia, ic, ib, mesh.lines)
		}

		mesh.normals = mesh.vertices
		mesh.compile()
		console.log('mesh.lines', mesh.lines, mesh.indexBuffers)
		return mesh

	}

	/**
	 * Creates a mesh from the JSON generated by the `convert/convert.py` script.
	 * Example usage:
	 *      var data = {
	 *          vertices: [[0, 0, 0], [1, 0, 0], [0, 1, 0]],
	 *          triangles: [[0, 1, 2]]
	 *      }
	 *      var mesh = GL.Mesh.load(data)
	 *
	 * @param {Object} json
	 * @param {Object} options
	 * @returns {GL.Mesh}
	 */
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

	/**
	 * Provides a convenient wrapper for WebGL shaders. A few uniforms and attributes,
	 * prefixed with `gl_`, are automatically added to all shader sources to make
	 * simple shaders easier to write.
	 * Headers for the following variables are automatically prepended to the passed source. The correct variables are
	 * also automatically passed to the shader when drawing.
	 *
	 * For vertex and fragment shaders:
	 uniform mat3 LGL_NormalMatrix;
	 uniform mat4 LGL_ModelViewMatrix;
	 uniform mat4 LGL_ProjectionMatrix;
	 uniform mat4 LGL_ModelViewProjectionMatrix;
	 uniform mat4 LGL_ModelViewMatrixInverse;
	 uniform mat4 LGL_ProjectionMatrixInverse;
	 uniform mat4 LGL_ModelViewProjectionMatrixInverse;

	 *
	 *
	 * Example usage:
	 *
	 *  const shader = new GL.Shader(
	 *      `void main() { gl_Position = LGL_ModelViewProjectionMatrix * LGL_Vertex; }`,
	 *      `uniform vec4 color; void main() { gl_FragColor = color; }`)
	 *
	 *  shader.uniforms({ color: [1, 0, 0, 1] }).draw(mesh)
	 *
	 * Compiles a shader program using the provided vertex and fragment shaders.
	 *
	 * @param vertexSource
	 * @param fragmentSource
	 * @alias GL.Shader
	 * @constructor
	 */
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
	uniform mat3 LGL_NormalMatrix;
	uniform mat4 LGL_ModelViewMatrix;
	uniform mat4 LGL_ProjectionMatrix;
	uniform mat4 LGL_ModelViewProjectionMatrix;
	uniform mat4 LGL_ModelViewMatrixInverse;
	uniform mat4 LGL_ProjectionMatrixInverse;
	uniform mat4 LGL_ModelViewProjectionMatrixInverse;
`
		var vertexHeader = header + `
	attribute vec4 LGL_Vertex;
	attribute vec4 LGL_TexCoord;
	attribute vec3 LGL_Normal;
	attribute vec4 LGL_Color;
`
		var fragmentHeader = `  precision highp float;` + header;

		const matrixNames = header.match(/\bLGL_\w+/g)

		// Compile and link errors are thrown as strings.
		function compileSource(type, source) {
			const shader = gl.createShader(type)
			gl.shaderSource(shader, source)
			gl.compileShader(shader)
			if (!gl.getShaderParameter(shader, gl.COMPILE_STATUS)) {
				throw new Error('compile error: ' + gl.getShaderInfoLog(shader))
			}
			return shader
		}
		this.program = gl.createProgram();
		gl.attachShader(this.program, compileSource(WebGLRenderingContext.VERTEX_SHADER, vertexHeader + vertexSource));
		gl.attachShader(this.program, compileSource(WebGLRenderingContext.FRAGMENT_SHADER, fragmentHeader + fragmentSource));
		gl.linkProgram(this.program);
		if (!gl.getProgramParameter(this.program, WebGLRenderingContext.LINK_STATUS)) {
			throw new Error('link error: ' + gl.getProgramInfoLog(this.program));
		}
		this.attributes = {};
		this.uniformLocations = {};

		// Sampler uniforms need to be uploaded using `gl.uniform1i()` instead of `gl.uniform1f()`.
		// To do this automatically, we detect and remember all uniform samplers
		this.isSampler = {}
		for (let i = gl.getProgramParameter(this.program, WebGLRenderingContext.ACTIVE_UNIFORMS); i-- > 0;) {
			let uniformInfo = gl.getActiveUniform(this.program, i)
			let ut = uniformInfo.type
			if (gl.SAMPLER_2D == ut || gl.SAMPLER_CUBE == ut) {
				this.isSampler[uniformInfo.name] = 1
			}
		}

		// Check for the use of built-in matrices that require expensive matrix
		// multiplications to compute, and record these in `activeMatrices`.
		this.activeMatrices = {};
		matrixNames.forEach(name => {
			if (gl.getUniformLocation(this.program, name)) {
				this.activeMatrices[name] = true
			}
		})

		if (NLA_DEBUG) {
			this.uniformInfo = {}
			for (var i = gl.getProgramParameter(this.program, gl.ACTIVE_UNIFORMS); i-- > 0;) {
				var info = gl.getActiveUniform(this.program, i)
				this.uniformInfo[info.name] = info
			}
		}
	}

	function isArray(obj) {
		return Array == obj.constructor || Float32Array == obj.constructor || Float64Array == obj.constructor
	}

	function isNumber(obj) {
		var str = Object.prototype.toString.call(obj);
		return str == '[object Number]' || str == '[object Boolean]';
	}

	Shader.prototype = {
		/**
		 * Set a uniform for each property of `uniforms`. The correct `gl.uniform*()` method is inferred from the value
		 * types and from the stored uniform sampler flags.
		 *
		 * @param {Object} uniforms
		 * @returns {_GL.Shader}
		 */
		uniforms: function(uniforms) {
			gl.useProgram(this.program);

			for (var name in uniforms) {
				var location = this.uniformLocations[name] || gl.getUniformLocation(this.program, name);
				assert(!!location, name + ' uniform is not used in shader')
				if (!location) continue;
				this.uniformLocations[name] = location;
				let value = uniforms[name];
				if (NLA_DEBUG) {
					var info = this.uniformInfo[name]
					assert(info.type != gl.FLOAT ||
						(1 == info.size && "number" == typeof value || isArray(value) && info.size == value.length && assertNumbers.apply(undefined, value)))
					assert(info.type != gl.INT ||
						(1 == info.size && "number" == typeof value && value % 1 == 0 ||
						isArray(value) && info.size == value.length && assertNumbers.apply(undefined, value) && value.every(x => x % 1 == 0)))
					assert(info.type != gl.FLOAT_VEC3 ||
						(1 == info.size && value instanceof V3 ||
						isArray(value) && info.size == value.length && assertVectors.apply(undefined, value)))
					assert(info.type != gl.FLOAT_VEC4 || isArray(value) && value.length == 4)
					assert(info.type != gl.FLOAT_MAT4 || value instanceof M4, () => value.toSource())
					assert(info.type != gl.FLOAT_MAT3 || value.length == 9)
				}
				if (value instanceof V3) {
					value = value.toArray()
				} else if (value instanceof M4) {
					value = value.m
				}
				if (isArray(value)) {
					switch (value.length) {
						case 1: gl.uniform1fv(location, value); break;
						case 2: gl.uniform2fv(location, value); break;
						case 3: gl.uniform3fv(location, value); break;
						case 4: gl.uniform4fv(location, value); break;
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
					(this.isSampler[name] || info.type == gl.INT ? gl.uniform1i : gl.uniform1f).call(gl, location, value);
				} else {
					throw new Error('attempted to set uniform "' + name + '" to invalid value ' + value);
				}
			}

			return this;
		},

		/**
		 * Sets all uniform matrix attributes, binds all relevant buffers, and draws the mesh geometry as indexed
		 * triangles or indexed lines. Set `mode` to `gl.LINES` (and either add indices to `lines` or call
		 * `computeWireframe()`) to draw the mesh in wireframe.
		 *
		 * @param {GL.Mesh} mesh
		 * @param {string=} mode Defaults to 'TRIANGLES'. Must be passed as string so the correct index buffer can be
		 *     automatically drawn.
		 * @param {number=} start int
		 * @param {number=} count int
		 */
		draw: function(mesh, mode, start, count) {
			assert(mesh.hasBeenCompiled, 'mesh.hasBeenCompiled')
			mode = mode || 'TRIANGLES'
			assert(_GL.DRAW_MODES.includes(mode), 'GL.DRAW_MODES.includes(mode) ' + mode)
			assert(mesh.indexBuffers[mode], `The index buffer ${mode} doesn't exist!`)
			this.drawBuffers(mesh.vertexBuffers, mesh.indexBuffers[mode], gl[mode], start, count);
		},

		/**
		 * Sets all uniform matrix attributes, binds all relevant buffers, and draws the
		 * indexed mesh geometry. The `vertexBuffers` argument is a map from attribute
		 * names to `Buffer` objects of type `gl.ARRAY_BUFFER`, `indexBuffer` is a `Buffer`
		 * object of type `gl.ELEMENT_ARRAY_BUFFER`, and `mode` is a WebGL primitive mode
		 * like `gl.TRIANGLES` or `gl.LINES`. This method automatically creates and caches
		 * vertex attribute pointers for attributes as needed.
		 *
		 * @param {Object} vertexBuffers An object which maps shader attribute names to the buffer they should be bound to
		 * @param {GL.Buffer|undefined} indexBuffer
		 * @param {number} mode
		 * @param {number=} start
		 * @param {number=} count
		 * @returns {Shader} this
		 */
		drawBuffers: function(vertexBuffers, indexBuffer, mode, start, count) {
			assert(_GL.DRAW_MODES.some(stringMode => gl[stringMode] == mode), 'GL.DRAW_MODES.some(stringMode => gl[stringMode] == mode) ' + mode)
			assertf(() => 1 <= Object.keys(vertexBuffers).length)
			Object.keys(vertexBuffers).forEach(key => assertInst(Buffer, vertexBuffers[key]))

			// Only construct up the built-in matrices that are active in the shader
			const on = this.activeMatrices
			let modelViewMatrixInverse = (on['LGL_ModelViewMatrixInverse'] || on['LGL_NormalMatrix'])
				&& gl.modelViewMatrix.inversed()
			let projectionMatrixInverse = on['LGL_ProjectionMatrixInverse'] && gl.projectionMatrix.inversed()
			let modelViewProjectionMatrix = (on['LGL_ModelViewProjectionMatrix'] || on['LGL_ModelViewProjectionMatrixInverse'])
				&& gl.projectionMatrix.times(gl.modelViewMatrix)

			let /** @dict */ uni = {} // Uniform Matrices
			on['LGL_ModelViewMatrix'] && (uni['LGL_ModelViewMatrix'] = gl.modelViewMatrix)
			on['LGL_ModelViewMatrixInverse'] && (uni['LGL_ModelViewMatrixInverse'] = modelViewMatrixInverse)
			on['LGL_ProjectionMatrix'] && (uni['LGL_ProjectionMatrix'] = gl.projectionMatrix)
			on['LGL_ProjectionMatrixInverse'] && (uni['LGL_ProjectionMatrixInverse'] = projectionMatrixInverse)
			on['LGL_ModelViewProjectionMatrix'] && (uni['LGL_ModelViewProjectionMatrix'] = modelViewProjectionMatrix)
			on['LGL_ModelViewProjectionMatrixInverse'] && (uni['LGL_ModelViewProjectionMatrixInverse'] = modelViewProjectionMatrix.inversed())
			if (on['LGL_NormalMatrix']) {
				const m = modelViewMatrixInverse.m;
				// transpose normal matrix
				uni['LGL_NormalMatrix'] = new Float32Array([m[0], m[4], m[8], m[1], m[5], m[9], m[2], m[6], m[10]])
			}
			this.uniforms(uni);

			// Create and enable attribute pointers as necessary.
			let minVertexBufferLength = Infinity
			for (let attribute in vertexBuffers) {
				var buffer = vertexBuffers[attribute];
				var location = this.attributes[attribute] || gl.getAttribLocation(this.program, attribute)
				if (location == -1 || !buffer.buffer) {
					//console.warn(`Vertex buffer ${attribute} was not bound because the attribute is not active.`)
					continue
				}
				this.attributes[attribute] = location;
				gl.bindBuffer(gl.ARRAY_BUFFER, buffer.buffer)
				gl.enableVertexAttribArray(location);
				//console.log(attribute)

				gl.vertexAttribPointer(location, buffer.spacing, gl.FLOAT, false, 0, 0);

				minVertexBufferLength = Math.min(minVertexBufferLength, buffer.count)
			}

			// Disable unused attribute pointers.
			for (let attribute in this.attributes) {
				if (!(attribute in vertexBuffers)) {
					gl.disableVertexAttribArray(this.attributes[attribute]);
				}
			}

			// Draw the geometry.
			if (minVertexBufferLength && (!indexBuffer || indexBuffer.buffer)) {
				count = count || (indexBuffer ? indexBuffer.count : minVertexBufferLength)
				assert(_GL.DRAW_MODE_CHECKS[mode](count), 'count ' + count + ' doesn\'t fulfill requirement '
					+ _GL.DRAW_MODE_CHECKS[mode].toSource() + ' for mode '+ _GL.DRAW_MODES.find(modeName => gl[modeName] == mode))

				if (indexBuffer) {
					assert((count || minVertexBufferLength) % indexBuffer.spacing == 0)
					assert((start || 0) % indexBuffer.spacing == 0)
					if (start + count > indexBuffer.count) {
						throw new Error("Buffer not long enough for passed parameters start/length/buffer length" +" "+ start +" "+ count + " "+ indexBuffer.buffer.length)
					}
					gl.bindBuffer(gl.ELEMENT_ARRAY_BUFFER, indexBuffer.buffer);
					// start parameter has to be multiple of sizeof(gl.UNSIGNED_SHORT)
					gl.drawElements(mode, count, WebGLRenderingContext.UNSIGNED_SHORT, 2 * start || 0)
				} else {
					if (start + count > minVertexBufferLength) {
						throw new Error("invalid")
					}
					gl.drawArrays(mode, start || 0, count)
				}
			}

			return this;
		}
	};
	/**
	 * These are all the draw modes usable in OpenGL ES
	 * @type {string[]}
	 */
	_GL.DRAW_MODES = ['POINTS', 'LINES', 'LINE_STRIP', 'LINE_LOOP', 'TRIANGLES', 'TRIANGLE_STRIP', 'TRIANGLE_FAN']
	_GL.SHADER_VAR_TYPES = ['FLOAT', 'FLOAT_MAT2', 'FLOAT_MAT3', 'FLOAT_MAT4', 'FLOAT_VEC2', 'FLOAT_VEC3', 'FLOAT_VEC4', 'INT', 'INT_VEC2', 'INT_VEC3', 'INT_VEC4', 'UNSIGNED_INT']
	_GL.DRAW_MODE_CHECKS = {
		POINTS: x => true,
		LINES: x => 0 == x % 2, // divisible by 2
		LINE_STRIP: x => x > 2, // need at least 2
		LINE_LOOP: x => x > 2, // more like > 3, but oh well
		TRIANGLES: x => 0 == x % 3, // divisible by 3
		TRIANGLE_STRIP: x => x > 3,
		TRIANGLE_FAN: x => x > 3 }
	_GL.DRAW_MODES.forEach(modeName => _GL.DRAW_MODE_CHECKS[WebGLRenderingContext[modeName]] = _GL.DRAW_MODE_CHECKS[modeName])

	/**
	 * Provides a simple wrapper around WebGL textures that supports render-to-texture.
	 *
	 * The arguments `width` and `height` give the size of the texture in texels.
	 * WebGL texture dimensions must be powers of two unless `filter` is set to
	 * either `gl.NEAREST` or `gl.LINEAR` and `wrap` is set to `gl.CLAMP_TO_EDGE`
	 * (which they are by default).
	 *
	 * Texture parameters can be passed in via the `options` argument.
	 * Example usage:
	 *
	 *      let tex = new GL.Texture(256, 256, {
	 *       magFilter: gl.NEAREST,
	 *       minFilter: gl.LINEAR,
	 *
	 *       wrapS: gl.REPEAT,
	 *       wrapT: gl.REPEAT,
	 *
	 *       format: gl.RGB, // Defaults to gl.RGBA
	 *       type: gl.FLOAT // Defaults to gl.UNSIGNED_BYTE
	 *     })
	 *
	 * @param {number} width
	 * @param {number} height
	 * @param {Object} options
	 *
	 * @param {number=} options.wrap Defaults to gl.CLAMP_TO_EDGE, or set wrapS and wrapT individually.
	 * @param {number=} options.wrapS
	 * @param {number=} options.wrapT
	 *
	 * @param {number=} options.filter Defaults to gl.LINEAR, or set minFilter and magFilter individually.
	 * @param {number=} options.minFilter
	 * @param {number=} options.magFilter
	 *
	 * @param {number=} options.format Defaults to gl.RGBA.
	 * @param {number=} options.type Defaults to gl.UNSIGNED_BYTE.
	 * @constructor
	 * @alias GL.Texture
	 */
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

	/**
	 * Return a new texture created from `image`, an `<img>` tag.
	 *
	 * @param {HTMLImageElement} image <img> element
	 * @param {Object} options See {@link Texture} for option descriptions.
	 * @returns {_GL.Texture}
	 */
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

	/**
	 * Returns a checkerboard texture that will switch to the correct texture when it loads.
	 *
	 * @param {string} url
	 * @param {Object} options Passed to {@link Texture.fromImage}
	 * @returns {_GL.Texture}
	 */
	Texture.fromURL = function(url, options) {
		checkerboardCanvas = checkerboardCanvas || (function() {
				var c = document.createElement('canvas').getContext('2d');
				c.canvas.width = c.canvas.height = 128;
				for (var y = 0; y < c.canvas.height; y += 16) {
					for (var x = 0; x < c.canvas.width; x += 16) {
						//noinspection JSBitwiseOperatorUsage
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


	return _GL;
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
		0, 1, 2,
		3, 2, 1
	)
	mesh.compile();
	return mesh;
}

/**
 *
 * @param {V3[]} vertices
 * @param {L3} lineAxis
 * @param {number} totalAngle in radians
 * @param {number} count
 * @param {boolean} close
 * @param {V3[]=} normals
 * @returns {GL.Mesh}
 */
GL.Mesh.rotation = function(vertices, lineAxis, totalAngle, count, close=true, normals) {
	var mesh = new GL.Mesh({normals: !!normals})
	var vc = vertices.length, vTotal = vc * count

	const rotMat = M4()
	const triangles = /** @type {number[]} */ mesh.triangles
	for (let i = 0; i < count; i++) {
		// add triangles
		var angle = totalAngle / count * i
		M4.rotationLine(lineAxis.anchor, lineAxis.dir1, angle, rotMat)
		Array.prototype.push.apply(mesh.vertices, rotMat.transformedPoints(vertices))

		normals && Array.prototype.push.apply(mesh.normals, rotMat.transformedVectors(normals))
		for (let j = 0; j < vc - 1; j++) {
			pushQuad(triangles, false,
				i * vc + j + 1, i * vc + j,
				((i + 1) * vc + j + 1) % vTotal, ((i + 1) * vc + j) % vTotal)
		}
	}

	mesh.compile()
	return mesh
}

/**
 *
 * Push two triangles:
 * c - d
 * | \ |
 * a - b
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
 * @param {number=} steps
 * @returns GL.Mesh
 */
GL.Mesh.offsetVertices = function(vertices, offset, close, normals, steps) {
	assertVectors.apply(undefined, vertices)
	assertVectors(offset)

	let mesh = new GL.Mesh({normals: !!normals})
	mesh.vertices = vertices.concat(vertices.map(v => v.plus(offset)))

	let triangles = mesh.triangles
	for (let i = 0; i < vertices.length - 1; i++) {
		pushQuad(triangles, false,
			i, i + 1,
			vertices.length + i, vertices.length + i + 1)
	}
	if (close) {
		pushQuad(triangles, false, 0, vertices.length - 1, vertices.length, vertices.length * 2 - 1)
	}
	if (normals) {
		mesh.normals = normals.concat(normals)
	}
	mesh.compile()
	return mesh
}

