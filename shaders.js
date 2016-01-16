/**
 * Created by aval on 15/01/2016.
 */
var fragmentShaderLighting= `
	uniform vec4 color;
	varying vec3 normal;
	void main() {
		vec3 lightDir = normalize(vec3(-1, -2, -3));
		float lightIntensity = 0.2 + 0.8 * max(0.0, -dot(lightDir, normalize(normal)));
		gl_FragColor = vec4(vec3(color) * lightIntensity, 1);
	}
`
var vertexShaderLighting = `
	uniform vec4 color;
	varying vec3 normal;
	void main() {
		gl_Position = gl_ModelViewProjectionMatrix * gl_Vertex;
		normal = gl_NormalMatrix * gl_Normal;
	}
`
var vertexShader = `
	varying vec4 pos;
	void main() {
		pos = vec4(position,1.0);
		gl_Position = projectionMatrix *
			modelViewMatrix *
			vec4(position,1.0);
	}
`
var fragmentShader= `
	uniform vec3 color;
	varying vec4 pos;
	void main() {
		float distance = pos.x * pos.x + pos.y * pos.y;
		if (distance <= 0.98) {
			gl_FragColor = vec4(color, 1.0);
		} else if (distance <= 1.0) {
			gl_FragColor = vec4(color, 0.5);
		} else {
			gl_FragColor = vec4(0.0, 0.0, 0.0, 0.0);
		}
	}
	/*
	 precision mediump float;

	 varying vec4 pos;


	 void main() {
	 float inside = pos.r * pos.r + pos.g * pos.g;
	 if (inside <= 1) {
	 gl_FragColor = vec4(0.0, 1.0, 0.0, 1.0);
	 } else {
	 gl_FragColor = vec4(0.0, 0.0, 0.0, 0.0);
	 }
	 }
	 */
`
var vertexShaderBasic= `
	void main() {
		gl_Position = gl_ModelViewProjectionMatrix * gl_Vertex;
	}
`
var vertexShaderArc= `
	uniform float step, offset;
	uniform float innerRadius, outerRadius;
	void main() {
		float radius = gl_Vertex.x == 1.0 ? outerRadius : innerRadius, angle = offset + gl_Vertex.y * step;
		vec4 p = vec4(radius * cos(angle), radius * sin(angle), 0, 1);
		gl_Position = gl_ModelViewProjectionMatrix * p;
}
`
var vertexShaderRing= `
	#define M_PI 3.1415926535897932384626433832795
	uniform float step;
	uniform float innerRadius, outerRadius;
	attribute float index;
	void main() {
		gl_Position = gl_ModelViewProjectionMatrix * vec4(index, index, index, 1);
		float id = atan(gl_Vertex.x, gl_Vertex.y) / M_PI  * 32.0;
		float radius = mod(id, 2.0) < 1.0 ? outerRadius : innerRadius;
		gl_Position = gl_ModelViewProjectionMatrix * vec4(radius * cos(index * step), radius * sin(index * step), 0, 1);
	}
`
var fragmentShaderColor= `
	uniform vec4 color;
	void main() {
		gl_FragColor = color;
	}
`
var fragmentShaderColorHighlight= `
	uniform vec4 color;
	void main() {
		float diagonal = (gl_FragCoord.x + 2.0 * gl_FragCoord.y);
		if (mod(diagonal, 50.0) > 40.0) { // mod(diagonal, 2.0) > 1.0
			discard;
			//gl_FragColor = color + vec4(0.2,0.2,0.2,0);
		} else {
			gl_FragColor = color - vec4(0.2,0.2,0.2,0);
		}
	}
`
var vertexShaderTextureColor= `
	varying vec2 texturePos;
	void main() {
		texturePos = gl_Vertex.xy;
		gl_Position = gl_ModelViewProjectionMatrix * gl_Vertex;
	}
`
var fragmentShaderTextureColor= `
	varying vec2 texturePos;
	uniform vec4 color;
	uniform sampler2D texture;
	void main() {
		gl_FragColor = texture2D(texture, texturePos) * color;
	}
`