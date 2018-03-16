import { ShaderType } from 'tsgl'

export const fragmentShaderLighting: ShaderType<{ color: 'FLOAT_VEC3'; camPos: 'FLOAT_VEC3' }> = `
	precision highp float;
	uniform vec4 color;
	uniform vec3 camPos;
	varying vec3 normal;
	varying vec4 vPosition;
	void main() {
		vec3 normal1 = normalize(normal);
		vec3 lightPos = vec3(1000, 2000, 4000);
		vec3 lightDir = normalize(vPosition.xyz - lightPos);
        vec3 reflectionDirection = reflect(lightDir, normal1);
        vec3 eyeDirection = normalize(camPos.xyz-vPosition.xyz);
        float uMaterialShininess = 256.0;
		float specularLightWeighting = pow(max(dot(reflectionDirection, eyeDirection), 0.0), uMaterialShininess);
		float lightIntensity = 0.6 + 0.2 * max(0.0, -dot(lightDir, normal1)) + 0.2*specularLightWeighting;
		gl_FragColor = vec4(vec3(color) * lightIntensity, 1);
	}
`
export const vertexShaderLighting: ShaderType<{ color: 'FLOAT_VEC4' }> = `
	uniform mat4 ts_ModelViewProjectionMatrix;
	uniform mat4 ts_ModelViewMatrix;
	attribute vec4 ts_Vertex;
	uniform mat3 ts_NormalMatrix;
	attribute vec3 ts_Normal;
	uniform vec4 color;
	varying vec3 normal;
	varying vec4 vPosition;
	void main() {
		gl_Position = ts_ModelViewProjectionMatrix * ts_Vertex;
        vPosition = ts_ModelViewMatrix * ts_Vertex;
		normal = normalize(ts_NormalMatrix * ts_Normal);
	}
`
export const vertexShaderWaves: ShaderType<{ color: 'FLOAT_VEC3' }> = `
	uniform mat4 ts_ModelViewProjectionMatrix;
	uniform mat4 ts_ModelViewMatrix;
	attribute vec4 ts_Vertex;
	uniform mat3 ts_NormalMatrix;
	attribute vec3 ts_Normal;
	uniform vec4 color;
	varying vec3 normal;
	varying vec4 vPosition;
	void main() {
		normal = normalize(ts_NormalMatrix * ts_Normal);
		float offset = mod  (((ts_Vertex.x + ts_Vertex.y + ts_Vertex.z) * 31.0), 20.0) - 10.0;
		vec4 modPos = ts_Vertex + vec4(normal * offset, 0);
		gl_Position = ts_ModelViewProjectionMatrix * modPos;
        vPosition = ts_ModelViewMatrix * modPos;
	}
`
export const vertexShader: ShaderType<{}> = `
	varying vec4 pos;
	void main() {
		pos = vec4(position,1.0);
		gl_Position = projectionMatrix *
			modelViewMatrix *
			vec4(position,1.0);
	}
`
export const fragmentShader: ShaderType<{ color: 'FLOAT_VEC3' }> = `
	precision highp float;
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
export const vertexShaderBasic: ShaderType<{}> = `
	uniform mat4 ts_ModelViewProjectionMatrix;
	attribute vec4 ts_Vertex;
	void main() {
		gl_Position = ts_ModelViewProjectionMatrix * ts_Vertex;
	}
`
export const vertexShaderColor: ShaderType<{}> = `
	uniform mat4 ts_ModelViewProjectionMatrix;
	attribute vec4 ts_Vertex;
	attribute vec4 color;
	varying vec4 fragColor;
	void main() {
		gl_Position = ts_ModelViewProjectionMatrix * ts_Vertex;
		fragColor = color;
	}
`
export const vertexShaderArc: ShaderType<{
	color: 'FLOAT_VEC4'
	step: 'FLOAT'
	offset: 'FLOAT'
	radius: 'FLOAT'
	width: 'FLOAT'
}> = `
	uniform mat4 ts_ModelViewProjectionMatrix;
	attribute vec4 ts_Vertex;
	uniform float step, offset;
	uniform float radius, width;
	void main() {
		float r = radius;
		float t = offset + ts_Vertex.x * step;
		float pRadius = r - ts_Vertex.y * width;
		vec4 p = vec4(pRadius * cos(t), pRadius * sin(t), 0, 1);
		gl_Position = ts_ModelViewProjectionMatrix * p;
}
`
export const vertexShaderConic3d: ShaderType<{
	startT: 'FLOAT'
	endT: 'FLOAT'
	scale: 'FLOAT'
	center: 'FLOAT_VEC3'
	f1: 'FLOAT_VEC3'
	f2: 'FLOAT_VEC3'
	mode: 'INT'
}> = `
	uniform mat4 ts_ModelViewProjectionMatrix;
	attribute vec4 ts_Vertex;
	uniform float startT, endT, scale;
	uniform vec3 center, f1, f2;
	uniform int mode;
	float sinh(float x) { return (exp(x) - exp(-x)) / 2.0; }
	float cosh(float x) { return (exp(x) + exp(-x)) / 2.0; }
	void main() {
		float t = startT + ts_Vertex.x * (endT - startT);

		vec3 normal = normalize(cross(f1, f2));

		vec3 p, tangent;
		if (0 == mode) { // ellipse
			p = center + f1 * cos(t) + f2 * sin(t);
			tangent = f1 * -sin(t) + f2 * cos(t);
		}
		if (1 == mode) { // parabola
			p = center + f1 * t + f2 * t * t;
			tangent = f1 + f2 * t;
		}
		if (2 == mode) { // hyperbola
			p = center + f1 * cosh(t) + f2 * sinh(t);
			tangent = f1 * sinh(t) + f2 * cosh(t);
		}
		vec3 outDir = normalize(cross(normal, tangent));
		vec3 p2 = p + scale * (outDir * ts_Vertex.y + normal * ts_Vertex.z);
		gl_Position = ts_ModelViewProjectionMatrix * vec4(p2, 1);
	}
`
export const vertexShaderBezier: ShaderType<{
	width: 'FLOAT'
	startT: 'FLOAT'
	endT: 'FLOAT'
	p0: 'FLOAT_VEC3'
	p1: 'FLOAT_VEC3'
	p2: 'FLOAT_VEC3'
	p3: 'FLOAT_VEC3'
}> = `
    // calculates a bezier curve using ts_Vertex.x as the (t) parameter of the curve
	uniform mat4 ts_ModelViewProjectionMatrix;
	attribute vec4 ts_Vertex;
	uniform float width, startT, endT;
	uniform vec3 p0, p1, p2, p3;
	void main() {
		// ts_Vertex.y is in [0, 1]
		float t = startT + ts_Vertex.x * (endT - startT), s = 1.0 - t;
		float c0 = s * s * s, c1 = 3.0 * s * s * t, c2 = 3.0 * s * t * t, c3 = t * t * t;
		vec3 pPos = p0 * c0 + p1 * c1 + p2 * c2 + p3 * c3;
		float c01 = 3.0 * s * s, c12 = 6.0 * s * t, c23 = 3.0 * t * t;
		vec3 pTangent = (p1 - p0) * c01 + (p2 - p1) * c12 + (p3 - p2) * c23;
		vec3 pNormal = normalize(vec3(pTangent.y, -pTangent.x, 0));
		vec4 p = vec4(pPos - ts_Vertex.y * width * pNormal, 1);
		gl_Position = ts_ModelViewProjectionMatrix * p;
	}
`
export const vertexShaderBezier3d: ShaderType<{
	scale: 'FLOAT'
	startT: 'FLOAT'
	endT: 'FLOAT'
	p0: 'FLOAT_VEC3'
	p1: 'FLOAT_VEC3'
	p2: 'FLOAT_VEC3'
	p3: 'FLOAT_VEC3'
	normal: 'FLOAT_VEC3'
}> = `
    // calculates a bezier curve using ts_Vertex.x as the (t) parameter of the curve
	uniform float scale, startT, endT;
	uniform vec3 ps[4];
	uniform vec3 p0, p1, p2, p3, normal;
	uniform mat4 ts_ModelViewProjectionMatrix;
	attribute vec4 ts_Vertex;
	void main() {
		// ts_Vertex.y is in [0, 1]
		vec3 p5 = ps[0];
		float t = startT + ts_Vertex.x * (endT - startT), s = 1.0 - t;
		float c0 = s * s * s, c1 = 3.0 * s * s * t, c2 = 3.0 * s * t * t, c3 = t * t * t;
		vec3 p = p0 * c0 + p1 * c1 + p2 * c2 + p3 * c3;
		float c01 = 3.0 * s * s, c12 = 6.0 * s * t, c23 = 3.0 * t * t;
		vec3 pTangent = (p1 - p0) * c01 + (p2 - p1) * c12 + (p3 - p2) * c23;
		vec3 outDir = normalize(cross(normal, pTangent));
		vec3 correctNormal = normalize(cross(pTangent, outDir));
		vec3 p2 = p + scale * (outDir * ts_Vertex.y + correctNormal * ts_Vertex.z);
		gl_Position = ts_ModelViewProjectionMatrix * vec4(p2, 1);
	}
`
export const vertexShaderGeneric: ShaderType<{ scale: 'FLOAT' }> = `
	uniform float scale;
	uniform mat4 ts_ModelViewProjectionMatrix;
	attribute vec4 ts_Vertex;
	uniform mat3 ts_NormalMatrix;
	attribute vec3 ts_Normal;
	void main() {
		vec3 normal = normalize(ts_NormalMatrix * ts_Normal);
		vec4 vertexPos = ts_Vertex + vec4(normal * scale, 0);
		gl_Position = ts_ModelViewProjectionMatrix * vertexPos;
	}
`
export const vertexShaderRing: ShaderType<{ step: 'FLOAT'; innerRadius: 'FLOAT'; outerRadius: 'FLOAT' }> = `
	#define M_PI 3.1415926535897932384626433832795
	uniform float step;
	uniform float innerRadius, outerRadius;
	attribute float index;
	uniform mat4 ts_ModelViewProjectionMatrix;
	attribute vec4 ts_Vertex;
	void main() {
		gl_Position = ts_ModelViewProjectionMatrix * vec4(index, index, index, 1);
		float id = atan(ts_Vertex.x, ts_Vertex.y) / M_PI  * 32.0;
		float radius = mod(id, 2.0) < 1.0 ? outerRadius : innerRadius;
		gl_Position = ts_ModelViewProjectionMatrix * vec4(radius * cos(index * step), radius * sin(index * step), 0, 1);
	}
`
export const fragmentShaderColor: ShaderType<{ color: 'FLOAT_VEC4' }> = `
	precision highp float;
	uniform vec4 color;
	void main() {
		gl_FragColor = color;
	}
`
export const fragmentShaderVaryingColor: ShaderType<{}> = `
	precision highp float;
	varying vec4 fragColor;
	void main() {
		gl_FragColor = fragColor;
	}
`
export const fragmentShaderColorHighlight: ShaderType<{}> = `
	precision highp float;
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
export const vertexShaderTexture: ShaderType<{}> = `
	varying vec2 texturePos;
	attribute vec4 ts_Vertex;
	uniform mat4 ts_ModelViewProjectionMatrix;
	void main() {
		texturePos = ts_Vertex.xy;
		gl_Position = ts_ModelViewProjectionMatrix * ts_Vertex;
	}
`
export const fragmentShaderTextureColor: ShaderType<{}> = `
	precision highp float;
	varying vec2 texturePos;
	uniform vec4 color;
	uniform sampler2D texture;
	void main() {
		gl_FragColor = texture2D(texture, texturePos) * color;
	}
`
export const fragmentShaderTexture: ShaderType<{}> = `
	precision highp float;
	varying vec2 texturePos;
	uniform sampler2D texture;
	void main() {
		gl_FragColor = texture2D(texture, texturePos);
	}
`
