#version 300 es
uniform mat4 ts_ModelViewProjectionMatrix;
uniform float startT, endT, scale;
uniform vec4 points[16];
uniform int pointCount, degree;
uniform float knots[21];
uniform vec3 normal;

in vec4 ts_Vertex;

const int MIN_DEGREE = 1;
const int MAX_DEGREE = 6;

int tInterval(float t) {
    for (int s = MIN_DEGREE; s < 21 - 1 - MIN_DEGREE; s++) {
        if (t >= knots[s] && t <= knots[s + 1]) {
            return s;
        }
    }
}

void main() {
    // ts_Vertex.x is in [0, 1]
    float t = startT + ts_Vertex.x * (endT - startT);

    int s = tInterval(t);

    vec4 v[16];
    for (int i = 0; i < 16; i++) {
        v[i] = points[i];
    }

    vec4 pTangent4;

    for (int l = 1; l <= MAX_DEGREE; l++) {
        if (!(l <= degree)) break;
        if (l == degree) {
            pTangent4 = (v[s - 1] - v[s]) *
            float(degree)
            / (knots[s] - knots[s + 1]);
        }
        // build level l of the pyramid
        for (int i = s; i > s - degree - 1 + l; i--) {
            float alpha = (t - knots[i]) / (knots[i + degree + 1 - l] - knots[i]);

            // interpolate each component
            v[i] = (1.0 - alpha) * v[i - 1] + alpha * v[i];
        }
    }

    vec4 p4 = v[s];

    vec3 p = p4.xyz / p4.w;
    vec3 pTangent = ((pTangent4.xyz * p4.w) - (p4.xyz * pTangent4.w)) / (p4.w * p4.w);

    vec3 outDir = normalize(cross(normal, pTangent));
    vec3 correctNormal = normalize(cross(pTangent, outDir));
    vec3 p2 = p + scale * (outDir * ts_Vertex.y + correctNormal * ts_Vertex.z);
    gl_Position = ts_ModelViewProjectionMatrix * vec4(p2, 1);
}
