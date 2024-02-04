#version 330
precision lowp float;

uniform vec4 split;
in vec4 p;
in vec4 hc;

out vec4 outColor;

void main() {
    if (dot(split, p) < 0)
        discard;

    outColor = hc;
}
