#version 330

uniform mat4 mvp;

in vec3 o;
//in int BC;
out vec4 p;
out vec4 hc;

void main() {
    p = vec4(o, 1.0);
	gl_Position = mvp * vec4(o, 1.0);
    //if( BC > 0)
    //    hc = vec4(1.0, 0.0, 0.0, 1.0);
    //else
        hc = vec4(1.0);
}
