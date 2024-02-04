#version 330

in vec3 position;
in vec4 q;
in vec3 scale_color;
in vec3 scale_value;


out vData {
   vec4 q; 
   vec3 scale_color;
   vec3 scale_value; 
} geo;

void main() {
	gl_Position = vec4(position, 1.0);
	geo.q = q;
	geo.scale_color = scale_color;
	geo.scale_value =scale_value * 2;
}
